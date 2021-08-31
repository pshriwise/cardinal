#include "NekRSProblem1d.h"
#include "Moose.h"
#include "AuxiliarySystem.h"
#include "TimeStepper.h"
#include "NekInterface.h"
#include "TimedPrint.h"

#include "nekrs.hpp"
#include "nekInterface/nekInterfaceAdapter.hpp"

registerMooseObject("CardinalApp", NekRSProblem1d);

bool NekRSProblem1d::_first = true;

template<>
InputParameters
validParams<NekRSProblem1d>()
{
  InputParameters params = validParams<ExternalProblem>();
  params.addParam<bool>("minimize_transfers_in", false, "Whether to only synchronize nekRS "
    "for the direction TO_EXTERNAL_APP on multiapp synchronization steps");
  params.addParam<bool>("minimize_transfers_out", false, "Whether to only synchronize nekRS "
    "for the direction FROM_EXTERNAL_APP on multiapp synchronization steps");

  params.addParam<bool>("nondimensional", false, "Whether nekRS is solved in non-dimensional form");
  params.addParam<bool>("moving_mesh", false, "Whether we have a moving mesh problem or not");
  params.addRangeCheckedParam<Real>("U_ref", 1.0, "U_ref > 0.0", "Reference velocity value for non-dimensional solution");
  params.addRangeCheckedParam<Real>("T_ref", 0.0, "T_ref >= 0.0", "Reference temperature value for non-dimensional solution");
  params.addRangeCheckedParam<Real>("dT_ref", 1.0, "dT_ref > 0.0", "Reference temperature range value for non-dimensional solution");
  params.addRangeCheckedParam<Real>("L_ref", 1.0, "L_ref > 0.0", "Reference length scale value for non-dimensional solution");
  params.addRangeCheckedParam<Real>("rho_0", 1.0, "rho_0 > 0.0", "Density parameter value for non-dimensional solution");
  params.addRangeCheckedParam<Real>("Cp_0", 1.0, "Cp_0 > 0.0", "Heat capacity parameter value for non-dimensional solution");

  params.addParam<PostprocessorName>("min_T", "If provided, postprocessor used to limit the minimum "
    "temperature (in dimensional form) in the nekRS problem");
  params.addParam<PostprocessorName>("max_T", "If provided, postprocessor used to limit the maximum "
    "temperature (in dimensional form) in the nekRS problem");

  params.addParam<bool>("SAMtoNek_interface", false, "inlet interface from SAM to nekRS");

  return params;
}

NekRSProblem1d::NekRSProblem1d(const InputParameters &params) : ExternalProblem(params),
    _serialized_solution(NumericVector<Number>::build(_communicator).release()),
    _moving_mesh(getParam<bool>("moving_mesh")),
    _minimize_transfers_in(_moving_mesh ? true : getParam<bool>("minimize_transfers_in")),
    _minimize_transfers_out(getParam<bool>("minimize_transfers_out")),
    _nondimensional(getParam<bool>("nondimensional")),
    _U_ref(getParam<Real>("U_ref")),
    _T_ref(getParam<Real>("T_ref")),
    _dT_ref(getParam<Real>("dT_ref")),
    _L_ref(getParam<Real>("L_ref")),
    _rho_0(getParam<Real>("rho_0")),
    _Cp_0(getParam<Real>("Cp_0")),
    _SAMtoNek_interface(getParam<bool>("SAMtoNek_interface")),
    _start_time(nekrs::startTime())
{
  // the NekRSProblem1d constructor is called right after building the mesh. In order
  // to have pretty screen output without conflicting with the timed print messages,
  // print diagnostic info related to the mesh here
  _nek_mesh = dynamic_cast<NekRSMesh*>(&mesh());

  if (!_nek_mesh)
    mooseError("Mesh for a 'NekRSProblem1d' must be of type 'NekRSMesh'! In your [Mesh] "
      "block, you should have 'type = NekRSMesh'");

  _nek_mesh->printMeshInfo();

  // if the mesh is moving, then we must minimize the incoming data transfers;
  // if the user set `minimize_transfers_in = false`, print a warning that we're overriding this setting
  if (_moving_mesh && params.isParamSetByUser("minimize_transfers_in"))
  {
    auto user_setting = getParam<bool>("minimize_transfers_in");
    if (!user_setting)
      mooseWarning("Overriding 'minimize_transfers_in' to 'true' for moving mesh problems!");
  }

  // if solving in nondimensional form, make sure that the user specified _all_ of the
  // necessary scaling quantities to prevent errors from forgetting one, which would take
  // a non-scaled default otherwise
  std::vector<std::string> scales = {"U_ref", "T_ref", "dT_ref", "L_ref", "rho_0", "Cp_0"};
  std::vector<std::string> descriptions = {"velocity", "temperature", "temperature range", "length",
    "density", "heat capacity"};
  for (std::size_t n = 0; n < scales.size(); ++n)
  {
    if (_nondimensional && !params.isParamSetByUser(scales[n]))
      paramError(scales[n], "When nekRS solves in non-dimensional form, a characterstic " + descriptions[n] + " must be provided!");
    else if (!_nondimensional && params.isParamSetByUser(scales[n]))
      mooseWarning("When nekRS solves in dimensional form, " + descriptions[n] + " is unused!");
  }

  // inform nekRS of the scaling that we are using if solving in non-dimensional form
  nekrs::solution::initializeDimensionalScales(_U_ref, _T_ref, _dT_ref, _L_ref, _rho_0, _Cp_0);

  // the way the data transfers are detected depend on nekRS being a sub-application,
  // so these settings are not invalid if nekRS is the master app (though you could
  // relax this in the future by reversing the synchronization step identification
  // from the nekRS-subapp case to the nekRS-master app case - it's just not implemented yet).
  if (_app.isUltimateMaster())
    if (_minimize_transfers_in || _minimize_transfers_out)
      mooseError("The 'minimize_transfers_in' and 'minimize_transfers_out' capabilities "
        "require that nekRS is receiving and sending data to a master application, but "
        "in your case nekRS is the master application.");

  // It's too complicated to make sure that the dimensional form _also_ works when our
  // reference coordinates are different from what MOOSE is expecting, so just throw an error
  if (_nondimensional && (std::abs(_nek_mesh->scaling() - _L_ref) > 1e-6))
    paramError("L_ref", "When solving in non-dimensional form, no capability exists to allow "
      "a nondimensional solution based on reference scales that are not in the same units as the "
      "coupled MOOSE application!\n\nIf solving nekRS in nondimensional form, you must choose "
      "reference dimensional scales in the same units as expected by MOOSE, i.e. 'L_ref' "
      "must match 'scaling' in 'NekRSMesh'.");

  // boundary-specific data
  _boundary = _nek_mesh->boundary();
  _n_surface_elems = _nek_mesh->numSurfaceElems();
  _n_vertices_per_surface = _nek_mesh->numVerticesPerSurface();

  // volume-specific data
  _volume = _nek_mesh->volume();
  _n_volume_elems = _nek_mesh->numVolumeElems();
  _n_vertices_per_volume = _nek_mesh->numVerticesPerVolume();

  // generic data
  _n_elems = _nek_mesh->numElems();
  _n_vertices_per_elem = _nek_mesh->numVerticesPerElem();

  // Depending on the type of coupling, initialize various problem parameters
  if (_boundary && !_volume) // only boundary coupling
  {
    _incoming = "boundary heat flux";
    _outgoing = "boundary temperature";
    _n_points = _n_surface_elems * _n_vertices_per_surface;
    _flux_face = (double *) calloc(_n_vertices_per_surface, sizeof(double));
  }
  else if (_volume && !_boundary) // only volume coupling
  {
    _incoming = "volume power density";
    _outgoing = "volume temperature";
    _n_points = _n_volume_elems * _n_vertices_per_volume;
    _source_elem = (double*) calloc(_n_vertices_per_volume, sizeof(double));
  }
  else // both volume and boundary coupling
  {
    _incoming = "boundary heat flux and volume power density";
    _outgoing = "volume temperature";
    _n_points = _n_volume_elems * _n_vertices_per_volume;
    _flux_elem = (double *) calloc(_n_vertices_per_volume, sizeof(double));
    _source_elem = (double*) calloc(_n_vertices_per_volume, sizeof(double));
  }

  if (_moving_mesh)
  {
    _incoming += " and mesh displacement";

    if (_boundary)
    {
      mooseError("Mesh displacement not supported in boundary coupling!");
      // pending release of mesh solver in nekRS...
    }
    else if (_volume)
    {
      nekrs::save_initial_mesh();
      _displacement_x = (double *) calloc(_n_vertices_per_volume, sizeof(double));
      _displacement_y = (double *) calloc(_n_vertices_per_volume, sizeof(double));
      _displacement_z = (double *) calloc(_n_vertices_per_volume, sizeof(double));
    }
  }


  // regardless of the boundary/volume coupling, we will always exchange temperature
  _T = (double*) calloc(_n_points, sizeof(double));

  nekrs::initializeInterpolationMatrices(_nek_mesh->numQuadraturePoints1D());

  // we can save some effort for the low-order situations where the interpolation
  // matrix is the identity matrix (i.e. for which equi-spaced libMesh nodes are an
  // exact subset of the nekRS GLL points). This will happen for any first-order mesh,
  // and if a second-order mesh is used with a polynomial order of 2 in nekRS. Because
  // we pretty much always use a polynomial order greater than 2 in nekRS, let's just
  // check the first case because this will simplify our code in the nekrs::boundarySolution
  // function. If you change this line, you MUST change the innermost if/else statement
  // in nekrs::boundarySolution!
  _needs_interpolation = _nek_mesh->numQuadraturePoints1D() > 2;
}

NekRSProblem1d::~NekRSProblem1d()
{
  // write nekRS solution to output if not already written for this step
  if (!isOutputStep())
    nekrs::outfld(_timestepper->nondimensionalDT(_time));

  if (_T) free(_T);
  if (_flux_face) free(_flux_face);
  if (_source_elem) free(_source_elem);
  if (_flux_elem) free(_flux_elem);
  if (_moving_mesh) { 
   free(_displacement_x) ;
   free(_displacement_y) ;
   free(_displacement_z) ;
  }

}

void
NekRSProblem1d::initialSetup()
{
  ExternalProblem::initialSetup();

  // While we don't require nekRS to actually _solve_ for the temperature, we should
  // print a warning if there is no temperature solve. For instance, the check in
  // NekApp makes sure that we have a [TEMPERATURE] block in the nekRS input file, but we
  // might still toggle the solver off by setting 'solver = none'. Warn the user if
  // the solve is turned off because this is really only a testing feature.
  bool has_temperature_solve = nekrs::hasTemperatureSolve();
  if (!has_temperature_solve)
    mooseWarning("By setting 'solver = none' for temperature in the .par file, nekRS "
      "will not solve for temperature.\n\nThe temperature transferred to MOOSE will remain "
      "fixed at its initial condition, and the heat flux and power transferred to nekRS will be unused.");

  // For boundary-based coupling, we should check that the correct flux boundary
  // condition is set on all of nekRS's boundaries. To avoid throwing this
  // error for test cases where we have a [TEMPERATURE] block but set its solve
  // to 'none', we also check whether we're actually computing for the temperature.
  auto boundary = _nek_mesh->boundary();
  if (boundary && has_temperature_solve)
  {
    for (const auto & b : *boundary)
      if (!nekrs::mesh::isHeatFluxBoundary(b))
      {
        const std::string type = nekrs::mesh::temperatureBoundaryType(b);
        mooseError("In order to send a boundary heat flux to nekRS, you must have a flux condition "
          "for each 'boundary' set in 'NekRSMesh'!\nBoundary " + std::to_string(b) + " is of type '" +
          type + "' instead of 'fixedGradient'.");
      }
  }

  // For volume-based coupling, we should check that there is a udf function providing
  // the source for the passive scalar equations (this is the analogue of the boundary
  // condition check for boundary-based coupling). NOTE: This check is imperfect, because
  // even if there is a source kernel, we cannot tell _which_ passive scalar equation that
  // it is applied to (we have source kernels for the RANS passive scalar equations, for instance).
  if (_nek_mesh->volume())
    if (has_temperature_solve && !nekrs::hasHeatSourceKernel())
      mooseError("In order to send a heat source to nekRS, you must have an OCCA kernel "
        "for the source in the passive scalar equations!");

  if (_moving_mesh && !nekrs::hasMovingMesh())
    mooseError("In order for MOOSE to compute a mesh deformation in NekRS, you "
      "must have 'solver = user' in the [MESH] block!");

  auto executioner = _app.getExecutioner();
  _transient_executioner = dynamic_cast<Transient *>(executioner);

  // nekRS only supports transient simulations - therefore, it does not make
  // sense to use anything except a Transient-derived executioner
  if (!_transient_executioner)
    mooseError("A Transient-type executioner should be used for nekRS!");

  // If the simulation start time is not zero, the app's time must be shifted
  // relative to its master app (if any). Until this is implemented, make sure
  // a start time of zero is used.
  const auto moose_start_time = _transient_executioner->getStartTime();
  if (moose_start_time != 0.0)
    mooseError("A non-zero start time is not yet available for 'NekRSProblem1d'!");

  // To get the correct time stepping information on the MOOSE side, we also
  // must use the NekTimeStepper
  TimeStepper * stepper = _transient_executioner->getTimeStepper();
  _timestepper = dynamic_cast<NekTimeStepper *>(stepper);
  if (!_timestepper)
    mooseError("The 'NekTimeStepper' stepper must be used with 'NekRSProblem1d'!");

  // Set the reference time for use in dimensionalizing/non-dimensionalizing the time
  _timestepper->setReferenceTime(_L_ref, _U_ref);

  // Also make sure that the start time is consistent with what MOOSE wants to use.
  // If different from what nekRS internally wants to use, use the MOOSE value.
  if (std::abs(moose_start_time - _start_time) > 1e-8)
  {
    mooseWarning("The start time set on 'NekRSProblem1d': " + Moose::stringify(moose_start_time) +
      " does not match the start time set in nekRS's .par file: " + Moose::stringify(_timestepper->dimensionalDT(_start_time)) + ". "
      "This may happen if you are using a restart file in nekRS.\n\n" +
      "Setting start time for 'NekRSProblem1d' to: " + Moose::stringify(moose_start_time));
    _start_time = moose_start_time;
  }

  // Then, dimensionalize the nekRS time so that all occurrences of _dt here are
  // in dimensional form
  _timestepper->dimensionalizeDT();

  if (_boundary)
  {
    if (_SAMtoNek_interface)
    {
    _SAM_mflow_inlet_interface = &getPostprocessorValueByName("SAM_mflow_inlet_interface");
    }
  }
  if (_minimize_transfers_in)
    _transfer_in = &getPostprocessorValueByName("transfer_in");

  if (isParamValid("min_T"))
  {
    auto name = getParam<PostprocessorName>("min_T");
    _min_T = &getPostprocessorValueByName(name);
  }

  if (isParamValid("max_T"))
  {
    auto name = getParam<PostprocessorName>("max_T");
    _max_T = &getPostprocessorValueByName(name);
  }

  // nekRS calls UDF_ExecuteStep once before the time stepping begins
  nekrs::udfExecuteStep(_start_time, _t_step, false /* not an output step */);

  // save initial mesh for moving mesh problems to match deformation in exodus output files
  if (_moving_mesh)
    nekrs::outfld(_timestepper->nondimensionalDT(_time));
}

bool
NekRSProblem1d::isOutputStep() const
{
  if (_app.isUltimateMaster())
  {
    bool last_step = nekrs::lastStep(_timestepper->nondimensionalDT(_time), _t_step, 0.0 /* dummy elapsed time */);

    // if NekApp is controlled by a master application, then the last time step
    // is controlled by that master application, in which case we don't want to
    // write at what nekRS thinks is the last step (since it may or may not be
    // the actual end step), especially because we already ensure that we write on the
    // last time step from MOOSE's perspective in NekRSProblem1d's destructor.
    if (last_step)
      return true;
  }

  // this routine does not check if we are on the last step - just whether we have
  // met the requested runtime or time step interval
  return nekrs::outputStep(_timestepper->nondimensionalDT(_time), _t_step);
}

void NekRSProblem1d::externalSolve()
{
  // The _dt member of NekRSProblem1d reflects the time step that MOOSE wants NekApp to
  // take. For instance, if NekApp is controlled by a master app and subcycling is used,
  // NekApp must advance to the time interval taken by the master app. If the time step
  // that MOOSE wants nekRS to take (i.e. _dt) is smaller than we'd like nekRS to take, error.
  if (_dt < _timestepper->minDT())
    mooseError("Requested time step of " + std::to_string(_dt) + " is smaller than the minimum "
      "time step allowed in nekRS!");

  // By using the _time object on the ExternalProblem base class (which represents the
  // time that we're simulating _to_, we need to pass sometimes slightly different
  // times into the nekRS routines, which assume that the "time" passed into their
  // routines is sometimes a different interpretation.
  double step_start_time = _time - _dt;
  double step_end_time = _time;

  bool is_output_step = isOutputStep();

  // Run a nekRS time step. After the time step, this also calls UDF_ExecuteStep,
  // evaluated at (step_end_time, _t_step)
  nekrs::runStep(_timestepper->nondimensionalDT(step_start_time),
    _timestepper->nondimensionalDT(_dt), _t_step);

  // Note: here, we copy to both the nrs solution arrays and to the Nek5000 backend arrays,
  // because it is possible that users may interact using the legacy usr-file approach.
  // If we move away from the Nek5000 backend entirely, we could replace this line with
  // direct OCCA memcpy calls. But we do definitely need some type of copy here for _every_
  // time step, even if we're not technically passing data to another app, because we have
  // postprocessors that touch the `nrs` arrays that can be called in an arbitrary fashion
  // by the user.
  nek::ocopyToNek(_timestepper->nondimensionalDT(step_end_time), _t_step);

  // limit the temperature based on user settings
  bool limit_temperature = _min_T || _max_T;
  std::string msg;
  if (_min_T && !_max_T)
    msg = "Limiting nekRS temperature to above minimum temperature of " + Moose::stringify(*_min_T);
  if (_max_T && !_min_T)
    msg = "Limiting nekRS temperature to below maximum temperature of " + Moose::stringify(*_max_T);
  if (_max_T && _min_T)
    msg = "Limiting nekRS temperature to within the range [" + Moose::stringify(*_min_T) + ", " +
      Moose::stringify(*_max_T);

  if (limit_temperature)
  {
    CONTROLLED_CONSOLE_TIMED_PRINT(0.0, 1.0, msg);
    nekrs::limitTemperature(_min_T, _max_T);
  }

  if (is_output_step)
    nekrs::outfld(_timestepper->nondimensionalDT(step_end_time));

  _time += _dt;
}

bool
NekRSProblem1d::synchronizeIn()
{
  bool synchronize = true;
  static bool first = true;

  if (_minimize_transfers_in)
  {
    // For the minimized incoming synchronization to work correctly, the value
    // of the incoming postprocessor must not be zero. We only need to check this for the very
    // first time we evaluate this function. This ensures that you don't accidentally set a
    // zero value as a default in the master application's postprocessor.
    if (first && *_transfer_in == false)
      mooseError("The default value for the 'transfer_in' postprocessor received by nekRS "
        "must not be false! Make sure that the master application's "
        "postprocessor is not zero.");

    if (*_transfer_in == false)
      synchronize = false;
    else
      setPostprocessorValueByName("transfer_in", false, 0);
  }

  first = false;
  return synchronize;
}

bool
NekRSProblem1d::synchronizeOut()
{
  bool synchronize = true;

  if (_minimize_transfers_out)
  {
    if (std::abs(_time - _dt - _transient_executioner->getTargetTime()) > _transient_executioner->timestepTol())
      synchronize = false;
  }

  return synchronize;
}

void
NekRSProblem1d::sendBoundaryVelocityCorrectedToNek()
{
  _console << "Sending corrected velocity to nekRS from first boundary specified in NekRSMesh ";

  auto & solution = _aux->solution();
  auto sys_number = _aux->number();

  if (_first)
  {
    _serialized_solution->init(_aux->sys().n_dofs(), false, SERIAL);
    _first = false;
  }

  solution.localize(*_serialized_solution);

  auto & mesh = _nek_mesh->getMesh();

  const double& _sam_mflow_inlet = *_SAM_mflow_inlet_interface;

  const vector<int>& vecRef = *_boundary; 
  
  double rhoArea_inlet =  nekrs::rhoArea_Direct(vecRef[0]); //TODO, assumes first BC listed in nekRSMesh is inlet Boundary

  const double _u_inlet_corrected = _sam_mflow_inlet / rhoArea_inlet ;

  for (unsigned int e = 0; e < _n_surface_elems; e++)
  {
    nekrs::u_inlet(e, _nek_mesh->order(), _u_inlet_corrected);
  }

  _console << _u_inlet_corrected << std::endl;

  _console << "done" << std::endl;

}

void
NekRSProblem1d::getBoundaryTemperatureFromNek()
{
  CONTROLLED_CONSOLE_TIMED_PRINT(0.0, 1.0, "Extracting nekRS temperature from boundary " + Moose::stringify(*_boundary));

  // Get the temperature solution from nekRS. Note that nekRS performs a global communication
  // here such that each nekRS process has all the boundary temperature information. That is,
  // every process knows the full boundary temperature solution
  nekrs::boundarySolution(_nek_mesh->order(), _needs_interpolation, field::temperature, _T);
}

void NekRSProblem1d::syncSolutions(ExternalProblem::Direction direction)
{
  switch(direction)
  {
    case ExternalProblem::Direction::TO_EXTERNAL_APP:
    {
      if (!synchronizeIn())
      {
        _console << "Skipping " << _incoming << " transfer to nekRS, not at synchronization step" << std::endl;
        return;
      }

      if (_boundary)
      {
        if (_SAMtoNek_interface)
        {  
          sendBoundaryVelocityCorrectedToNek();
        }
      }

      nekrs::copyScratchToDevice();

      break;
    }

    case ExternalProblem::Direction::FROM_EXTERNAL_APP:
    {
      if (!synchronizeOut())
      {
        _console << "Skipping " << _outgoing << " transfer out of nekRS, not at synchronization step" << std::endl;
        return;
      }

      if (!_volume)
        getBoundaryTemperatureFromNek();

      _console << " Interpolated temperature min/max values: " <<
        minInterpolatedTemperature() << ", " << maxInterpolatedTemperature() << std::endl;

      break;
    }
    default:
      mooseError("Unhandled 'Transfer::DIRECTION' enum!");
  }
}

double
NekRSProblem1d::maxInterpolatedTemperature() const
{
  double maximum = std::numeric_limits<double>::min();

  for (int i = 0; i < _n_points; ++i)
    maximum = std::max(maximum, _T[i]);

  return maximum;
}

double
NekRSProblem1d::minInterpolatedTemperature() const
{
  double minimum = std::numeric_limits<double>::max();

  for (int i = 0; i < _n_points; ++i)
    minimum = std::min(minimum, _T[i]);

  return minimum;
}

void
NekRSProblem1d::addExternalVariables()
{
  auto var_params = _factory.getValidParams("MooseVariable");
  var_params.set<MooseEnum>("family") = "LAGRANGE";

  switch (_nek_mesh->order())
  {
    case order::first:
      var_params.set<MooseEnum>("order") = "FIRST";
      break;
    case order::second:
      var_params.set<MooseEnum>("order") = "SECOND";
      break;
    default:
      mooseError("Unhandled 'NekOrderEnum' in 'NekRSProblem1d'!");
  }

  // Because this temperature represents the reconstruction of nekRS's temperature
  // onto the NekRSMesh, we set the order to match the desired order of the mesh.
  // Note that this does _not_ imply anything about the order of the temperature
  // variable in the MOOSE app (such as BISON) coupled to nekRS. This is just the
  // variable that nekRS writes into, and then MOOSE's transfer classes can handle
  // any additional interpolations needed from 'temp' into the receiving-app's fields.
  addAuxVariable("MooseVariable", "temp", var_params);
  _temp_var = _aux->getFieldVariable<Real>(0, "temp").number();

  if (_boundary)
  {
    // Likewise, because this flux represents the reconstruction of the flux variable
    // that becomes a boundary condition in the nekRS model, we set the order to match
    // the desired order of the surface. Note that this does _not_ imply anything
    // about the order of the surface flux in the MOOSE app (such as BISON) coupled
    // to nekRS. This is just the variable that nekRS reads from - MOOSE's transfer
    // classes handle any additional interpolations needed from the flux on the
    // sending app (such as BISON) into 'avg_flux'.
    addAuxVariable("MooseVariable", "avg_flux", var_params);
    _avg_flux_var = _aux->getFieldVariable<Real>(0, "avg_flux").number();

    // add the postprocessor that receives the flux integral for normalization
    auto pp_params = _factory.getValidParams("Receiver");
    addPostprocessor("Receiver", "flux_integral", pp_params);
  }

  if (_volume)
  {
    addAuxVariable("MooseVariable", "heat_source", var_params);
    _heat_source_var = _aux->getFieldVariable<Real>(0, "heat_source").number();

    // add the postprocessor that receives the source integral for normalization
    auto pp_params = _factory.getValidParams("Receiver");
    addPostprocessor("Receiver", "source_integral", pp_params);
  }

  // add the displacement aux variables from the solid mechanics solver; these will
  // be needed regardless of whether the displacement is boundary- or volume-based
  if (_moving_mesh)
  {
    addAuxVariable("MooseVariable", "disp_x", var_params);
    _disp_x_var = _aux->getFieldVariable<Real>(0, "disp_x").number();

    addAuxVariable("MooseVariable", "disp_y", var_params);
    _disp_y_var = _aux->getFieldVariable<Real>(0, "disp_y").number();

    addAuxVariable("MooseVariable", "disp_z", var_params);
    _disp_z_var = _aux->getFieldVariable<Real>(0, "disp_z").number();
  }

  if (_minimize_transfers_in)
  {
    auto pp_params = _factory.getValidParams("Receiver");
    addPostprocessor("Receiver", "transfer_in", pp_params);
  }
}
