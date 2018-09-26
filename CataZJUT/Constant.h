#ifndef CONSTANTS_H
#define CONSTANTS_H


namespace CATAZJUT{

namespace constants {

/// Pi (\f$\pi\f$). The ratio of the circumference of a circle
/// to its diameter.
const double Pi = 3.14159265358979323846264338327;

/// Tau (\f$\tau\f$). Twice the number pi (\f$2\pi\f$). The ratio
/// of the circumference of a circle to its radius.
const double Tau = 2.0 * Pi;

/// Conversion factor from degrees to radians.
///
/// Equal to \f$\frac{\pi}{180}\f$.
const double DegreesToRadians = Pi / 180.0;

/// Conversion factor from radians to degrees.
///
/// Equal to \f$\frac{180}{\pi}\f$.
const double RadiansToDegrees = 180.0 / Pi;

/// Conversion factor from calories to joules.
const double CaloriesToJoules = 4.184;

/// Conversion factor from joules to calories.
const double JoulesToCalories = 1.0 / CaloriesToJoules;

/// Avogadro constant (\f$N_{A}\f$).
const double AvogadroConstant = 6.02214179e23;

/// Boltzmann constant (\f$k_{B}\f$)
const double BoltzmannConstant = 1.3806503e-23;

/// Gas constant.
///
/// Equal to \f$N_{A} k_{B}\f$.
const double GasConstant = AvogadroConstant * BoltzmannConstant;

/// Planck constant (\f$h\f$).
const double PlanckConstant = 6.62606896e-34;

/// Reduced Planck constant (\f$\hbar\f$).
const double ReducedPlanckConstant = PlanckConstant / (2.0 * Pi);

/// Mass of a proton (\f$m_{p}\f$).
const double ProtonMass = 1.672621637e-27;

/// Mass of an electron (\f$m_{e}\f$).
const double ElectronMass = 9.10938215e-31;

/// Elementary charge (\f$e\f$).
const double ElementaryCharge = 1.602176487e-19;

/// Charge of a proton (\f$+e\f$).
const double ProtonCharge = ElementaryCharge;

/// Charge of an electron (\f$-e\f$).
const double ElectronCharge = -ElementaryCharge;

/// Faraday constant (\f$F\f$).
///
/// Equal to \f$e N_{A}\f$.
const double FaradayConstant = AvogadroConstant * ElementaryCharge;

/// Speed of light (\f$c\f$) in meters per second (\f$\frac{m}{s}\f$).
const double SpeedOfLight = 299792458;

/// Vacuum permeability (\f$\mu_{0}\f$).
const double VacuumPermeability = 1.2566370614e-6;

/// Vacuum permittivity (\f$\epsilon_{0}\f$).
///
/// Equal to \f$\frac{1}{\mu_{0} c^{2}}\f$.
const double VacuumPermittivity = 1.0 / (VacuumPermeability * SpeedOfLight * SpeedOfLight);

/// Fine structure constant (\f$\alpha\f$).
///
/// Equal to \f$\frac{e^{2} c \epsilon_{0}}{2 h}\f$.
const double FineStructureConstant = (ElementaryCharge * ElementaryCharge * SpeedOfLight * VacuumPermittivity) / (2.0 * PlanckConstant);

/// Conversion factor between Bohr units and Angstroms.
const double BohrToAngstroms = 0.52918;

} // end constants namespace

}

#endif
