
// Physical Unit & Constant
//
// GEANT4\source\externals\clhep\include\CLHEP\Units  

#pragma once
  
 // SystemOfUnits.h
 //
 // Basic Unit : mm, rad, ns, eplus, MeV
 //
 // Length [L]
 //
 static const double millimeter  = 1.;                        
 static const double millimeter2 = millimeter*millimeter;
 static const double millimeter3 = millimeter*millimeter*millimeter;
 
 static const double centimeter  = 10.*millimeter;   
 static const double centimeter2 = centimeter*centimeter;
 static const double centimeter3 = centimeter*centimeter*centimeter;

 static const double meter  = 1000.*millimeter;                  
 static const double meter2 = meter*meter;
 static const double meter3 = meter*meter*meter;

 static const double kilometer = 1000.*meter;                   
 static const double kilometer2 = kilometer*kilometer;
 static const double kilometer3 = kilometer*kilometer*kilometer;

 static const double parsec = 3.0856775807e+16*meter;

 static const double micrometer = 1.e-6 *meter;             
 static const double  nanometer = 1.e-9 *meter;
 static const double  angstrom  = 1.e-10*meter;
 static const double  fermi     = 1.e-15*meter;
 
 static const double      barn = 1.e-28*meter2;
 static const double millibarn = 1.e-3 *barn;
 static const double microbarn = 1.e-6 *barn;
 static const double  nanobarn = 1.e-9 *barn;
 static const double  picobarn = 1.e-12*barn;
 
 // symbols
 static const double nm  = nanometer;                        
 static const double um  = micrometer;                        
 
 static const double mm  = millimeter;                        
 static const double mm2 = millimeter2;
 static const double mm3 = millimeter3;

 static const double cm  = centimeter;   
 static const double cm2 = centimeter2;
 static const double cm3 = centimeter3;

 static const double m  = meter;                  
 static const double m2 = meter2;
 static const double m3 = meter3;

 static const double km  = kilometer;                   
 static const double km2 = kilometer2;
 static const double km3 = kilometer3;
 
 static const double pc = parsec;
 
 //
 // Angle
 //
 static const double radian      = 1.;                  
 static const double milliradian = 1.e-3*radian;
 static const double degree = (3.14159265358979323846/180.0)*radian;
 
 static const double   steradian = 1.;
 
 // symbols
 static const double rad  = radian;
 static const double mrad = milliradian;
 static const double sr   = steradian;
 static const double deg  = degree;

 //
 // Time [T]
 //
 static const double nanosecond  = 1.;
 static const double second      = 1.e+9 *nanosecond;
 static const double millisecond = 1.e-3 *second;
 static const double microsecond = 1.e-6 *second;
 static const double  picosecond = 1.e-12*second;
 
 static const double hertz = 1./second;
 static const double kilohertz = 1.e+3*hertz;
 static const double megahertz = 1.e+6*hertz;
 
 // symbols
 static const double ns = nanosecond;
 static const double  s = second;
 static const double ms = millisecond;

 //
 // Electric charge [Q]
 //
 static const double eplus = 1. ;// positron charge
 static const double e_SI  = 1.602176487e-19;// positron charge in coulomb
 static const double coulomb = eplus/e_SI;// coulomb = 6.24150 e+18 * eplus
 
 //
 // Energy [E]
 //
 static const double megaelectronvolt = 1. ;
 static const double     electronvolt = 1.e-6*megaelectronvolt;
 static const double kiloelectronvolt = 1.e-3*megaelectronvolt;
 static const double gigaelectronvolt = 1.e+3*megaelectronvolt;
 static const double teraelectronvolt = 1.e+6*megaelectronvolt;
 static const double petaelectronvolt = 1.e+9*megaelectronvolt;

 static const double joule = electronvolt/e_SI;// joule = 6.24150 e+12 * MeV
 
 // symbols
 static const double MeV = megaelectronvolt;
 static const double  eV = electronvolt;
 static const double keV = kiloelectronvolt;
 static const double GeV = gigaelectronvolt;
 static const double TeV = teraelectronvolt;
 static const double PeV = petaelectronvolt;

// PhysicalConstants.h
//
static const double     pi  = 3.14159265358979323846;
static const double  twopi  = 2*pi;
static const double halfpi  = pi/2;
static const double     pi2 = pi*pi;

//
// c   = 299.792458 mm/ns
// c^2 = 898.7404 (mm/ns)^2 
//
static const double c   = 2.99792458e+8 * m/s;
static const double c2 = c * c;

//
// h     = 4.13566e-12 MeV*ns
// hbar  = 6.58212e-13 MeV*ns
// hbarc = 197.32705e-12 MeV*mm
//
static const double h_Planck      = 6.62606896e-34 * joule*s;
static const double hbar_Planck   = h_Planck/twopi;
static const double hbarc         = hbar_Planck * c;
static const double hbarc2 = hbarc * hbarc;

//
//
//
static const double electron_charge = - eplus; // see SystemOfUnits.h
static const double e_squared = eplus * eplus;

//
// amu_c2 - atomic equivalent mass unit
//        - AKA, unified atomic mass unit (u)
// amu    - atomic mass unit
//
static const double electron_mass_c2 = 0.510998910 * MeV;
static const double   proton_mass_c2 = 938.272013 * MeV;
static const double  neutron_mass_c2 = 939.56536 * MeV;
static const double           amu_c2 = 931.494028 * MeV;
static const double              amu = amu_c2/c2;