/* const/test.c
 * 
 * Copyright (C) 2003, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const.h>
#include <gsl/gsl_test.h>

#include <gsl/gsl_ieee_utils.h>

int
main (void)
{
  gsl_ieee_env_setup ();

  /* Basic check to make sure the header files are functioning */

  {
    double c = GSL_CONST_MKS_SPEED_OF_LIGHT;
    double eps = GSL_CONST_MKS_VACUUM_PERMITTIVITY;
    double mu = GSL_CONST_MKS_VACUUM_PERMEABILITY;

    gsl_test_rel (c, 1.0/sqrt(eps*mu), 1e-6, "speed of light (mks)");
  }

  {
    double ly = GSL_CONST_CGS_LIGHT_YEAR;
    double c = GSL_CONST_CGS_SPEED_OF_LIGHT;
    double y = 365.2425 * GSL_CONST_CGS_DAY;
    
    gsl_test_rel (ly, c * y, 1e-6, "light year (cgs)");
  }

  {
    double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
    double eps = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
    double mu = GSL_CONST_MKSA_VACUUM_PERMEABILITY;

    gsl_test_rel (c, 1.0/sqrt(eps*mu), 1e-6, "speed of light (mksa)");
  }

  {
    double ly = GSL_CONST_CGSM_LIGHT_YEAR;
    double c = GSL_CONST_CGSM_SPEED_OF_LIGHT;
    double y = 365.2425 * GSL_CONST_CGSM_DAY;
    
    gsl_test_rel (ly, c * y, 1e-6, "light year (cgsm)");
  }

  {
    double micro = GSL_CONST_NUM_MICRO;
    double mega = GSL_CONST_NUM_MEGA;
    double kilo = GSL_CONST_NUM_KILO;

    gsl_test_rel (mega/kilo, 1/(micro*kilo), 1e-10, "kilo (mega/kilo, 1/(micro*kilo))");
  }

  {
    double d = GSL_CONST_MKSA_DEBYE;
    double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
    double desu = d * c * 1000.0;
    
    gsl_test_rel (desu, 1e-18, 1e-10, "debye (esu)");
  }

  {
    double k = GSL_CONST_MKSA_BOLTZMANN;
    double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
    double h = GSL_CONST_MKSA_PLANCKS_CONSTANT_H;
    double s = 2 * pow(M_PI, 5.0) * pow(k, 4.0) / (15 * pow(c, 2.0) * pow(h, 3.0));
    double sigma = GSL_CONST_MKSA_STEFAN_BOLTZMANN_CONSTANT;
    
    gsl_test_rel(s, sigma, 1e-10, "stefan boltzmann constant");
  }

  /* CODATA 2022 consistency checks (see const/gsl_const_mksa.h) */

  {
    double h = GSL_CONST_MKSA_PLANCKS_CONSTANT_H;
    double hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;

    gsl_test_rel (hbar, h / (2.0 * M_PI), 1e-14, "hbar = h / 2 pi");
  }

  {
    double e = GSL_CONST_MKSA_ELECTRON_CHARGE;
    double F = GSL_CONST_MKSA_FARADAY;

    gsl_test_rel (F, GSL_CONST_NUM_AVOGADRO * e, 1e-9, "faraday = N_A e");
  }

  {
    double R = GSL_CONST_MKSA_MOLAR_GAS;
    double V = GSL_CONST_MKSA_STANDARD_GAS_VOLUME;

    gsl_test_rel (V, R * 273.15 / 100000.0, 1e-9, "standard gas volume (273.15 K, 100 kPa)");
  }

  {
    double u = GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS;

    gsl_test_rel (u * GSL_CONST_NUM_AVOGADRO, 1e-3, 1e-8, "unified atomic mass * N_A = 1 g/mol");
  }

  {
    gsl_test_rel (GSL_CONST_MKSA_ELECTRON_VOLT, GSL_CONST_MKSA_ELECTRON_CHARGE, 1e-15, "electron volt = e J");
  }

  {
    /* masses relative to the unified atomic mass unit */
    gsl_test_rel (GSL_CONST_MKSA_MASS_ELECTRON / GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS, 5.485799090441e-4, 1e-9, "electron mass in u");
    gsl_test_rel (GSL_CONST_MKSA_MASS_PROTON / GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS, 1.0072764665789, 1e-9, "proton mass in u");
    gsl_test_rel (GSL_CONST_MKSA_MASS_NEUTRON / GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS, 1.00866491605, 1e-9, "neutron mass in u");
    gsl_test_rel (GSL_CONST_MKSA_MASS_MUON / GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS, 0.1134289257, 1e-8, "muon mass in u");
  }

  {
    double a = GSL_CONST_MKSA_RADIATION_DENSITY_CONSTANT;
    double sigma = GSL_CONST_MKSA_STEFAN_BOLTZMANN_CONSTANT;
    double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;

    gsl_test_rel (a, 4 * sigma / c, 1e-12, "radiation density constant (4 sigma / c)");
  }

  {
    /* IUGG mean Earth radii for the WGS-84 ellipsoid */
    double a = GSL_CONST_MKSA_EARTH_EQUATORIAL_RADIUS_A;
    double b = GSL_CONST_MKSA_EARTH_POLAR_RADIUS_B;
    double R1 = GSL_CONST_MKSA_EARTH_MEAN_RADIUS_R1;

    gsl_test_rel (R1, (2.0 * a + b) / 3.0, 1e-15, "earth mean radius R1 = (2a + b)/3");
    gsl_test_rel (R1, 6371008.7714, 1e-9, "earth mean radius R1 (6371.0087714 km)");
    gsl_test_rel (GSL_CONST_MKSA_EARTH_AUTHALIC_RADIUS_R2, 6371007.1809, 1e-9, "earth authalic radius R2");
    gsl_test_rel (GSL_CONST_MKSA_EARTH_VOLUMETRIC_RADIUS_R3, 6371000.79, 1e-9, "earth volumetric radius R3");
  }

  /* the cgs and cgsm constants must agree with their mks counterparts */

  {
    gsl_test_rel (GSL_CONST_CGS_SPEED_OF_LIGHT, 1e2 * GSL_CONST_MKS_SPEED_OF_LIGHT, 1e-15, "speed of light (cgs vs mks)");
    gsl_test_rel (GSL_CONST_CGS_GRAVITATIONAL_CONSTANT, 1e3 * GSL_CONST_MKS_GRAVITATIONAL_CONSTANT, 1e-15, "gravitational constant (cgs vs mks)");
    gsl_test_rel (GSL_CONST_CGS_BOHR_RADIUS, 1e2 * GSL_CONST_MKS_BOHR_RADIUS, 1e-15, "bohr radius (cgs vs mks)");
    gsl_test_rel (GSL_CONST_CGS_STEFAN_BOLTZMANN_CONSTANT, 1e3 * GSL_CONST_MKS_STEFAN_BOLTZMANN_CONSTANT, 1e-15, "stefan boltzmann constant (cgs vs mks)");
    gsl_test_rel (GSL_CONST_CGS_EARTH_MEAN_RADIUS_R1, 1e2 * GSL_CONST_MKS_EARTH_MEAN_RADIUS_R1, 1e-15, "earth mean radius R1 (cgs vs mks)");
    gsl_test_rel (GSL_CONST_CGSM_BOHR_MAGNETON, 1e3 * GSL_CONST_MKS_BOHR_MAGNETON, 1e-15, "bohr magneton (cgsm vs mks)");
    gsl_test_rel (GSL_CONST_CGSM_ELECTRON_CHARGE, 1e-1 * GSL_CONST_MKS_ELECTRON_CHARGE, 1e-15, "electron charge (cgsm vs mks)");
    gsl_test_rel (GSL_CONST_CGSM_RADIATION_DENSITY_CONSTANT, 1e1 * GSL_CONST_MKS_RADIATION_DENSITY_CONSTANT, 1e-15, "radiation density constant (cgsm vs mks)");
  }

  /* pin representative values to their CODATA 2022 published values, so
     that an accidental change or a partial update is caught */

  {
    gsl_test_rel (GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT, 6.67430e-11, 1e-9, "CODATA 2022 gravitational constant");
    gsl_test_rel (GSL_CONST_NUM_FINE_STRUCTURE, 7.2973525643e-3, 1e-9, "CODATA 2022 fine structure constant");
    gsl_test_rel (GSL_CONST_MKSA_MASS_ELECTRON, 9.1093837139e-31, 1e-9, "CODATA 2022 electron mass");
    gsl_test_rel (GSL_CONST_MKSA_RYDBERG, 2.1798723611030e-18, 1e-9, "CODATA 2022 Rydberg constant (R_inf h c)");
    gsl_test_rel (GSL_CONST_MKSA_VACUUM_PERMITTIVITY, 8.8541878188e-12, 1e-9, "CODATA 2022 vacuum permittivity");
    gsl_test_rel (GSL_CONST_MKSA_VACUUM_PERMEABILITY, 1.25663706127e-6, 1e-9, "CODATA 2022 vacuum permeability");
    gsl_test_rel (GSL_CONST_MKSA_BOHR_RADIUS, 5.29177210544e-11, 1e-9, "CODATA 2022 Bohr radius");
    gsl_test_rel (GSL_CONST_MKSA_THOMSON_CROSS_SECTION, 6.6524587051e-29, 1e-9, "CODATA 2022 Thomson cross section");
  }

  exit (gsl_test_summary ());
}

