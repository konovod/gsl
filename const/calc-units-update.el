;; calc-units-update.el -- update for calc-units.el, a calculator for Emacs
;;
;; Copyright (C) 2002,2008 Jochen K�pper
;;
;; This file is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published
;; by the Free Software Foundation; either version 3, or (at your
;; option) any later version.
;;
;; It is distributed in the hope that it will be useful, but WITHOUT
;; ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
;; or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
;; License for more details.
;;
;; These values are taken from the following sources:
;; * CODATA 2022
;;   - https://physics.nist.gov/cuu/Constants/
;;   - https://physics.nist.gov/cuu/Constants/Table/allascii.txt
;;   - CODATA 2006 (superseded values kept above for reference)
;; * NIST
;;   - http://physics.nist.gov/Pubs/SP811/appenB9.html
;; * NASA JPL
;;   - http://neo.jpl.nasa.gov/glossary/au.html

(require 'calc)

(setq math-additional-units
  '(;; length
    ( au      "149597870691. m"        "Astronomical Unit" ) ;; NASA JPL 
    ;; mass
    ( amu     "1.66053906892e-27 kg"   "Unified atomic mass" ) ;; CODATA 2022
    ;; pressure
    ( inH2O   "2.490889e2 Pa"          "Inch of water" ) ;; NIST 
    ( ech     "1.602176634e-19 C"      "Elementary charge" ) ;; exact (SI 2019)
    ( e       "ech"                    "Elementary charge" )
    ;; other physical quantities
    ( h       "6.62607015e-34 J s"     "*Planck's constant" ) ;; exact (SI 2019)
    ( hbar    "h / 2 pi"               "Planck's constant" )
    ( mu0     "1.25663706127e-6 H/m"  "Permeability of vacuum" ) ;; CODATA 2022
    ( G       "6.67430e-11 m^3/kg^1/s^2" "Gravitational constant" ) ;; CODATA 2022
    ( Nav     "6.02214076e23 / mol"    "Avagadro's constant" ) ;; exact (SI 2019)
    ( me      "9.1093837139e-31 kg"    "Electron rest mass" ) ;; CODATA 2022
    ( mp      "1.67262192595e-27 kg"   "Proton rest mass" ) ;; CODATA 2022
    ( mn      "1.67492750056e-27 kg"   "Neutron rest mass" ) ;; CODATA 2022
    ( mmu     "1.883531627e-28 kg"     "Muon rest mass" ) ;; CODATA 2022
    ( Ryd     "10973731.568157 h c/m"  "Rydberg's constant (energy)" ) ;; CODATA 2022 R_inf
    ( k       "1.380649e-23 J/K"      "Boltzmann's constant" ) ;; exact (SI 2019)
    ( fsc     "7.2973525643e-3"        "Fine structure constant" ) ;; CODATA 2022
    ( muB     "9.2740100657e-24 J/T"   "Bohr magneton" ) ;; CODATA 2022
    ( muN     "5.0507837393e-27 J/T"   "Nuclear magneton" ) ;; CODATA 2022
    ( mue     "9.2847646917e-24 J/T"   "Electron magnetic moment (absolute value)" ) ;; CODATA 2022
    ( mup     "1.41060679545e-26 J/T"  "Proton magnetic moment" ) ;; CODATA 2022
    ( R0      "8.314462618 J/mol/K"    "Molar gas constant" ) ;; exact (SI 2019)
    ( V0      "22.71095464e-3 m^3/mol" "Standard volume of ideal gas" ) ;; exact (SI 2019)
    ( flam    "1.07639104e-3 lam"      "Footlambert" )
    ( Torr    "atm/760"                "Torr" )
    ( fc      "10.76 lx"               "Footcandle" ) ;; m^2/(ft^2), not full accurate
 ))


(setq math-standard-units
      (append math-standard-units
              '(( abamp       "10 A"                   "*Abampere" ))))

(setq math-standard-units-systems 
      (append math-standard-units-systems
              '((cgsm ((m '(* (var cm var-cm) 100))
                       (A '(/ (var abamp var-abamp) 10))))
                (mksa ((g '(* (var kg var-kg) (float 1 -3))))))))

(setq math-units-table nil)

(provide 'calc-units-update)
