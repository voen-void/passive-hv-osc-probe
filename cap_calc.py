# this is a very, *very* rough estimation of the internal capacitance between the voltage divider, and the shielding tube.
# link to github repo of the project that this script has been written for: https://github.com/voen-void/passive-hv-osc-probe/tree/main
# formula from: https://testbook.com/physics-formulas/cylindrical-capacitor-formula
# this calculation divides the voltage divider into n (seg_num) segments which are assumed to be of equal length.
# each such segment is assumed to be a pair of coaxial cylinders, representing a single resistor - capacitor pair of the divider
# script written by Jan Klimas, 2025

import math

e0 = 8.854187817E-12 # [ F / m ]   dielectric constant of vacuum
er = 6               # [ no unit ] relative dielectric constant of the resin used for casting
l  = 0.026           # [ m ]       length of a single segment must be meters!
b  = 37.6            # [ mm ]      inner diameter of the shielding tube, a and b must share the same unit.
a  = 10              # [ mm ]      outer diameter of the inner cylinder a and b must share the same unit.
seg_num = 5          # [ no unit ] number of segments - 5 in the case of this project

# note that the inner cylinder (a) is not physical but a way to approximate the total effective area of the inner conductor, this value has been guesstimated.


# calculating the capacitance of a single segment
C_seg = ( 2 * math.pi * e0 * er * l ) / math.log( b / a )
# the capacitance is multiplied by 1E+12 before printing to print out the value in pF which is more readable
print( "Capacitance of a single segment of the divider:" , C_seg * 1E+12 , "pF" )

# total capacitance - these capacitances add up, as they are connected in parallel
C_tot = C_seg * seg_num
# the capacitance is multiplied by 1E+12 before printing to print out the value in pF which is more readable
print( "Capacitance of" , seg_num , "consecutive segments:" , C_tot * 1E+12 , "pF" )