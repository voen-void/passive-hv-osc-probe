import math

# this is a very, *very* rough estimation of the internal capacitance between the voltage divider, and the shielding tube.
# formula from: https://testbook.com/physics-formulas/cylindrical-capacitor-formula

e0 = 8.854187817E-12 # [ F / m ]
er = 6
l  = 0.026           # [ m ]  must be meters!
b  = 37.6            # [ mm ] a, b must share the same unit.
a  = 7               # [ mm ] a, b must share the same unit.

C_seg = ( 2 * math.pi * e0 * er * l ) / math.log( b / a )
print( C_seg )

seg_num = 5

C_tot = C_seg * seg_num
print( C_tot )