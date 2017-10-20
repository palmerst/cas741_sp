from Species import *
from ChemEq import *
from Convert import *

from scipy.optimize import *

r1 = "(H2O)l = (H)+ + (OH)- , -14"
r2 = "(Fe)3+ + (H2O)l = (Fe(OH))2+ + (H)+ , -2.19"
r3 = "(Fe)3+ + 2(H2O)l = (Fe(OH)2)+ + 2(H)+ , -5.67"
r4 = "(Fe)3+ + 4(H2O)l = (Fe(OH)4)- + 4(H)+ , -21.6"

ChemEq(r1)
ChemEq(r2)
ChemEq(r3)
ChemEq(r4)

f = makeRootFunc("Fe", 0.000010)

res = [fsolve(f, (1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,0), ph) for ph in range(0, 15)]
#res = [differential_evolution(f, ((1e-50,1),(1e-50,1),(1e-50,1),(1e-50,1),(1e-50,1),(1e-50,1),(1e-50,1)), (ph,)) for ph in range(1, 15)]

#res = [leastsq(f, (1e-4,1e-4,1e-4,1e-4,1e-4,1e-4), ph) for ph in range(0, 15)]

for (h, oh, fe3, feoh, feoh2, feoh3, _) in res:
  print(h, oh, fe3, feoh, feoh2, feoh3)

'''
ChemEq(r1)
ChemEq(r2)
ChemEq(r3)
ChemEq(r4)

for r in rxns:
  print(r.reaction)
  print(r.K)
  print(r.r)
  print(r.p)
  print("\n\n")
  
for (_ , s) in species.items():
  print(s.formula)
  print(s.state)
  print(s.charge)
  print(s.composition)
  print("\n\n")
  
eqns = [equilibrium(r) for r in rxns]

for e in eqns:
  print(dump(e))
  
print("\n")
print(dump(massBalance("Fe", 0.000010)))

print("\n")
print(dump(chargeBalance()))
'''