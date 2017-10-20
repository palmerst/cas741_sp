from Species import *
from ChemEq import *

r1 = "(H2O)l = (H)+ + (OH)- , -14"
r2 = "(Fe)3+ + (H2O)l = (Fe(OH))2+ + (H)+ , -2.19"
r3 = "(Fe)3+ + 2(H2O)l = (Fe(OH)2)+ + 2(H)+ , -5.67"
r4 = "(Fe)3+ + 4(H2O)l = (Fe(OH)4)- + 4(H)+ , -21.6"

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