from SpecGen.Convert import *
from SpecGen.ChemSys import *

from scipy.optimize import *

def calcSpec(syst):
  (f , l) = makeRootFunc(syst)

  phRange = [ph/100 for ph in range(0, 1401)]

  results = [fsolve(f, tuple([-5.0]*len(l)), (-ph, -(14-ph)), xtol=1e-03) for ph in phRange]

  outData = [[] for _ in range(len(l))]

  for res in results:
    for i in range(0,len(res)):
      outData[i].append(10**res[i])
   
  return (phRange, outData, l)

  
  