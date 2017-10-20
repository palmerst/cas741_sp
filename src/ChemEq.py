import re
from Species import *

## reactions list
rxns = []

class ChemEq:
  def __init__(self, rString):
    rxn , logK = rString.split(",")
    left , right = rxn.split(" = ")
    reactants = left.split(" + ")
    products = right.split(" + ")
    
    ## strip any whitespace
    logK = logK.strip()
    reactants = [s.strip() for s in reactants]
    products = [s.strip() for s in products]
    
    self.r = {}
    self.p = {}
    
    ## break down reaction species
    def parseRxnSpecies(rsp):
      m = re.match('^([1-9]?|[1-9][0-9]*)(.*)$', rsp)    
      
      if m is None:
        raise SystemExit("Reaction species %s is ill-formed" % rsp)
      
      num = m.group(1)
      if num == "":
        num = 1
      
      sp = m.group(2)
      
      ## parse and store the species
      Species(sp)
    
      return (sp, int(num))
    
    ## store reactants + numbers in r dict
    for reactant in reactants:
      sp, num = parseRxnSpecies(reactant)      
      self.r[sp] = num
    
    ## store products + numbers in p dict    
    for product in products:
      sp, num = parseRxnSpecies(product)      
      self.p[sp] = num
      
    self.K = 10**float(logK)
    self.reaction = rxn  # store reaction string for printing
    
    rxns.append(self)