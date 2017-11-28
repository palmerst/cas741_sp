import re

from Species import *
from ChemEq import *
from Calculations import *
from Plot import *

class ChemSys:
  
  def __init__(self):
    self.sN = 1
    self.species = {}
    self.rxns = []
    self.totals = {}
    
    self.species["(H)+"] = Species("(H)+", "aq", 1, { "H" : 1 }, -1)
    self.species["(OH)-"] = Species("(OH)-", "aq", -1, { "H" : 1, "O" : 1 }, -2)
    self.species["(H2O)l"] = Species("(H2O)l", "l", 0, { "H" : 2, "O" : 1 }, -3)

    
  def specGen(self, plotName):
    (phRange, outData, labels) = calcSpec(self)
    genDiagram(plotName, phRange, "pH", outData, "conc", labels)  

   
  def registerTotal(self, elt, tot):
    if tot <= 0:
      raise SystemExit("Element total concentration must be greater than 0.")    
    self.totals[elt] = tot
    
    
  def registerRxn(self, rString):
    rxn , logK = rString.split(",")
    left , right = rxn.split(" = ")
    reactants = left.split(" + ")
    products = right.split(" + ")
    
    ## strip any whitespace
    logK = logK.strip()
    reactants = [s.strip() for s in reactants]
    products = [s.strip() for s in products]
    
    r = {}
    p = {}
    
    ## store reactants + numbers in r dict
    for reactant in reactants:
      sp, num = self.__parseRxnSpecies(reactant)      
      r[sp] = num
    
    ## store products + numbers in p dict    
    for product in products:
      sp, num = self.__parseRxnSpecies(product)      
      p[sp] = num
      
    logK = float(logK)
    
    self.rxns.append(ChemEq(rxn, logK, r, p))


    
  ## break down reaction species
  def __parseRxnSpecies(self, rsp):
    m = re.match('^([1-9]?|[1-9][0-9]*)(.*)$', rsp)    
    
    if m is None:
      raise SystemExit("Reaction species %s is ill-formed" % rsp)
    
    num = m.group(1)
    if num == "":
      num = 1
    
    sp = m.group(2)
    
    ## parse and store the species
    self.__addSpecies(sp)

    return (sp, int(num))
    
  
  def __addSpecies(self, spString):
    ## ignore H+, OH-, H2O
    if spString in ["(H)+", "(OH)-", "(H2O)l"]:
      return
  
    ## if we've stored this species already, just retrieve it
    if spString in self.species:
      return
    
    ## otherwise, parse and store    
    formula = spString   
    (state, charge, composition) = self.__parseSpecies(spString)
   
    self.species[spString] = Species(formula, state, charge, composition, self.sN)
    self.sN += 1
    
    return
 
 
  def __parseSpecies(self, str):
    m = re.match('^\((.*)\)(?:(([1-9]?|[1-9][0-9]*)([\+-]))|([lgs]|aq))$', str)
    
    if m is None:
      raise SystemExit("Species %s is ill-formed (state)" % str)
    
    sp = m.group(1)
    
    ## if there is a charge, get charge
    if m.group(2) is not None: 
      (mag, sign) = (m.group(3), m.group(4))
    
      if mag == "":
        charge = 1
      else:
        charge = int(mag)
    
      if sign == "-":
        charge *= -1
      
      ## charged species must be in solution        
      state = "aq"
  
    ## if no charge, then get state
    else:
      charge = 0
      state = m.group(5)
      
    comp = self.__parseComposition(sp, 1)
        
    return (state, charge, comp)
    
 
  ## parse the formula to get element totals
  def __parseComposition(self, component, n):
    if component == "":
      return {}

    ## flag for parenthesized component
    pFlag = False

    ## if next component is in brackets
    ##   we must manually find what is between brackets
    ##   (regex can't do this as it has no state)
    if component[0] == '(':
      parens = 1
      l = 0
      
      for c in component[1:]:
        if c == '(':
          parens += 1
        elif c == ')':
          parens -= 1
        
        if parens == 0:
          break
        
        l += 1
      
      ## report error if unbalanced parens, or nothing between parens
      if parens != 0 or l == 0:
        raise SystemExit("Species %s is ill-formed (parentheses)" % str)
      
      ## take what is between parens          
      c = component[1:1 + l] 
      
      ## regex to find number of occurrences and rest of formula
      m = re.match('^([1-9][0-9]*)?(.*)$', component[2 + l:]) 
      (num, rest) = (m.group(1), m.group(2))
      
      pFlag = True
    
    ## otherwise just use full regex
    else:                
      m = re.match('^([A-Z][a-z]?)([1-9][0-9]*)?(.*)$', component)
    
      if m is None:
        raise SystemExit("Species %s is ill-formed (formula)" % str)
        
      (c, num, rest) = (m.group(1), m.group(2), m.group(3))
    
    
    dictRest = {}
    
    ## if there is more to the formula, 
    ##   then parse it and store totals in dictRest
    if rest is not None:
      dictRest = self.__parseComposition(rest, n)
    
    
    if num is not None:
      dictComp = self.__parseComposition(c, n * int(num))
    ## if there is no number and component was parenthesized
    ##   then there is still parsing to do
    elif pFlag: 
      dictComp = self.__parseComposition(c, n)
    else:
      dictComp = { c: n }
    
    ## add the total of the current component 
    ##   into the dictionary  
    for key, value in dictComp.items():
      dictRest[key] = dictRest.get(key, 0) + value
    
    return dictRest