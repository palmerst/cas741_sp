from ast import *
from types import *
from functools import reduce

from Species import *
from ChemEq import *

def makeRootFunc(elt, tot):  # element totals need a home
  def toAbs(spec):
    if spec.state != "aq":
      return None
    return Name("x" + str(spec.number), Store())
    
  specs = species.values()
  vars = list(map(toAbs, specs))
  vars = [x for x in vars if x is not None]

  equil = [equilibrium(r) for r in rxns]
  
  phEqn = BinOp(Name("x" + str(species["(H)+"].number), Load()), Sub(), BinOp(Num(10), Pow(), UnaryOp(USub(), Name('ph', Load()))))
   
  results = equil + [massBalance(elt, tot)] + [chargeBalance()] + [phEqn]
  
  extra = len(results) - len(vars)
  vars = vars + [Name('_', Store())] * extra  
  
  varList = [Tuple(vars, Store())]  
  resultList = Tuple(results, Load())   
    
  eqn_ast = FunctionDef(
    name='rootFunc',
    args=arguments(args=[arg('p', None)], vararg=arg('arg', None), kwonlyargs=[], kw_defaults=[], kwarg=None, defaults=[]),
    body=[Assign([Tuple([Name('ph', Store())],Store())], Name('arg', Load())), Assign(varList, Name('p', Load())), Return(value=resultList)],
    decorator_list=[],
    lineno=1,
    col_offset=0
  )

  eqn_ast = fix_missing_locations(eqn_ast)

  eqn_mod_ast = Module(body=[eqn_ast])

  ## compile the AST into bytecode so we can use the function
  eqn_code = compile(eqn_mod_ast, "<not_a_file>", "exec")
  ## get the function bytecode
  rootFunc_code = [c for c in eqn_code.co_consts if isinstance(c, CodeType)][0]
  
  ## create function from bytecode
  return FunctionType(rootFunc_code, {})
    
def equilibrium(r):
  def getTerm(sStr, n):
    spec = species[sStr]
    if spec.state != "aq":
      return Num(1)  # or None and filter?
    term = Name("x" + str(spec.number), Load())
    if n > 1:
      term = BinOp(term, Pow(), Num(n))
    return term
  
  def combine(rxnPairs):
   # l = [(k, v) for (k, v) in rxnPairs]
    l = list(map(lambda x: getTerm(*x), rxnPairs))
    return reduce(lambda x, y: BinOp(x, Mult(), y), l, Num(1))
  
  numer = combine(r.p.items())
  denom = combine(r.r.items())
   
  K = BinOp(numer, Div(), denom)
  
  return BinOp(K, Sub(), Num(r.K))
  
def massBalance(elt, tot):
  def getTerm(spec):
    if spec.state != "aq":
      return Num(0)  # or None and filter?
    n = spec.composition.get(elt)
    if n is None:
      return Num(0)
    term = Name("x" + str(spec.number), Load())
    if n > 1:
      term = BinOp(Num(n), Mult(), term)
    return term
    
  specs = species.values()
  massTerms = list(map(getTerm, specs))
  massSum = reduce(lambda x, y: BinOp(x, Add(), y), massTerms, Num(0)) 
  
  return BinOp(massSum, Sub(), Num(tot))

def chargeBalance():
  def getTerm(spec):
    if spec.state != "aq":
      return Num(0)  # or None and filter?
    n = spec.charge
    if n == 0:
      return Num(0)
    term = Name("x" + str(spec.number), Load())
    return BinOp(Num(n), Mult(), term)
  
  specs = species.values()
  chargeTerms = list(map(getTerm, specs))
  return reduce(lambda x, y: BinOp(x, Add(), y), chargeTerms, Num(0)) 
 

def positiveConc():
  def getTerm(spec):
    if spec.state != "aq":
      return None
    n = spec.charge
    if n == 0:
      return Num(0)
    term = Name("x" + str(spec.number), Load())
    return BinOp(Num(n), Mult(), term)
  
  specs = species.values()
  chargeTerms = list(map(getTerm, specs))
