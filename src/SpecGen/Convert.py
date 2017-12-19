from ast import *
from types import *
from functools import reduce
from math import log10

from SpecGen.ChemSys import *
from SpecGen.Species import *
from SpecGen.ChemEq import *

def makeRootFunc(syst):  
  def toInput():
    specs = syst.species.values()
    vars = []
    vlist = []
    for spec in specs:
      if spec.state != "aq":
        continue
      if spec.number < 0:
        continue
      vars.append(Name("x" + str(spec.number), Store()))
      vlist.append(spec.formula)
    return (vars , vlist)
  
  
  (vars , vlist) = toInput()
  
  if not vlist:
    raise RuntimeError("Chemical system contains no free aqueous species.")  

  equil = [equilibrium(r, syst.species) for r in syst.rxns]
  
  mBal = [massBalance(elt, tot, syst.species) for (elt, tot) in syst.totals.items()]
   
  results = equil + mBal 
  
  varList = [Tuple(vars, Store())]  
  resultList = Tuple(results, Load())   
    
  eqn_ast = FunctionDef(
    name='rootFunc',
    args=arguments(args=[arg('p', None)], vararg=arg('arg', None), kwonlyargs=[], kw_defaults=[], kwarg=None, defaults=[]),
    body=[Assign([Tuple([Name("x" + str(syst.species["(H)+"].number), Store()),Name("x" + str(syst.species["(OH)-"].number), Store())],Store())], Name('arg', Load())), Assign(varList, Name('p', Load())), Return(value=resultList)],
    decorator_list=[],
    lineno=1,
    col_offset=0
  )

  eqn_mod_ast = Module(body=[ImportFrom(module='math', names=[alias(name='log10', asname=None)]), eqn_ast])

  eqn_mod_ast = fix_missing_locations(eqn_mod_ast)
  
  ## compile the AST into bytecode so we can use the function
  eqn_code = compile(eqn_mod_ast, "<not_a_file>", "exec")
  ## get the function bytecode
  rootFunc_code = [c for c in eqn_code.co_consts if isinstance(c, CodeType)][0]
  
  ## create function from bytecode
  return (FunctionType(rootFunc_code, {'log10': log10}) , vlist)
    
def equilibrium(r, species):
  def getTerm(sStr, n):
    spec = species[sStr]
    if spec.state != "aq":
      return Num(0)  # or None and filter?
    term = Name("x" + str(spec.number), Load())
    term = BinOp(Num(n), Mult(), term)
    return term
  
  def combine(rxnPairs):
   # l = [(k, v) for (k, v) in rxnPairs]
    l = list(map(lambda x: getTerm(*x), rxnPairs))
    return reduce(lambda x, y: BinOp(x, Add(), y), l, Num(0))
  
  pdt = combine(r.p.items())
  rct = combine(r.r.items())
   
  logK = BinOp(pdt, Sub(), rct)
  
  return BinOp(logK, Sub(), Num(r.logK))
  
def massBalance(elt, tot, species):
  def getTerm(spec):
    if spec.state != "aq":
      return Num(0)  # or None and filter?
    n = spec.composition.get(elt)
    if n is None:
      return Num(0)
    term = BinOp(Num(10), Pow(), Name("x" + str(spec.number), Load()))
    term = BinOp(Num(n), Mult(), term)
    return term
    
  specs = species.values()
  massTerms = list(map(getTerm, specs))
  massSum = reduce(lambda x, y: BinOp(x, Add(), y), massTerms, Num(0)) 
  
  logMassSum = Call(func=Name(id='log10', ctx=Load()), args=[massSum], keywords=[])
  logTot = Call(func=Name(id='log10', ctx=Load()), args=[Num(tot)], keywords=[])
  
  return BinOp(logMassSum, Sub(), logTot)
