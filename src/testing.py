from ast import *
from types import *

from scipy.optimize import fsolve
from numpy import float64

## test:  solving
##  [h+][oh-] - 10^-14 = 0
##  [h+] - 10 ** -ph = 0

class Species:
  def __init__(self):
    self.symb = ""
    self.charge = 0

class ChemEq:
  def __init__(self):
    self.K = 0
  
H = Species()
H.symb = "H"
H.charge = 1

OH = Species()
OH.symb = "OH"
OH.charge = -1

waterDissoc = ChemEq()
waterDissoc.K = 10**-14


##############################################
## START GENERATING FUNCTION TO BE SOLVED
##
## Using python AST to build equations -- this can be used to generate functions on the fly
##   generated functions can be passed to non-linear solver

## ph equation
phEqn = BinOp(Name(H.symb, Load()), Sub(), BinOp(Num(10), Pow(), UnaryOp(USub(), Name('ph', Load()))))

## will generate these from ChemEq instances
waterDissocEqn = BinOp(BinOp(Name(H.symb, Load()), Mult(), Name(OH.symb, Load())), Sub(), Num(waterDissoc.K))
  
results = Tuple([phEqn, waterDissocEqn], Load())
  
varList = [Tuple([Name(H.symb, Store()), Name(OH.symb, Store())], Store())] 
  
eqn_ast = FunctionDef(
    name='equations',
    args=arguments(args=[arg('p', None), arg('ph', None)], vararg=None, kwonlyargs=[], kw_defaults=[], kwarg=None, defaults=[]),
    body=[Assign(varList, Name('p', Load())), Return(value=results)],
    decorator_list=[],
    lineno=1,
    col_offset=0
)

eqn_ast = fix_missing_locations(eqn_ast)

eqn_mod_ast = Module(body=[eqn_ast])

## compile the AST into bytecode so we can use the function
eqn_code = compile(eqn_mod_ast, "<not_a_file>", "exec")
## get the function bytecode
feq_code = [c for c in eqn_code.co_consts if isinstance(c, CodeType)][0]
## create function from bytecode
f = FunctionType(feq_code, {})

## we have f now;  f can be passed to fsolve

## END GENERATION 
##############################################


res = [fsolve(f, (10**-7,10**-7), ph) for ph in range(0, 15)]

for (h, oh) in res:
  print(h, oh)
  
print("\n\n")

a = parse("\
from math import log10\n\n\
def f():\n\
  return log10(3)")
  
a = "\
from math import log10\n\n\
def f():\n\
  return log10(3)"

exec(a, globals())  
  
#a_code = compile(a, "<not_a_file>", "exec")
## get the function bytecode
#f_code = [c for c in a_code.co_consts if isinstance(c, CodeType)][0]

#print(a_code.co_consts)

## create function from bytecode
#f = FunctionType(f_code, {})

print(f())
