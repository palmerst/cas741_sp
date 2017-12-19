import pytest
import filecmp
import os

from SpecGen.ChemSys import *
from SpecGen.Convert import *
from SpecGen.Calculations import *
from SpecGen.Plot import *


# systems for testing
emptySys = ChemSys()

feSys = ChemSys()
feSys.registerRxn("(Fe)3+ + (H2O)l = (Fe(OH))2+ + (H)+ , -2.19")
feSys.registerRxn("(Fe)3+ + 2(H2O)l = (Fe(OH)2)+ + 2(H)+ , -5.67")
feSys.registerRxn("(Fe)3+ + 4(H2O)l = (Fe(OH)4)- + 4(H)+ , -21.6")
feSys.registerTotal("Fe", 0.000010)

co3Sys = ChemSys()
co3Sys.registerRxn("(CO2)aq + (H2O)l = (H2CO3)aq , -2.77")
co3Sys.registerRxn("(H2CO3)aq = (HCO3)- + (H)+ , -3.6")
co3Sys.registerRxn("(HCO3)- = (CO3)2- + (H)+ , -10.33")
co3Sys.registerTotal("C", 0.000010)



###########################
# Plotting module testing #
###########################
class TestPlot(object):
  
  # test plotting line vs known diagram
  def test_plot(self):
    genDiagram("test_plot", range(0, 11), "x", [range(0, 11)], "y", ["y1"])
    assert filecmp.cmp("test_plot.png", "./testref/test_plot_ref.png")
    os.remove("test_plot.png")


###################################
# Input Conversion module testing #
###################################
class TestConvert(object):
  
  # for use in testing function equivalence
  # this is not an actual test.
  def expectedFe(self, p, *arg):
    (h, oh) = arg
    (x1, x2, x3, x4) = p
    return ( x2 + h - x1 + 2.19, 
           x3 + 2*h - x1 + 5.67,
           x4 + 4*h - x1 + 21.6,
           log10(10**x1 + 10**x2 + 10**x3 + 10**x4) - log10(0.000010) )
           
  # for use in testing function equivalence
  # this is not an actual test.
  def expectedCo3(self, p, *arg):
    (h, oh) = arg
    (x1, x2, x3, x4) = p
    return ( x2 - x1 + 2.77, 
           x3 + h - x2 + 3.6,
           x4 + h - x3 + 10.33,
           log10(10**x1 + 10**x2 + 10**x3 + 10**x4) - log10(0.000010) )

  # the best we can do here is test for extensional equality for some inputs
  # testing over the range -14 to 0 with step 2 for all inputs
  # this corresponds to concentrations between 10^-14 M and 1 M  
  def test_convert_fe(self):
    (fn, _) = makeRootFunc(feSys)
    
    for a in range(-14, 1, 2):
      for b in range(-14, 1, 2):
        for c in range(-14, 1, 2):
          for d in range(-14, 1, 2):
            for e in range(-14, 1, 2):
              assert fn((a, b, c, d), e, -(14+e)) == self.expectedFe((a, b, c, d), e, -(14+e)) 
    
  # the best we can do here is test for extensional equality for some inputs
  # testing over the range -14 to 0 with step 2 for all inputs
  # this corresponds to concentrations between 10^-14 M and 1 M  
  def test_convert_co3(self):
    (fn, _) = makeRootFunc(co3Sys)
    
    for a in range(-14, 1, 2):
      for b in range(-14, 1, 2):
        for c in range(-14, 1, 2):
          for d in range(-14, 1, 2):
            for e in range(-14, 1, 2):
              assert fn((a, b, c, d), e, -(14+e)) == self.expectedCo3((a, b, c, d), e, -(14+e))  
  
  # conversion for empty system should throw  
  def test_convert_exc(self):
    with pytest.raises(RuntimeError):
      makeRootFunc(emptySys) 
    
    
    
###############################   
# Calculations module testing #
###############################  
class TestCalculations(object):

  # calculating speciation with empty system should throw
  def test_calcZero(self):
    with pytest.raises(RuntimeError):
      calcSpec(emptySys)  

  # calculation of system with 1 species, tied to H+
  def test_calcSimple(self):
    cs = ChemSys()
    cs.registerRxn("(T)aq = (H)+ , 0")
    (_, out, _) = calcSpec(cs)
    expect = [[10**(-x/100) for x in range(0, 1401)]]
    assert out == expect
    

##################################   
# Chemical System module testing #
##################################
class TestChemSys(object):

  # register empty string -> ValueError
  def test_register_rxn_bad1(self):
    cs = ChemSys()
    with pytest.raises(ValueError):
      cs.registerRxn("")  

  # register missing logK -> ValueError
  def test_register_rxn_bad2(self):
    cs = ChemSys()
    with pytest.raises(ValueError):
      cs.registerRxn("(H)+")
      
  # register missing = -> ValueError
  def test_register_rxn_bad3(self):
    cs = ChemSys()
    with pytest.raises(ValueError):
      cs.registerRxn("(H)+ , 0")
 
  # register with ill-formed state -> RuntimeError
  def test_register_rxn_bad4(self):
    cs = ChemSys()
    with pytest.raises(RuntimeError):
      cs.registerRxn("(H)bad = (H)bad , 0") 

  # register with ill-formed formula (non-alpha symbol) -> RuntimeError
  def test_register_rxn_bad5(self):
    cs = ChemSys()
    with pytest.raises(RuntimeError):
      cs.registerRxn("(1)+ = (1)+ , 0")  

  # register with ill-formed formula (beginning with lower case) -> RuntimeError
  def test_register_rxn_bad6(self):
    cs = ChemSys()
    with pytest.raises(RuntimeError):
      cs.registerRxn("(h)+ = (h)+ , 0")       
      
  # register with ill-formed formula (no parentheses) -> RuntimeError
  def test_register_rxn_bad7(self):
    cs = ChemSys()
    with pytest.raises(RuntimeError):
      cs.registerRxn("H+ = H+ , 0")       
      
  # register with unbalanced parentheses -> RuntimeError
  def test_register_rxn_bad8(self):
    cs = ChemSys()
    with pytest.raises(RuntimeError):
      cs.registerRxn("((H)+ = (H)+ , 0")       

  # register with superfluous -> should not throw
  def test_register_rxn_superfluous(self):
    cs = ChemSys()
    cs.registerRxn("((H))+ = (H)+ , 0") 
    assert True
    
  # register with a lot of parentheses -> should not throw
  def test_register_rxn_parens(self):
    cs = ChemSys()
    cs.registerRxn("((H(OH)2(H)))l = (H)+ , 0") 
    assert True
    
  # register total with neg. conc -> RuntimeError
  def test_register_total_neg(self):
    cs = ChemSys()
    with pytest.raises(RuntimeError):
      cs.registerTotal("H", -1) 
      
  # register total with zero conc -> RuntimeError
  def test_register_total_zero(self):
    cs = ChemSys()
    with pytest.raises(RuntimeError):
      cs.registerTotal("H", 0) 
    
  # register total with pos. conc -> should not throw
  def test_register_total_pos(self):
    cs = ChemSys()
    cs.registerTotal("H", 1)

  # generate diagram for fe system, test against known
  def test_gen_fe(self):
    feSys.specGen("fe")
    assert filecmp.cmp("fe.png", "./testref/fe_syst.png")
    os.remove("fe.png")
    
  # generate diagram for co3 system, test against known
  def test_gen_co3(self):
    co3Sys.specGen("co3")
    assert filecmp.cmp("co3.png", "./testref/co3_syst.png")
    os.remove("co3.png")