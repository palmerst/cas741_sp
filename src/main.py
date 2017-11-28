from ChemSys import *

sys = ChemSys()
sys.registerRxn("(Fe)3+ + (H2O)l = (Fe(OH))2+ + (H)+ , -2.19")
sys.registerRxn("(Fe)3+ + 2(H2O)l = (Fe(OH)2)+ + 2(H)+ , -5.67")
sys.registerRxn("(Fe)3+ + 4(H2O)l = (Fe(OH)4)- + 4(H)+ , -21.6")
sys.registerTotal("Fe", 0.000010)
sys.specGen("fe_syst")

sys2 = ChemSys()
sys2.registerRxn("(CO2)aq + (H2O)l = (H2CO3)aq , -2.77")
sys2.registerRxn("(H2CO3)aq = (HCO3)- + (H)+ , -3.6")
sys2.registerRxn("(HCO3)- = (CO3)2- + (H)+ , -10.33")
sys2.registerTotal("C", 0.000010)
sys2.specGen("co2_syst")