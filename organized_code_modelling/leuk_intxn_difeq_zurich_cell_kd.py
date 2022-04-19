from pysb import *                 # http://docs.pysb.org/en/latest/tutorial.html
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.core import *
from sympy import sympify
from pysb.bng import *
from pysb.integrate import odesolve
from pysb.bng import run_ssa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

kinmod_tablebound = pd.read_csv("kineticmodel_v8_zurich_resting_leukmanual_tablebound.csv")

'''
# filter out rows on kinmod_tablebound that want to simulate blocking for perturbation tests
mut_gene = "APLP2" #['APLP2', 'CD209', 'CD320', 'CD58', 'CD80', 'CDH2', 'CHL1', 'CNTN1', 'EFNB1', 'ICAM1', 'IFNGR1', 'IL1RAP', 'JAG1', 'L1CAM', 'MCAM', 'NRCAM', 'OLR1', 'PDL1', 'PI16', 'PTPRF', 'SEMA4D', 'SEMA7A', 'SIRPA', 'SLITRK4', 'TNR16', 'TNR21', 'VASN', 'VISTA']
kinmod_tablebound = kinmod_tablebound[~kinmod_tablebound['Pair_label'].str.contains(" "+mut_gene+"$|^"+mut_gene+" ")]  # regex match to avoid "PVR" matching all genes includnig that string like pvrl1
'''

kinmod_cellpairs = kinmod_tablebound.drop(['Interaction', 'Pair_label'], axis = 1).sum()
new_indx = kinmod_cellpairs.index.str.replace(" \\+ ", "_")
new_indx = pd.Index(["k_"]*len(new_indx)).str.cat(new_indx)

kinmod_cellpairs = kinmod_cellpairs.rename(dict(zip(kinmod_cellpairs.index, new_indx)))
kdarg_list = ['k_T4_T4','k_T4_T8','k_T4_NK','k_T4_DC','k_T4_M14',
              'k_T4_M16','k_T4_B','k_T8_T8','k_T8_NK','k_T8_DC',
              'k_T8_M14','k_T8_M16','k_T8_B','k_NK_NK','k_NK_DC',
              'k_NK_M14','k_NK_M16','k_NK_B','k_DC_DC','k_DC_M14',
              'k_DC_M16','k_DC_B','k_M14_M14','k_M14_M16','k_M14_B',
              'k_M16_M16','k_M16_B','k_B_B']
kinmod_cellpairs = kinmod_cellpairs.filter(items = kdarg_list, axis = 'index')

kinmod_celldict = dict(zip(kinmod_cellpairs.index, (1000 / kinmod_cellpairs).values)) 



Model()                            # start of definition for model

Monomer('T4', ['b'])
Monomer('T8', ['b'])
Monomer('NK', ['b'])
Monomer('DC', ['b'])
Monomer('M14', ['b'])
Monomer('M16', ['b'])
Monomer('B', ['b'])

Initial(T4(b = None), Parameter('T4_0', 0.2*0.65)) # initialize at experimentally-observed cell type abundance values
Initial(T8(b = None), Parameter('T8_0', 0.2*0.35))
Initial(NK(b = None), Parameter('NK_0', 0.06))
Initial(DC(b = None), Parameter('DC_0', 0.006))
Initial(M14(b = None), Parameter('M14_0', 0.06))
Initial(M16(b = None), Parameter('M16_0', 0.01))
Initial(B(b = None), Parameter('B_0', 0.05))

for kd_arg in kinmod_celldict.keys():
    Parameter(kd_arg, kinmod_celldict[kd_arg])
#end

Parameter('k_bump', 1)  # only fixing the K_off value of the overall Kd, so assumes relative association rate of cells bumping into each other is constant

Rule('bind_T4_T4', T4(b = None) + T4(b = None) | T4(b = 1) % T4(b = 1), k_bump, k_T4_T4)
Rule('bind_T4_T8', T4(b = None) + T8(b = None) | T4(b = 1) % T8(b = 1), k_bump, k_T4_T8)
Rule('bind_T4_NK', T4(b = None) + NK(b = None) | T4(b = 1) % NK(b = 1), k_bump, k_T4_NK)
Rule('bind_T4_DC', T4(b = None) + DC(b = None) | T4(b = 1) % DC(b = 1), k_bump, k_T4_DC)
Rule('bind_T4_M14', T4(b = None) + M14(b = None) | T4(b = 1) % M14(b = 1), k_bump, k_T4_M14)
Rule('bind_T4_M16', T4(b = None) + M16(b = None) | T4(b = 1) % M16(b = 1), k_bump, k_T4_M16)
Rule('bind_T4_B', T4(b = None) + B(b = None) | T4(b = 1) % B(b = 1), k_bump, k_T4_B)
Rule('bind_T8_T8', T8(b = None) + T8(b = None) | T8(b = 1) % T8(b = 1), k_bump, k_T8_T8)
Rule('bind_T8_NK', T8(b = None) + NK(b = None) | T8(b = 1) % NK(b = 1), k_bump, k_T8_NK)
Rule('bind_T8_DC', T8(b = None) + DC(b = None) | T8(b = 1) % DC(b = 1), k_bump, k_T8_DC)
Rule('bind_T8_M14', T8(b = None) + M14(b = None) | T8(b = 1) % M14(b = 1), k_bump, k_T8_M14)
Rule('bind_T8_M16', T8(b = None) + M16(b = None) | T8(b = 1) % M16(b = 1), k_bump, k_T8_M16)
Rule('bind_T8_B', T8(b = None) + B(b = None) | T8(b = 1) % B(b = 1), k_bump, k_T8_B)
Rule('bind_NK_NK', NK(b = None) + NK(b = None) | NK(b = 1) % NK(b = 1), k_bump, k_NK_NK)
Rule('bind_NK_DC', NK(b = None) + DC(b = None) | NK(b = 1) % DC(b = 1), k_bump, k_NK_DC)
Rule('bind_NK_M14', NK(b = None) + M14(b = None) | NK(b = 1) % M14(b = 1), k_bump, k_NK_M14)
Rule('bind_NK_M16', NK(b = None) + M16(b = None) | NK(b = 1) % M16(b = 1), k_bump, k_NK_M16)
Rule('bind_NK_B', NK(b = None) + B(b = None) | NK(b = 1) % B(b = 1), k_bump, k_NK_B)
Rule('bind_DC_DC', DC(b = None) + DC(b = None) | DC(b = 1) % DC(b = 1), k_bump, k_DC_DC)
Rule('bind_DC_M14', DC(b = None) + M14(b = None) | DC(b = 1) % M14(b = 1), k_bump, k_DC_M14)
Rule('bind_DC_M16', DC(b = None) + M16(b = None) | DC(b = 1) % M16(b = 1), k_bump, k_DC_M16)
Rule('bind_DC_B', DC(b = None) + B(b = None) | DC(b = 1) % B(b = 1), k_bump, k_DC_B)
Rule('bind_M14_M14', M14(b = None) + M14(b = None) | M14(b = 1) % M14(b = 1), k_bump, k_M14_M14)
Rule('bind_M14_M16', M14(b = None) + M16(b = None) | M14(b = 1) % M16(b = 1), k_bump, k_M14_M16)
Rule('bind_M14_B', M14(b = None) + B(b = None) | M14(b = 1) % B(b = 1), k_bump, k_M14_B)
Rule('bind_M16_M16', M16(b = None) + M16(b = None) | M16(b = 1) % M16(b = 1), k_bump, k_M16_M16)
Rule('bind_M16_B', M16(b = None) + B(b = None) | M16(b = 1) % B(b = 1), k_bump, k_M16_B)
Rule('bind_B_B', B(b = None) + B(b = None) | B(b = 1) % B(b = 1), k_bump, k_B_B)

Observable('oT4', T4(b=None))
Observable('oT8', T8(b=None))
Observable('oNK', NK(b=None))
Observable('oDC', DC(b=None))
Observable('oM14', M14(b=None))
Observable('oM16', M16(b=None))
Observable('oB', B(b=None))
Observable('oT4_T4', T4(b = 1) % T4(b = 1))
Observable('oT4_T8', T4(b = 1) % T8(b = 1))
Observable('oT4_NK', T4(b = 1) % NK(b = 1))
Observable('oT4_DC', T4(b = 1) % DC(b = 1))
Observable('oT4_M14', T4(b = 1) % M14(b = 1))
Observable('oT4_M16', T4(b = 1) % M16(b = 1))
Observable('oT4_B', T4(b = 1) % B(b = 1))
Observable('oT8_T8', T8(b = 1) % T8(b = 1))
Observable('oT8_NK', T8(b = 1) % NK(b = 1))   
Observable('oT8_DC', T8(b = 1) % DC(b = 1))
Observable('oT8_M14', T8(b = 1) % M14(b = 1))
Observable('oT8_M16', T8(b = 1) % M16(b = 1))
Observable('oT8_B', T8(b = 1) % B(b = 1))
Observable('oNK_NK', NK(b = 1) % NK(b = 1))
Observable('oNK_DC', NK(b = 1) % DC(b = 1))
Observable('oNK_M14', NK(b = 1) % M14(b = 1))
Observable('oNK_M16', NK(b = 1) % M16(b = 1))
Observable('oNK_B', NK(b = 1) % B(b = 1))
Observable('oDC_DC', DC(b = 1) % DC(b = 1))
Observable('oDC_M14', DC(b = 1) % M14(b = 1))
Observable('oDC_M16', DC(b = 1) % M16(b = 1))
Observable('oDC_B', DC(b = 1) % B(b = 1))
Observable('oM14_M14', M14(b = 1) % M14(b = 1))
Observable('oM14_M16', M14(b = 1) % M16(b = 1))
Observable('oM14_B', M14(b = 1) % B(b = 1))
Observable('oM16_M16', M16(b = 1) % M16(b = 1))
Observable('oM16_B', M16(b = 1) % B(b = 1))
Observable('oB_B', B(b = 1) % B(b = 1))

#python -m pysb.export .\leuk_intxn_difeq_zurich_cell_kd.py mathematica > zurich_difeq_mathematica.txt