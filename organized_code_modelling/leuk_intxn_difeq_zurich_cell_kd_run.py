
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import pylab as pl
from os import chdir
chdir("./organized_code_modelling")

import leuk_intxn_difeq_zurich_cell_kd as m

t = pl.linspace(0, 10)
simres = ScipyOdeSimulator(m.model, tspan=t).run()
yout = simres.all

np.savetxt("./difeq_outfiles/difeq_out.txt", yout, 
            header = " ".join([tup[0] for tup in yout.dtype.fields.items()]))


#### simple plot to check output ####

import matplotlib.pyplot as plt
if_plot_pairs_only = False

plt.figure()
if not(if_plot_pairs_only):
    plt.plot(t, yout['oT4'], label = "T4")
    plt.plot(t, yout['oT8'], label = "T8")
    plt.plot(t, yout['oNK'], label = "NK")
    plt.plot(t, yout['oDC'], label = "DC")
    plt.plot(t, yout['oM14'], label = "M14")
    plt.plot(t, yout['oM16'], label = "M16")
    plt.plot(t, yout['oB'], label = "B")

plt.plot(t, yout['oT4_T4'], label = "T4_T4")
plt.plot(t, yout['oT4_T8'], label = "T4_T8")
plt.plot(t, yout['oT4_NK'], label = "T4_NK")
plt.plot(t, yout['oT4_DC'], label = "T4_DC")
plt.plot(t, yout['oT4_M14'], label = "T4_M14")
plt.plot(t, yout['oT4_M16'], label = "T4_M16")
plt.plot(t, yout['oT4_B'], label = "T4_B")
if if_plot_pairs_only:
    plt.plot(t, yout['oT8_T8'], label = "T8_T8")
    plt.plot(t, yout['oT8_NK'], label = "T8_NK")
    plt.plot(t, yout['oT8_DC'], label = "T8_DC")
    plt.plot(t, yout['oT8_M14'], label = "T8_M14")
    plt.plot(t, yout['oT8_M16'], label = "T8_M16")
    plt.plot(t, yout['oT8_B'], label = "T8_B")
    plt.plot(t, yout['oNK_NK'], label = "NK_NK")
    plt.plot(t, yout['oNK_DC'], label = "NK_DC")
    plt.plot(t, yout['oNK_M14'], label = "NK_M14")
    plt.plot(t, yout['oNK_M16'], label = "NK_M16")
    plt.plot(t, yout['oNK_B'], label = "NK_B")
    plt.plot(t, yout['oDC_DC'], label = "DC_DC")
    plt.plot(t, yout['oDC_M14'], label = "DC_M14")
    plt.plot(t, yout['oDC_M16'], label = "DC_M16")
    plt.plot(t, yout['oDC_B'], label = "DC_B")
    plt.plot(t, yout['oM14_M14'], label = "M14_M14")
    plt.plot(t, yout['oM14_M16'], label = "M14_M16")
    plt.plot(t, yout['oM14_B'], label = "M14_B")
    plt.plot(t, yout['oM16_M16'], label = "M16_M16")
    plt.plot(t, yout['oM16_B'], label = "M16_B")
    plt.plot(t, yout['oB_B'], label = "B_B")

plt.legend()
plt.xlabel("Time (s)")
plt.ylabel("Cell type proportion")
plt.show()
