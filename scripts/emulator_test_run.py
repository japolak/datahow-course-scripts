import pandas as pd
import numpy as np
import scipy
import pyDOE2
import importlib
import scipy.integrate
from modules import emulator
from modules.transformations import owu_from_csv

# mu_g_max = [0.01, 0.10]
# mu_d_max = [0.025, 0.075]
# K_g_Glc = [1, 2]
# K_g_Lac = [50, 80]
# K_I_Lac = [30, 50]
# K_d_Lac = [50, 80]
# k_Glc = [0.04, 0.12]
# k_Lac = [0.06, 0.15]
# k_Prod = [1, 3]

mu_g_max = 0.01
mu_d_max = 0.025
K_g_Glc = 1
K_g_Lac = 50
K_I_Lac = 30
K_d_Lac = 50
k_Glc = 0.04
k_Lac = 0.06
k_Prod = 1

# Process parameters
feed_start = [3, 6]
feed_end = [12.0, 16]
Glc_feed_rate = [0.50, 2.0]
Glc_0 = [12.0, 15.0]
VCD_0 = [0.3, 0.6]

var_lims = {
    "mu_g_max": mu_g_max,
    "mu_d_max": mu_d_max,
    "K_g_Glc": K_g_Glc,
    "K_g_Lac": K_g_Lac,
    "K_I_Lac": K_I_Lac,
    "K_d_Lac": K_d_Lac,
    "k_Glc": k_Glc,
    "k_Lac": k_Lac,
    "k_Prod": k_Prod,
    "feed_start": feed_start,
    "feed_end": feed_end,
    "Glc_feed_rate": Glc_feed_rate,
    "Glc_0": Glc_0,
    "VCD_0": VCD_0,
}

num_runs = 5
filename = "mytable.csv"
owu_df = emulator.generate_data(var_lims, num_runs, filename)
owu = owu_from_csv(filename)
print(owu)
