import pandas as pd
import numpy as np
import pyDOE2
from emulator_model import predict_chrom_phase

def generate_doe(num_exp: int, var_lims: dict, num_center_points=1, seed=123):
    rng = np.random.default_rng(seed)
    num_vars = len(var_lims)
    num_center_points = 1
    num_samples =num_exp-num_center_points
    # determine which vars are part of DOE
    doe_var = [1 if type(v) is list and v[0]<v[1] else 0 for v in var_lims.values()]
    doe_var_idx = np.cumsum(doe_var)-1

    # sample points in the latin hypercube
    doe_plan = pyDOE2.lhs(sum(doe_var),samples=num_samples,criterion="c")

    # fill remaining unscaled vars
    doe_unscaled = np.ones([num_exp,num_vars])*0.5
    for i in range(num_vars):
      if doe_var[i] == 1:
        doe_unscaled[num_center_points:,i] = doe_plan[:,doe_var_idx[i]]

    # scale all vars according to var_lims
    doe_scaled = doe_unscaled
    for i,k in enumerate(var_lims.keys()):
      if doe_var[i] ==1:
        doe_scaled[:,i] = doe_unscaled[:,i] * (var_lims[k][1]-var_lims[k][0]) + var_lims[k][0] 
      
      
    return doe_scaled


def generate_data(var_lims, num_runs, filename="mytable_data.csv"):
    
    num_center_points = 1
    model_param_combinations = generate_doe(num_runs, var_lims, num_center_points)
    doe_design = pd.DataFrame(model_param_combinations, columns = [k for k in var_lims.keys()])

    col_names = ["timestamps", "X:VCD", "X:Glc", "X:Lac", "X:Titer"]
    owu_df = pd.DataFrame(columns=col_names)

    i=0
    for i in range(num_runs):
        mu_g_max = model_param_combinations[i, 0]
        mu_d_max = model_param_combinations[i, 1]
        K_g_Glc = model_param_combinations[i, 2]
        K_g_Lac = model_param_combinations[i, 3]
        K_I_Lac = model_param_combinations[i, 4]
        K_d_Lac = model_param_combinations[i, 5]
        k_Glc = model_param_combinations[i, 6]
        k_Lac = model_param_combinations[i, 7]
        k_Prod = model_param_combinations[i, 8]
        feed_start = model_param_combinations[i, 9]
        feed_end = model_param_combinations[i, 10]
        Glc_feed_rate = model_param_combinations[i, 11]
        Glc_0 = model_param_combinations[i, 12]
        VCD_0 = model_param_combinations[i, 13]


        model_param = (mu_g_max, mu_d_max, K_g_Glc, K_I_Lac, K_d_Lac, k_Glc, k_Lac, k_Prod)
        process_param = (feed_start, feed_end, Glc_feed_rate, Glc_0, VCD_0)

        t, y = predict_chrom_phase(model_param, process_param)

        t = np.array([t])/24
        res = np.hstack([t.T, y])

        owu_df = pd.concat([owu_df, pd.DataFrame(res, columns=col_names)], ignore_index=True)


    owu_df.to_csv(filename, index=False)
    return owu_df
