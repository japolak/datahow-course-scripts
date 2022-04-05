import pandas as pd
import numpy as np
import scipy.stats
import scipy.integrate


def chrome_ode(t, y, p):
    mu_g_max, mu_d_max, K_g_Glc, K_I_Lac, K_d_Lac, k_Glc, k_Lac, k_Prod, feed_start, feed_end, Glc_feed_rate = p

    VCD, Glc, Lac, titer  = y[0], y[1], y[2], y[3]
    MM_Glc = Glc / (K_g_Glc + Glc)
    mu_g = mu_g_max * MM_Glc * K_I_Lac / (K_I_Lac + Lac)
    phi = np.exp(0.1 * (Glc - 75.0))
    mu_d = mu_d_max * (1.0 + phi / (1.0 + phi)) * Lac / (K_d_Lac + Lac)
    growth_ratio = mu_g / mu_g_max

    # compute mass balances
    Glc_Min = Glc / (0.05 + Glc)
    dVCDdt = (mu_g - mu_d) * VCD
    dGlcdt = -k_Glc * Glc_Min * VCD
    dLacdt = k_Lac * VCD
    dTiterdt = k_Prod * MM_Glc * ((1.0 - growth_ratio) ** 2.0) * VCD

    # add feed rate
    if feed_end >= t >= feed_start:
        dGlcdt += Glc_feed_rate

    return [dVCDdt, dGlcdt, dLacdt, dTiterdt]


def generate_doe(num_exp: int, var_lims: dict, num_center_points: int, seed=123):
    rng = np.random.default_rng(seed)
    num_vars = len(var_lims)

    sampler = scipy.stats.qmc.LatinHypercube(d=num_vars, centered=True, seed=rng)

    doe_plan = sampler.random(n=num_exp - num_center_points)
    doe_plan = np.vstack([doe_plan, 0.5 * np.ones((num_center_points, num_vars))])

    # scales columns
    doe_plan[:, 0] = doe_plan[:, 0] * (var_lims["mu_g_max"][1] - var_lims["mu_g_max"][0]) + var_lims["mu_g_max"][0]
    doe_plan[:, 1] = doe_plan[:, 1] * (var_lims["mu_d_max"][1] - var_lims["mu_d_max"][0]) + var_lims["mu_d_max"][0]
    doe_plan[:, 2] = doe_plan[:, 2] * (var_lims["K_g_Glc"][1] - var_lims["K_g_Glc"][0]) + var_lims["K_g_Glc"][0]
    doe_plan[:, 3] = doe_plan[:, 3] * (var_lims["K_g_Lac"][1] - var_lims["K_g_Lac"][0]) + var_lims["K_g_Lac"][0]
    doe_plan[:, 4] = doe_plan[:, 4] * (var_lims["K_I_Lac"][1] - var_lims["K_I_Lac"][0]) + var_lims["K_I_Lac"][0]
    doe_plan[:, 5] = doe_plan[:, 5] * (var_lims["K_d_Lac"][1] - var_lims["K_d_Lac"][0]) + var_lims["K_d_Lac"][0]
    doe_plan[:, 6] = doe_plan[:, 6] * (var_lims["k_Glc"][1] - var_lims["k_Glc"][0]) + var_lims["k_Glc"][0]
    doe_plan[:, 7] = doe_plan[:, 7] * (var_lims["k_Lac"][1] - var_lims["k_Lac"][0]) + var_lims["k_Lac"][0]
    doe_plan[:, 8] = doe_plan[:, 8] * (var_lims["k_Prod"][1] - var_lims["k_Prod"][0]) + var_lims["k_Prod"][0]
    doe_plan[:, 9] = doe_plan[:, 9] * (var_lims["feed_start"][1] - var_lims["feed_start"][0]) + var_lims["feed_start"][0]
    doe_plan[:, 10] = doe_plan[:, 10] * (var_lims["feed_end"][1] - var_lims["feed_end"][0]) + var_lims["feed_end"][0]
    doe_plan[:, 11] = doe_plan[:, 11] * (var_lims["Glc_feed_rate"][1] - var_lims["Glc_feed_rate"][0]) + var_lims["Glc_feed_rate"][0]
    doe_plan[:, 12] = doe_plan[:, 12] * (var_lims["Glc_0"][1] - var_lims["Glc_0"][0]) + var_lims["Glc_0"][0]
    doe_plan[:, 13] = doe_plan[:, 13] * (var_lims["VCD_0"][1] - var_lims["VCD_0"][0]) + var_lims["VCD_0"][0]

    return doe_plan


def predict_chrom_phase(model_param, feed_start, feed_end, Glc_feed_rate, Glc_0, VCD_0):
    y0 = [VCD_0, Glc_0, 0, 0]
    t_start, t_end = 0, 24 * 14
    t_span = np.arange(t_start, t_end + 24, 24)

    mu_g_max, mu_d_max, K_g_Glc, K_I_Lac, K_d_Lac, k_Glc, k_Lac, k_Prod = model_param
    p = (mu_g_max, mu_d_max, K_g_Glc, K_I_Lac, K_d_Lac, k_Glc, k_Lac, k_Prod, 24.0 * feed_start, 24.0 * feed_end, Glc_feed_rate)

    # integrates equations
    sol = scipy.integrate.solve_ivp(
        chrome_ode,
        t_span=[t_start, t_end],
        y0=y0,
        t_eval=t_span,
        method="BDF",
        args=(
            [p]
        ),
        rtol=1e-6,
        atol=1e-6,
    )

    t = sol.t.tolist()
    y = sol.y.T.tolist()

    return t, y


def generate_data(var_lims, amm_runs, filename="mytable_data.csv"):
    num_center_points = 2
    model_param_combinations = generate_doe(amm_runs, var_lims, num_center_points)

    col_names = ["timestamps", "X:VCD", "X:Glc", "X:Lac", "X:Titer"]
    final_df = pd.DataFrame(columns=col_names)

    for i in range(amm_runs):
        mu_g_max = model_param_combinations[i, 0]
        mu_d_max = model_param_combinations[i, 1]
        K_g_Glc = model_param_combinations[i, 2]
        K_I_Lac = model_param_combinations[i, 3]
        K_d_Lac = model_param_combinations[i, 4]
        k_Glc = model_param_combinations[i, 5]
        k_Lac = model_param_combinations[i, 6]
        k_Prod = model_param_combinations[i, 7]
        feed_start = model_param_combinations[i, 8]
        feed_end = model_param_combinations[i, 9]
        Glc_feed_rate = model_param_combinations[i, 10]
        Glc_0 = model_param_combinations[i, 11]
        VCD_0 = model_param_combinations[i, 12]

        model_param = (mu_g_max, mu_d_max, K_g_Glc, K_I_Lac, K_d_Lac, k_Glc, k_Lac, k_Prod)

        t, y = predict_chrom_phase(model_param, feed_start, feed_end, Glc_feed_rate, Glc_0, VCD_0)
        t = np.array([t])
        res = np.hstack([t.T, y])

        final_df = pd.concat([final_df, pd.DataFrame(res, columns=col_names)], ignore_index=True)

    final_df.to_csv(filename, index=False)
    return final_df


if __name__ == '__main__':
    var_lims = {"mu_g_max": [0, 1],
                "mu_d_max": [0, 1],
                "K_g_Glc": [0, 1],
                "K_g_Lac": [0, 1],
                "K_I_Lac": [0, 1],
                "K_d_Lac": [0, 1],
                "k_Glc": [0, 1],
                "k_Lac": [0, 1],
                "k_Prod": [0, 1],
                "feed_start": [0, 1],
                "feed_end": [4, 5],
                "Glc_feed_rate": [0, 1],
                "Glc_0": [1, 2],
                "VCD_0": [1, 2]}

    df = generate_data(var_lims, 5)
    print(df)