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


def predict_chrom_phase(model_param, process_param):
    feed_start, feed_end, Glc_feed_rate, Glc_0, VCD_0 = process_param
    mu_g_max, mu_d_max, K_g_Glc, K_I_Lac, K_d_Lac, k_Glc, k_Lac, k_Prod = model_param

    y0 = [VCD_0, Glc_0, 0, 0]
    t_start, t_end = 0, 24 * 14
    t_span = np.arange(t_start, t_end + 24, 24)
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
