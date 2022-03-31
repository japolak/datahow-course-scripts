# declare modules
using DifferentialEquations
using Plots
using DataFrames
using CSV
using LatinHypercubeSampling

function generate_data(feed_start, feed_end, Glc_feed_rate, Glc_0, VCD_0, mu_g_max, mu_d_max, K_g_Glc, K_I_Lac, K_d_Lac, k_Glc, k_Lac, k_Prod)

    # function to compute ODE derivatives
    function chrom_ode!(dy, y, p, t)

        mu_g_max, mu_d_max, K_g_Glc, K_I_Lac, K_d_Lac, k_Glc, k_Lac, k_Prod, feed_start, feed_end, Glc_feed_rate = p

        VCD = @view y[1:1]
        Glc = @view y[2:2]
        Lac = @view y[3:3]
        dVCDdt = @view dy[1:1]
        dGlcdt = @view dy[2:2]
        dLacdt = @view dy[3:3]
        dTiterdt = @view dy[4:4]

        # compute growth and death rates
        MM_Glc = Glc[1] / (K_g_Glc + Glc[1])
        mu_g = mu_g_max * MM_Glc * K_I_Lac / (K_I_Lac + Lac[1])
        phi = exp(0.1 * (Glc[1] - 75.0))
        mu_d = mu_d_max * (1.0 + phi / (1.0 + phi)) * Lac[1] / (K_d_Lac + Lac[1])
        growth_ratio = mu_g / mu_g_max

        # compute mass balances
        Glc_Min = Glc[1] / (0.05 + Glc[1])
        dVCDdt[1] = (mu_g - mu_d) * VCD[1]
        dGlcdt[1] = -k_Glc * Glc_Min * VCD[1]
        dLacdt[1] = k_Lac * VCD[1]
        dTiterdt[1] = k_Prod * MM_Glc * ((1.0 - growth_ratio)^2.0) * VCD[1]

        # add feed rate
        if (t >= feed_start) && (t <= feed_end)
            dGlcdt[1] += Glc_feed_rate
        end

        return nothing
    end

    function generate_doe(num_exp, var_lims, num_center_points)


        # define the number of experiments in the latin hypercube
        Nlhc = num_exp - num_center_points
        num_vars = length(keys(var_lims))

        doe_plan, _ = LHCoptim(Nlhc, num_vars, 100)

        scaled_doe_plan = scaleLHC(doe_plan, [
            (var_lims["mu_g_max"][1], var_lims["mu_g_max"][2]),
            (var_lims["mu_d_max"][1], var_lims["mu_d_max"][2]),
            (var_lims["K_g_Glc"][1], var_lims["K_g_Glc"][2]),
            (var_lims["K_g_Lac"][1], var_lims["K_g_Lac"][2]),
            (var_lims["K_I_Lac"][1], var_lims["K_I_Lac"][2]),
            (var_lims["K_d_Lac"][1], var_lims["K_d_Lac"][2]),
            (var_lims["k_Glc"][1], var_lims["k_Glc"][2]),
            (var_lims["k_Lac"][1], var_lims["k_Lac"][2]),
            (var_lims["k_Prod"][1], var_lims["k_Prod"][2]),
            (var_lims["feed_start"][1], var_lims["feed_start"][2]),
            (var_lims["feed_end"][1], var_lims["feed_end"][2]),
            (var_lims["Glc_feed_rate"][1], var_lims["Glc_feed_rate"][2]),
            (var_lims["Glc_0"][1], var_lims["Glc_0"][2]),
            (var_lims["VCD_0"][1], var_lims["VCD_0"][2]),
        ])
        for i = 1:num_center_points
            c_point = [
                var_lims["mu_g_max"][1] * 0.5 + var_lims["mu_g_max"][2] * 0.5,
                var_lims["mu_d_max"][1] * 0.5 + var_lims["mu_d_max"][2] * 0.5,
                var_lims["K_g_Glc"][1] * 0.5 + var_lims["K_g_Glc"][2] * 0.5,
                var_lims["K_g_Lac"][1] * 0.5 + var_lims["K_g_Lac"][2] * 0.5,
                var_lims["K_I_Lac"][1] * 0.5 + var_lims["K_I_Lac"][2] * 0.5,
                var_lims["K_d_Lac"][1] * 0.5 + var_lims["K_d_Lac"][2] * 0.5,
                var_lims["k_Glc"][1] * 0.5 + var_lims["k_Glc"][2] * 0.5,
                var_lims["k_Lac"][1] * 0.5 + var_lims["k_Lac"][2] * 0.5,
                var_lims["k_Prod"][1] * 0.5 + var_lims["k_Prod"][2] * 0.5,
                var_lims["feed_start"][1] * 0.5 + var_lims["feed_start"][2] * 0.5,
                var_lims["feed_end"][1] * 0.5 + var_lims["feed_end"][2] * 0.5,
                var_lims["Glc_feed_rate"][1] * 0.5 + var_lims["Glc_feed_rate"][2] * 0.5,
                var_lims["Glc_0"][1] * 0.5 + var_lims["Glc_0"][2] * 0.5,
                var_lims["VCD_0"][1] * 0.5 + var_lims["VCD_0"][2] * 0.5
            ]
            scaled_doe_plan = [scaled_doe_plan; c_point']
        end
        print("breakpoint")
        return scaled_doe_plan

    end


    # function to run the model over a phase
    function predict_chrom_phase(model_param, feed_start, feed_end, Glc_feed_rate, Glc_0, VCD_0)

        y0 = zeros(Float64, 4)
        y0[1] = VCD_0
        y0[2] = Glc_0
        tau_span = (0.0, 24.0 * 14.0)
        tau_output = 24.0 .* range(0.0, 14.0, 15)

        mu_g_max, mu_d_max, K_g_Glc, K_I_Lac, K_d_Lac, k_Glc, k_Lac, k_Prod = model_param
        p = mu_g_max, mu_d_max, K_g_Glc, K_I_Lac, K_d_Lac, k_Glc, k_Lac, k_Prod, 24.0 * feed_start, 24.0 * feed_end, Glc_feed_rate

        prob = ODEProblem(chrom_ode!, y0, tau_span, p)

        sol = solve(prob, BS3(), reltol=1e-4, abstol=1e-4, saveat=tau_output)

        t = sol.t ./ 24.0
        y = Array(sol)

        return t, y
    end

    model_param = (mu_g_max, mu_d_max, K_g_Glc, K_I_Lac, K_d_Lac, k_Glc, k_Lac, k_Prod)

    i = 1
    col_names = ["timestamps", "X:VCD", "X:Glc", "X:Lac", "X:Titer"]

    local t, y = predict_chrom_phase(model_param, feed_start, feed_end, Glc_feed_rate, Glc_0, VCD_0)
    res = [t'; y]'
    final_df = DataFrame(res, col_names)

    efg = open("mytable_1_run.csv", "w")
    CSV.write("mytable_1_run.csv", final_df)


    return final_df
    # plot(t, y[1, :], vars = (0, 2))
end


# var_lims = Dict("mu_g_max" => [0.01, 0.10],
#     "mu_d_max" => [0.025, 0.075],
#     "K_g_Glc" => [1, 2],
#     "K_g_Lac" => [50, 80],
#     "K_I_Lac" => [30, 50],
#     "K_d_Lac" => [50, 80],
#     "k_Glc" => [0.04, 0.12],
#     "k_Lac" => [0.06, 0.15],
#     "k_Prod" => [1, 3],
#     "feed_start" => [3, 6],
#     "feed_end" => [12.0, 16],
#     "Glc_feed_rate" => [12.0, 15.0],
#     "Glc_0" => [12.0, 15.0],
#     "VCD_0" => [0.3, 0.6])

# generate_data(feed_start, feed_end, Glc_feed_rate, Glc_0, VCD_0, mu_g_max, mu_d_max, K_g_Glc, K_I_Lac, K_d_Lac, k_Glc, k_Lac, k_Prod, var_lims, amm_runs)
# generate_data(var_lims, amm_runs)