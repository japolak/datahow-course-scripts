# declare modules
using DifferentialEquations
using DataFrames
using CSV
using LatinHypercubeSampling

function generate_data(var_lims, amm_runs, filename="mytable_data.csv")

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

    num_center_points = 2
    model_param_combinations = generate_doe(amm_runs, var_lims, num_center_points)

    i = 1
    # amm_runs = 8
    col_names = ["timestamps", "X:VCD", "X:Glc", "X:Lac", "X:Titer"]
    for i = 1:amm_runs


        mu_g_max = model_param_combinations[i, 1]
        mu_d_max = model_param_combinations[i, 2]
        K_g_Glc = model_param_combinations[i, 3]
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

        local t, y = predict_chrom_phase(model_param, feed_start, feed_end, Glc_feed_rate, Glc_0, VCD_0)
        res = [t'; y]'
        new_df = DataFrame(res, col_names)

        if i == 1
            final_df = new_df
        else
            global final_df = vcat(final_df, new_df)
        end
        i += 1
    end
    efg = open(filename, "w")
    CSV.write(filename, final_df)

    return final_df
    
end


