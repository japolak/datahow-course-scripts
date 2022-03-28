function [conc_profile,BWU,OWU] = run_experiment(mdl_param,proc_param,varargin)
% Function to run an in-silico experiment
%
%   Call examples:
%   (see Process_Characterization.mlx)
%
% Mandatory input variables
%========================================
% - mdl_param: structure variable with the parameters of the in-silico model
%   + mdl_param.mu_g_max = maximum growth rate
%   + mdl_param.mu_d_max = maximum death rate
%   + mdl_param.K_g_Glc  = Michaelis-Menten (MM) coefficient for glucose
%   + mdl_param.K_I_Lac  = MM coefficient for lactate (inhibition)
%   + mdl_param.K_d_Lac  = MM coefficient for lactate (death)
%   + mdl_param.k_Glc    = specific glucose uptake rate
%   + mdl_param.k_Lac    = specific lactate production rate
%   + mdl_param.k_Prod   = specific product production rate
%   + mdl_param.k_Aggr   = product aggregation rate
% - proc_param: table with process parameters 
%   + proc_param = table(Glc_start,feed_start,feed_end,feed_rate,VCD_start, ...
%     'VariableNames',{'GlcConc0','GlcFeedStart','GlcFeedEnd','GlcFeedRate','VCD0'});
%   where:
%   + GlcConc0 = initial glucose concentration (t = 0 days)
%   + GlcFeedStart = start day for glucose feed
%   + GlcFeedEnd = end day for glucose feed
%   + GlcFeedRate = glucose feed rate
%   + VCD0 = initial VCD (t = 0 days)
%
% Optional input variables
%========================================
% - inputs_std_error: 2x1 vector of absolute errors for the input values of
%   VCD and glucose. The initial (input) value of VCD and Glc, respectively,
%   as well as the glucose feed rate, are randomized using these two values
%   as measured of the error standard deviation. If not supplied, no
%   randomization is applied.
% - outputs_std_error: 5x1 vector of absolute errors for the output values
%   of VCD, Glc, Lac, Titer and aggregates, respectively. The output values
%   of these quantities are randomized using these 5 values as measured of
%   the error standard deviation. If not supplied, no randomization is
%   applied.
%
% Output variables
%========================================
% - conc_profile: 15x6 array with the concentrations at each time (from day
%   zero to day 14). The first column is time (days). The other columns are
%   VCD, Glucose, Lactate, titer, aggregates. Since aggregates are measured
%   only at the end, only the last value is available (at day 14)
% - BWU: corresponding BWU row vector (aggregates are excluded)
% - OWU: corresponding OWU matrix (aggregates are excluded)
%

%% Input parser

p = inputParser;
addRequired(p,'mdl_param');
addRequired(p,'proc_param');
addOptional(p,'inputs_std_error',zeros(2,1))
addOptional(p,'outputs_std_error',zeros(5,1))
parse(p,mdl_param,proc_param,varargin{:});

%% Run experiment

% define feed schedule
Glc_feed_schedule = zeros(14,2,"double");
Glc_feed_schedule(:,1) = 0:13;
Glc_feed_schedule(proc_param.GlcFeedStart+1:proc_param.GlcFeedEnd,2) = ...
    proc_param.GlcFeedRate;

% run experiment
[conc_profile,BWU,OWU] = run_general_experiment(mdl_param, ...
    Glc_feed_schedule,proc_param.VCD0,proc_param.GlcConc0, ...
    p.Results.inputs_std_error,p.Results.outputs_std_error);

end
