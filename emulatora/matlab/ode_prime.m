function dXdt = ode_prime(t,X,Glc_feed_schedule,mdl_param)
% This function integrates the mass balance

% update integration function call counter
global integr_fcn_calls
integr_fcn_calls = integr_fcn_calls+1;

% get concentrations
VCD = X(1);
Glc = X(2);
Lac = X(3);
Titer = X(4);

% initialize vector of derivatives
% X(1) -> Viable Cell Density (VCD [1e6 cells/ml])
% X(2) -> Glucose (Glc [mmol])
% X(3) -> Lactate (Lac [mmol])
% X(4) -> Titer (Product concentration [g/L])
% X(5) -> Aggregates [g/L]
dXdt = zeros(5,1);

% compute growth and death rates
MM_Glc = Glc/(mdl_param.K_g_Glc+Glc);
mu_g = mdl_param.mu_g_max*MM_Glc* ...
       mdl_param.K_I_Lac/(mdl_param.K_I_Lac+Lac);
phi = exp(0.1*(Glc-75));
mu_d = mdl_param.mu_d_max*(1+phi/(1+phi))* ...
       Lac/(mdl_param.K_d_Lac+Lac);
growth_ratio = mu_g/mdl_param.mu_g_max;
   
% define Glc feed rate
if t == 0
    n = 1;
else
    n = find(Glc_feed_schedule(:,1) < t,1,'last');
end
Glc_feed_rate = Glc_feed_schedule(n,2);

% compute mass balances
Glc_Min = Glc/(0.05+Glc);
dXdt(1) = (mu_g-mu_d)*VCD;
dXdt(2) = -mdl_param.k_Glc*Glc_Min*VCD+Glc_feed_rate;
dXdt(3) = mdl_param.k_Lac*VCD;
dXdt(5) = mdl_param.k_Aggr*Titer^2;
dXdt(4) = mdl_param.k_Prod*MM_Glc*(1-growth_ratio)^2*VCD-2*dXdt(5);

end
