function [conc_profile,BWU,OWU] = run_general_experiment(mdl_param, ...
    Glc_feed_schedule,VCD_start,Glc_start,inputs_std_err,outputs_std_err)
% Function to run an in-silico experiment with a general Glc feed schedule

% define time steps for integration (transform all in hours)
time_steps = 24*(0:0.1:14);
Glc_feed_schedule(:,2) = Glc_feed_schedule(:,2)+inputs_std_err(2)*randn(1);
Glc_feed_schedule(:,1) = Glc_feed_schedule(:,1)*24;
Glc_feed_schedule(:,2) = Glc_feed_schedule(:,2)/24;

% define initial conditions
X0 = zeros(5,1);
X0(1) = VCD_start+inputs_std_err(1)*randn(1);
X0(2) = Glc_start+inputs_std_err(2)*randn(1);

% define integration options
global integr_fcn_calls
integr_fcn_calls = 0;
max_integration_calls = 500000;
ode_prime_output_an = @(t,X,flag) ...
    ode_prime_output(t,X,flag,max_integration_calls);
options = odeset('RelTol',1e-4,'AbsTol',1e-4,'MaxStep',0.01, ...
    'OutputFcn',ode_prime_output_an,'BDF','on');
warning('off')

% run step integration
ode_prime_an = @(t,X) ode_prime(t,X,Glc_feed_schedule,mdl_param);
try
    [t_out,X_out] = ode15s(ode_prime_an,time_steps,X0,options);
catch ME
    throw(ME)
end

% define output
if size(X_out,1) ~= length(time_steps)
    proc_param
    error('Integration halted')
end
t_out = t_out/24;
v = find(mod(t_out,1) == 0);
conc_profile = [t_out(v),X_out(v,:)];

% apply output errors
outputs_std_err = outputs_std_err(:)';
abs_err = repmat(outputs_std_err,14,1).*randn(14,5);
conc_profile(2:end,2:6) = conc_profile(2:end,2:6)+abs_err;

% define BWU vector
if nargout >= 2
    v = find(mod(conc_profile(:,1),1) == 0);
    BWU = zeros(1,75,'double');
    for idx = 1:14
        n = (idx-1)*5;
        BWU(n+1:n+4) = conc_profile(v(idx),2:5);
        BWU(n+5) = Glc_feed_schedule(idx,2)*24;
    end
    BWU(71:74) = conc_profile(v(15),2:5);
end

% define OWU vector
if nargout >= 3
    OWU = zeros(15,5,'double');
    OWU(:,1:4) = conc_profile(:,2:5);
    OWU(1:14,5) = Glc_feed_schedule(:,2)*24;
end

end