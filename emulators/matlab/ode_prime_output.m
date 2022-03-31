function status = ode_prime_output(~,~,~,max_integration_calls)
% function to control process of integration

global integr_fcn_calls

status = 0;
if integr_fcn_calls > max_integration_calls
    status = 1;
end

end
