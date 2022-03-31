function [DoE,f_DoE,BWU,OWU,Conc] = run_DoE(mdl_param,DoE_nondim, ...
    Glc_start_range,feed_start_range,feed_end_range,feed_rate_range, ...
    VCD_start_range,inputs_std_error,outputs_std_error,NCpus)

% define dimensional DoE
DoE = DoE_nondim;
DoE(:,1) = mean(Glc_start_range)+0.5*(Glc_start_range(2)-Glc_start_range(1))*DoE(:,1);
DoE(:,2) = round(mean(feed_start_range)+0.5*(feed_start_range(2)-feed_start_range(1))*DoE(:,2));
DoE(:,3) = round(mean(feed_end_range)+0.5*(feed_end_range(2)-feed_end_range(1))*DoE(:,3));
DoE(:,4) = mean(feed_rate_range)+0.5*(feed_rate_range(2)-feed_rate_range(1))*DoE(:,4);
DoE(:,5) = mean(VCD_start_range)+0.5*(VCD_start_range(2)-VCD_start_range(1))*DoE(:,5);

% define local cluster
if NCpus > 0
    p = gcp('nocreate');
    if isempty(p) || (p.NumWorkers ~= NCpus)
        if ~isempty(p), delete(p); end
        c = parcluster('local');
        c.NumWorkers = NCpus;
        parpool(c);
    end
end

% run experiments
N = size(DoE,1);
f_DoE = zeros(N,1,'double');
BWU = zeros(N,75,'double');
OWU_cell = cell(N,1);
Conc = cell(N,1);
if NCpus > 1
    parfor nexp = 1:N
        [BWU(nexp,:),OWU_cell{nexp},f_DoE(nexp),Conc{nexp}] = ...
            run_single_exp(DoE(nexp,:),mdl_param, ...
            inputs_std_error,outputs_std_error);
    end
else
    for nexp = 1:N
        [BWU(nexp,:),OWU_cell{nexp},f_DoE(nexp),Conc{nexp}] = ...
            run_single_exp(DoE(nexp,:),mdl_param, ...
            inputs_std_error,outputs_std_error);
    end
end

% construct the OWU matrix
OWU = zeros(15*N,5,'double');
for nexp = 1:N
    OWU((nexp-1)*15+1:nexp*15,:) = OWU_cell{nexp};
end

end

function [BWU,OWU_cell,f_DoE,Conc] = run_single_exp(DoE, ...
    mdl_param,inputs_std_error,outputs_std_error)

proc_param_DoE = struct;
proc_param_DoE.GlcConc0 = DoE(1);
proc_param_DoE.GlcFeedStart = DoE(2);
proc_param_DoE.GlcFeedEnd = DoE(3);
proc_param_DoE.GlcFeedRate = DoE(4);
proc_param_DoE.VCD0 = DoE(5);
[conc_profile,BWU,OWU_cell] = run_experiment(mdl_param,proc_param_DoE, ...
    'inputs_std_error',inputs_std_error,'outputs_std_error',outputs_std_error);
f_DoE = conc_profile(end,5);
Conc = conc_profile;

end

