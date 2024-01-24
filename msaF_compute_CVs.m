function result = msaF_compute_CVs(msaData)
% Computes the Shapley contribution values when full information 
% of the multi-perturbations performance is available.
%
% Syntax: result = msaF_compute_CVs(msaData)
%
% Input:
%   msaData     - Data in the configuration-wise format, with the full set
%                 of 2^n configurations available.
%
% Output:
%   result      - A structure containing the computed CVs. Using the standard
%                 format, this struct contains a single field, 'result.sh',
%                 which is a vector of CVs per element. In the full
%                 information situation there is no estimation error, hence
%                 no '*.sd' field.

% The MSA matlab package, written by David Deutscher, June 2004.

% GENERAL STUFF
error(nargchk(1,1,nargin));
msa_internal_global_consts;
[num_configs, dummy, num_tasks] = msa_internal_checks(msaData, msa_fmt_cfgF_wise);

% shorthand
N = msaData.num_elements;

% number of intact elements per configuration
intact_count = sum(msaData.configs,2);

% accessories
all_factorials = [0 1 cumprod(1:N)]; % all_factorials(i) = (i-2)! [where (-1)! is a fictitious unused value]
count_plus  = all_factorials(intact_count-1 +2) .* all_factorials(N-intact_count   +2);
count_minus = all_factorials(intact_count   +2) .* all_factorials(N-intact_count-1 +2);
totals      = count_plus * msaData.configs;

% fast computation of CVs
result.sh = ((repmat(count_plus,num_tasks,1)  .* msaData.perfs') * msaData.configs - ... 
      (repmat(count_minus,num_tasks,1) .* msaData.perfs') * (~msaData.configs)) ...
     ./ repmat(totals,num_tasks,1);
