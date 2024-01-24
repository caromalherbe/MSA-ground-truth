function result = msaE_compute_CVs_cfg(msaData)
% Computes the Shapley contribution values when only partial information 
% of the multi-perturbations performance is available. This function
% expects the data in the configuration-wise format.
%
% Syntax: result = msaE_compute_CVs_cfg(msaData)
%
% Input:
%   msaData      - Data in the configuration-wise, estimated format (i.e.,
%                  including the field msaData.perms).
%
% Output:
%   result      - A structure containing the computed CVs. Using the standard
%                 format, this struct contains the fields, 'result.sh' and
%                 'result.se', vectors of the CVs and standard estimation error 
%                 per element.

% The MSA matlab package, written by David Deutscher, June 2004.

% GENERAL STUFF
error(nargchk(1,1,nargin));
msa_internal_global_consts;
[num_configs, num_perms, num_tasks] = msa_internal_checks(msaData, msa_fmt_cfgE_wise);

% get the permutation matrix
if ~isfield(msaData.perms,'matrix')
    msaData = msaE_recreate_samples(msaData);
end

% display
if (msa_cnst_display_level >= 1), disp('Computing'), end

% get the performance matrix, ordered by permutations 
consecutive_perf = msaE_cfg2prm_perfs(msaData, 0); % 0 for skipping the input checks

% compute deltas:
consecutive_deltas = -diff(consecutive_perf,1,2); % 1st order diff, dimension 2 (column-wisw)

% allocate memory
deltas = zeros(num_perms, msaData.num_elements, num_tasks);

% permute deltas into place:
[rows, dummy, pages] = ndgrid(1:num_perms,1:msaData.num_elements,1:num_tasks);
deltas( sub2ind(size(deltas),rows,msaData.perms.matrix,pages) ) = consecutive_deltas;
    
% and shapley
result.sh = squeeze(mean(deltas,1));
result.se = squeeze(std(deltas,0,1)) / sqrt(num_perms);