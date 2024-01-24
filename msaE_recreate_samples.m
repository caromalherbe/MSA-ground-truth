function msaData = msaE_recreate_samples(msaData, compact_storage)
% Recreates a random sample of permutations, sampled at some earlier time,
% given the seeds for random-number-generator.
%
% Syntax: result = msaE_recreate_samples(msaData, [compact_storage=1])
%
% Input:
%   msaData      - Data in the permutation-wise format, including
%                  the fields 'seeds' and 'sizes' that allow the 
%                  recomputation of the random sample of permutations.
%   [compact_storage] - optional. If compact_storage==1, to save storage, 
%                  only the needed prefix of each permutation is saved, 
%                  of size: min(num_elements ; max_perturbed ; size(perfs,2)-1)
%                  if the appropriate fields exist in the input.
%                  Set compact_storage=0 to compute the full permutations. 
%
% Output: 
%   result       - A copy of the input structure, but with the field
%                  result.perms.matrix filled with the recomputed
%                  permutations.
                
% The MSA matlab package, written by David Deutscher, June 2004.

% GENERAL STUFF
error(nargchk(1,2,nargin));
msa_internal_global_consts;

if nargin < 2
    compact_storage = 1;
end

% compute only if seeds/sizes exist
msa_internal_checks(msaData, msa_fmt_prmSS_only);

% display
if (msa_cnst_display_level >= 1), disp('Re-creating the permutations'), end

% prefix of permutation to save = minimum(n, k, size(perfs,2)-1)
prefix_size = msaData.num_elements;
if compact_storage
    if isfield(msaData,'max_perturbed')
        prefix_size = msaData.max_perturbed;
    end
    if isfield(msaData,'perfs')
        prefix_size = min(prefix_size, size(msaData.perfs,2)-1);
    end
end

% allocate memory
msaData.perms.matrix = zeros(sum(msaData.perms.sizes),prefix_size);

% recreate
next = 1;
for i = 1:length(msaData.perms.seeds);
    batch = msaE_get_samples(msaData.perms.sizes(i), ...
                             msaData.num_elements, ...
                             prefix_size, ...
                             msaData.perms.seeds(i));
    msaData.perms.matrix(next : next+msaData.perms.sizes(i)-1, 1:prefix_size) = batch.perms.matrix;
    next = next+msaData.perms.sizes(i);
end