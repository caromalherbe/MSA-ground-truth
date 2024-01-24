function result = msaE_get_samples(num_perms, n, k, seed)
% Computes a random sample of permutations.
% Returns the data in the permutation-wise format.
%
% Syntax: result = msaE_get_samples(num_perms, n, [k], [seed])
%
% Input:
%   num_perms       - number of requested permutations.
%   num_elements    - number of elements in the system.
%   [max_perturbed] - (optional) return only permutation prefixes of
%                     size 'max_perturbed', used for kpCVs; 
%                     The default is max_perturbed==num_elements.
%   [seed]          - (optional) random seed to use. The default is to use a
%                     random seed based on the current time.
%
% Output: 
%   result          - in the permutation-wise format, includes the
%                     fields: perms.seeds, perms.sizes, perms.matrix;
%
% result.perms.matrix can be rmfield() to save 
% space at the expense of longer running times.
% (see package manual for details).

% The MSA matlab package, written by David Deutscher, June 2004.

error(nargchk(2,4,nargin));

% handle inputs:
if nargin < 3
    k = n;
elseif k > n | k < 1
    error('Invalid k: must be in [1,num_elements]');
end

if nargin < 4
    result.perms.seeds = sum(100*clock);
else
    result.perms.seeds = seed;
end

% 
result.perms.sizes = num_perms;
result.num_elements = n;
result.max_perturbed = k;

% init random number generator
rand('state',result.perms.seeds);

% compute random permutations, keeping only the first k elements in each
result.perms.matrix = zeros(num_perms, k);
for pidx = 1:num_perms
    perm = randperm(n);
    result.perms.matrix(pidx,:) = perm(1:k);
end
