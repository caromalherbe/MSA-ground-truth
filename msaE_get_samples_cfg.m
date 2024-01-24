function result = msaE_get_samples_cfg(num_perms, num_elements, seed)
% Computes a random sample of permutations, and creates the appropriate
% multi-perturbation configurations whose performance should be measured.
% Returns the data in the configuration-wise format.
%
% Syntax: result = msaE_get_samples_cfg(num_perms, num_elements, [seed])
%
% Input:
%   num_perms    - number of requested permutations.
%   num_elements - number of elements in the system.
%   [seed]       - (optional) random seed to use. The default is to use a
%                  random seed based on the current time.
%
% Output: 
%   result       - in the configuration-wise format, includes the
%                  fields: perms.seeds, perms.sizes, perms.matrix;
%                     and: configs
%
%                  result.perms.matrix can be rmfield() to save 
%                  space at the expense of longer running times.
%                  (see reference manual for details).

% The MSA matlab package, written by David Deutscher, June 2004.

% GENERAL STUFF
error(nargchk(2,3,nargin));
msa_internal_global_consts;

% handle inputs:
if nargin < 3
    seed = sum(100*clock);
end

% shorthand
N = num_elements;

% 1. create permutations
% ----------------------
result = msaE_get_samples(num_perms,N,N,seed);

% 2. create configurations
% ------------------------
% the 'all lesioned' configuration, many times
result.configs = zeros(2 + num_perms*(N-1) , N);
% the 'all intact' as the second config
result.configs(2, :) = 1;

% basic_conf will be the (N-1 x N) matrix of the form:
%   0 1 1 1 1
%   0 0 1 1 1
%   0 0 0 1 1
%   0 0 0 0 1
basic_conf = triu( ones(N-1,N), 1 );

for loop = 1:num_perms
    if msa_cnst_display_level >= 1
      	disp(sprintf('computing permutation number %d',loop));
    end
    
    % add the N-1 configurations: without the 'all lesioned' and 'all intact'
    place = 3 + (loop-1)*(N-1);
    result.configs(place:place+N-2, result.perms.matrix(loop,:)) = basic_conf;
end

% 3. eliminate duplicate configurations
% -------------------------------------
if msa_cnst_display_level >= 1
    disp('Removing duplicates...');
end
% sort
result.configs = sortrows(result.configs);
% find rows that are the same as the previous row:
dup = [0 ; all( diff(result.configs) == 0 , 2)];
% keep only unique configurations:
result.configs = result.configs(~dup,:);