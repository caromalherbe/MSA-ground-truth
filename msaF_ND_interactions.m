function [interactions, dividends] = msaF_ND_interactions(msaData)
% Computes the interactions of all subsets of elements (i.e., interactions
% in all dimensions/orders, 1..num_elements), given full information 
% of all multi-perturbations performance.
%
% Syntax: [interactions, dividends] = msaF_ND_interactions(msaData)
%
% Input:
%   msaData    - Data in the configuration-wise format, with the full set
%                of 2^n configurations available.
%
% Output:
%   interactions - 
%       A struct array of size (1 x num_elements). 
%       interactions(d) describes the interactions of order d (interactions 
%       between element groups of size d).
%       Each cell contains two fields:
%        1. interactions(d).elements   - is the list of elements in each group.
%           Each group is a row with d columns, for a total of (n choose d)
%           rows.
%        2. interactions(d).inter_size - is a column vector of the interaction 
%           sizes for the corresponding groups.
%           Note that the singleton interactions (d==1) are exactly the CVs 
%           returned by msaF_compute_CVs().
%
%   dividends -
%       A column vector of length 2^n containing the dividends, ordered
%       according the configruations matrix.

% The MSA matlab package, written by David Deutscher, June 2004.
% The algorithm for this function is based on ideas from:
% Grabisch M., Marichal J.L., Roubens M.
% "Equivalent Representations of Set Functions"

% GENERAL STUFF
error(nargchk(1,1,nargin));
msa_internal_global_consts;
[num_configs, num_perms, num_tasks] = msa_internal_checks(msaData, msa_fmt_cfgF_wise);

if num_tasks > 1
    error('Currently not supporting multiple tasks');
end

% shorthand
N = msaData.num_elements;

% allocate memory
dividends = zeros(num_configs,1);

% calculate the dividends
counts = sum( msaData.configs, 2 );
for i = 1:num_configs
    S = msaData.configs(i,:);
    % find configs contained in S
    Tidx = ~any(msaData.configs(:,S==0),2);
    % and compute
    dividends(i) = sum( ((-1) .^ (sum(S)-counts(Tidx))) .* msaData.perfs(Tidx) );
end

% for each dimension
for dim = 1:N
    % find all groups of size dim
    groups = find(counts == dim);
    
    % allocate memory
    interactions(dim).elements = zeros(length(groups), dim);
    interactions(dim).inter_size = zeros(length(groups), 1);
    % for each group
	for gidx = 1:length(groups)
       	group = msaData.configs(groups(gidx),:);
        % find configs containing the group
    	Tidx = all(msaData.configs(:,group~=0), 2);
        % and compute
        interactions(dim).elements(gidx,:) = find(group);  
        interactions(dim).inter_size(gidx) = sum( dividends(Tidx) ./ (counts(Tidx) - dim + 1) );
    end
end