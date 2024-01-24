function result = msaE_compute_CVs(msaData)
% Computes the Shapley contribution values when only partial information 
% of the multi-perturbations performance is available. This function
% expects the data in the permutation-wise format.
% This function (and format) supports the computation of k-limited
% perturbations CVs (kpCVs), when concurrently perturbing up to k elements.
%
% Syntax: result = msaE_compute_CVs_cfg(msaData)
%
% Input:
%   msaData     - Data in the permutation-wise format (note the mandatory
%                 field 'max_perturbed', indicating the value of k for kpCVs).
%
% Output:
%   result      - A structure containing the computed CVs. Using the standard
%                 format, this struct contains the fields, 'result.sh' and
%                 'result.se', vectors of the CVs and standard estimation error 
%                 per element.
%                 The struct might include an extra field, 'result.counts', 
%                 which is a vector indicating the number of marginal 
%                 contributions averaged per element. 
%                 This field appears only when computing kpCVs with k~=n
%                 (since if k==n, each element appears exactly once in
%                 each permutation, and thus 'counts' will contain the 
%                 number of sampled permutations for all elements).

% The MSA matlab package, written by David Deutscher, June 2004.

% GENERAL STUFF
error(nargchk(1,2,nargin));
msa_internal_global_consts;
[dummy1, num_perms, dummy2] = msa_internal_checks(msaData, msa_fmt_prm_wise);

% shorthand
N = msaData.num_elements;

% get the permutation matrix
if ~isfield(msaData.perms,'matrix')
    msaData = msaE_recreate_samples(msaData);
end

%%%%%%%%%%%
% PAHSE 1 %
%%%%%%%%%%%

% display
if (msa_cnst_display_level >= 1), disp('Initializing'), end

% use only the first k perturbations, if there are more
if size(msaData.perms.matrix,2) > msaData.max_perturbed
    msaData.perms.matrix = msaData.perms.matrix(:, 1 : msaData.max_perturbed );
end
if size(msaData.perfs,2) > msaData.max_perturbed+1
    msaData.perfs        = msaData.perfs(:, 1 : msaData.max_perturbed+1 );
end

% compute unordered deltas
deltas = -diff(msaData.perfs,1,2);   % 1st order along dimension 2 (column);

% DEBUG ASSERTION:
if any(size(deltas) ~= size(msaData.perms.matrix))
    error('Assertion failed: deltas & matrix have unequal dimensions');
end

% change to column vectors
deltas = deltas(:);
msaData.perms.matrix = msaData.perms.matrix(:);

%%%%%%%%%%%
% PAHSE 2 %
%%%%%%%%%%%

% display
if (msa_cnst_display_level >= 1), disp('Counting appearances'), end

% sort: make all deltas of each element consecutive in the ordering
[element,order]=sort(msaData.perms.matrix);
deltas = deltas(order);

% trick #1: quickly count how many deltas each element have
notAsThePrevious = diff(double(element));
lastPositions = [find(notAsThePrevious) ; length(element)];
counts = diff( [0 ; lastPositions] )';

% display
if (msa_cnst_display_level >= 1), disp('Reordering deltas'), end

% trick #2, BAD version: sum the appropriate deltas for each element.
% problem: we loose accuracy by using this large cumsum !
%cumulative_sum = cumsum(deltas);
%sums = diff( [0 ; cumulative_sum(lastPositions)] );

% trick #2 (revised): order the deltas into columns, one column for each
% element, by creating the proper indice vectors:
col_indices = double(element);
row_indices = ones(size(deltas));
row_indices(find(notAsThePrevious)+1) = (-counts(1:end-1)+1);
row_indices = cumsum(row_indices);

ordered_deltas = zeros(max(counts), N);
ordered_deltas( sub2ind(size(ordered_deltas),row_indices,col_indices) ) = deltas;

%%%%%%%%%%%
% PAHSE 3 %
%%%%%%%%%%%

% display
if (msa_cnst_display_level >= 1), disp('Computing'), end

% handle the case where not all elements has deltas
if length(counts) ~= N
    warning('Some of the elements have no sampled marginal contributions. Their CVs cannot be estimated');
    % find which elements do have something
    participating_elemets = element(lastPositions);
    
    % insert zeros into proper places
    temp = zeros(1, N);
    temp(participating_elemets) = counts;
    counts = temp;
    clear temp;
end

if any(counts == 1)
    warning('Some of the elements have only a single marginal contribution. SE was set to NaN');
end

% If k == n, the matrices perfs & perms might have been truncated
% to save the zeros in the end - but all counts, for all elements 
% are KNOWN to be exactly num_perms
if msaData.max_perturbed == msaData.num_elements
    [result.sh, result.se] = mean_and_se(ordered_deltas,num_perms);
else
    [result.sh, result.se] = mean_and_se(ordered_deltas,counts);
    result.counts = counts;
end

% display
if (msa_cnst_display_level >= 1), disp('Finished'), end

% END OF FUNCTION %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m,s] = mean_and_se(data,counts)
% Compute mean/se for a given sample, with given counts per element.
% Assume any missing or extra data is zero.
% If count==0  =>  m = s = NaN;
% If count==1  =>  s = NaN;

% stats
[rows,elements] = size(data);
if prod(size(counts)) == 1
    % scalr expansion
    counts = repmat(counts,1,elements);
end
zeroed = (counts == 0);
two_and_more = (counts >= 2);

% MEAN:
m = repmat(NaN, 1, elements);
m(~zeroed) = sum(data(:,~zeroed),1) ./ counts(~zeroed);

% STD: (better accuracy than vectorized versions)
s = repmat(NaN, 1, elements);
for e=find(two_and_more)
    if counts(e) > rows
        s(e) = std([data(:,e); zeros(counts(e)-rows,1)]);
    else
        s(e) = std(data(1:counts(e),e));
    end
end
s(two_and_more) = s(two_and_more) ./ sqrt(counts(two_and_more));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS IS THE CLEAN, BUT SLOW, CODE. IT'S HERE AS REFERENCE:               %
%                                                                          %
% % compute shapley + std for each element seperately                      %
% % (since they have different counts)                                     %
% for elem = 1:N                                                           %
%     hits = (perms == elem);                                              %
%     count = nnz(hits);                                                   %
%                                                                          %
%     if nargout >= 3                                                      %
%         all_hits(elem) = count;                                          %
%     end                                                                  %
%                                                                          %
%     if  count == 0                                                       %
%         warning(sprintf('Element #%d did not appear at all',elem));      %
%         sh(elem) = nan;                                                  %
%     else                                                                 %
%         sh(elem) = mean(deltas(hits));                                   %
%         % second arg of std(): divide by sqrt(count) and not sqrt(N-1)   %
%         sd(elem) = std(deltas(hits),1) / sqrt(count);                    %
%                                                                          %
%     end                                                                  %
% end                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%