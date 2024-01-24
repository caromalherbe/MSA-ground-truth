function [found, place] = msa_internal_find_row(config, matrix, max)
% FOR INTERNAL USE by the MSA matlab package.

% The MSA matlab package, written by David Deutscher, June 2004.

% SYNTAX: [found, place] = find_row(config, matrix, [last_row])
%         Tries to find the row vector 'config' in 'matrix', using binary
%         search.
%         'last_row' can be specified to limit the search to rows
%         1 to 'last_row' of 'matrix'.
% ASSUMES: 1. 'matrix' is sortrow()-ed, 2D matrix.
%          2. 'config' is a row vector.
%          3. both have the same number of columns.
% returns: found   - 1 if the row is found, 0 otherwise.
%          place   - The row number in 'matrix'.

if nargin < 3
    max = size(matrix,1);
end
min=1;
while max-min >= 0
    mid = floor((min+max)/2);
    c = compare(config,matrix(mid,:));
    switch(c)
    case 0
        found = 1;
        place = mid;
        return;
    case -1
        max = mid-1;
    case 1
        min = mid+1;
    end
end
found = 0;
if c < 0
    place = mid;
else
    place = mid+1;
end


function c=compare(a,b)
i = find(a~=b);
if isempty(i)
    c = 0;
elseif a(i(1)) > b(i(1))
    c = 1;
else
    c = -1;
end