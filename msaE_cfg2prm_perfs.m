function cperf = msaE_cfg2prm_perfs(msaData,skipchecks);
% Reorders the performance scores into the matrix format needed for the
% permutation-wise format. The output is also useful for plotting the
% performance along consecutive lesions, from all-intact to all-lesioned.
%
% Syntax: cperf = msaE_cfg2prm_perfs(msaData)
%
% Input:
%   msaData      - Data in the configuration-wise, estimated format (i.e.,
%                  including the field msaData.perms).
%
% Output: 
%   cperf        - The performance matrix in the permutation-wise format:
%                  one row per permutation, including the performance
%                  of the all-intact system and then scores for consecutive 
%                  lesions of elements.

% The MSA matlab package, written by David Deutscher, June 2004.

% GENERAL STUFF
error(nargchk(1,2,nargin));

% DO INPUT CHECKS
if nargin < 2 | ~skipchecks
    msa_internal_global_consts;
    [num_configs, num_perms, num_tasks] = msa_internal_checks(msaData, msa_fmt_cfgE_wise);

    % get the permutation matrix
    if ~isfield(msaData.perms,'matrix')
        msaData = msaE_recreate_samples(msaData);
    end
    
% OR DON'T DO INPUT CHECKS
else
    num_configs = size(msaData.configs, 1);
    num_perms   = size(msaData.perms.matrix, 1);
    num_tasks   = size(perfs, 2);
end

% shorthand
num_elements = msaData.num_elements;

% display
if (msa_cnst_display_level >= 1), disp('Ordering performance scores'), end

% get the all intact scores (per task)
[u,all_intact_index] = msa_internal_find_row(ones(1,num_elements), msaData.configs);
if ~u
    error('configuration matrix must include the ''all intact'' config.');
end
all_intact_perfs = msaData.perfs(all_intact_index,:);

% put tasks in the 3rd dimension
all_intact_perfs = shiftdim(all_intact_perfs,1); 

% allocate memory
cperf = zeros(num_perms, num_elements+1, num_tasks);

% duplicate the 'all_intact' scores to the first column of all rows:
cperf(:,1,:) = repmat(all_intact_perfs,[num_perms,1,1]);

% for every permutation, arrange all other scores
for permutation_no = 1:num_perms
    
    if msa_cnst_display_level >= 2
        disp(sprintf('Ordering permutation #%d',permutation_no));
    end
    
    config = ones(1,num_elements);
    % walk the permutation element-by-element
    for j = 1:num_elements
        
        % the next element
        element = msaData.perms.matrix(permutation_no,j);
        
        % lesion this element
        config(element) = 0;
        
        % get the corresponding scores
        [u,r] = msa_internal_find_row(config,msaData.configs);
        if ~u
            error(['Missing configuration: ' num2str(config)]);
        end

        % and put them in place
        cperf(permutation_no, j+1, :) = shiftdim(msaData.perfs(r,:),1);
    end
end