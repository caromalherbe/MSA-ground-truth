function [num_configs, num_perms, num_tasks] = msa_internal_checks(x, fmt)
% FOR INTERNAL USE by the MSA matlab package.

% The MSA matlab package, written by David Deutscher, June 2004.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OVERVIEW, BY FORMAT
% =-=-=-=-=-=-=-=-=-=
% CFG/FULL
% 1. mandatory num_elements
% 2. mandatory configs + sort rows
% 3. mandatory perfs
% 4. consistency: size(configs) = [2^num_elements num_elements]
%                 size(perfs,1) = 2^num_elements
% 
% CFG/ESTIMATED
% 1. mandatory num_elements
% 2. mandatory configs + sort rows
% 3. mandatory perfs
% 4. mandatory perms + content validity + size(..,2) = num_elements
% 5. consistency: size(configs) = [size(perfs,1) num_elements]
%                 
% PRM
% 1. num_elements:  existance of mandatory field
% 2. max_perturbed: existance of mandatory field 
%     + 1 <= max_perturbed <= num_elements
% 3. perfs:         existance of mandatory field 
% 4. perms:         existance of mandatory field
%     + content validity: (a) either matrix or seeds/sizes 
%                         (b) dimensions of seeds/sizes
% 5. consistency: size(perfs,1) = num_perms
% 6. consistency:   size(perfs,2) = size(perms,2)+1
% 7. IF max_perturbed < num_elements
%     (a) size(perfs,2) >= max_perturbed+1
%     (b) size(perms,2) >= max_perturbed
%
% SS-only
% 1. perms exist and contain sizes + seeds
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


msa_internal_global_consts;

% arg check
error(nargchk(2,2,nargin));


% mandatory field
if ~isfield(x, 'num_elements')
    error('Expecting data to include the mandatory field ''num_elements''');
end

% optional fields
switch fmt

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case msa_fmt_cfgF_wise  % Full CFG_WISE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isfield(x, 'configs') | ~isfield(x, 'perfs')
            error('Expecting data in configuration-wise format');
        end

        if msa_cnst_extra_checks
            % check row order of configs
            [sorted, idx] = sortrows(x.configs);
            if any(any(x.configs ~= sorted))
                warning('Configs matrix must be row-sorted. Sorting...');
                x.configs = sorted;
                x.perfs = x.perfs(idx,:);
            end

            % check they're all binary
            if any(any(x.configs ~= 0 & x.configs ~= 1))
                error('Expecting a binary configs matrix');
            end
        end
        
        % return values
        num_configs = size(x.configs,1);
        num_perms   = nan;
        num_tasks   = size(x.perfs,2);
        
        % consistency
        if num_configs ~= 2^x.num_elements | num_configs ~= size(x.perfs,1)
            error('Inconsistent matrix sizes: ''perfs'' and ''configs'' must have 2^num_elements rows');
        end
        if size(x.configs,2) ~= x.num_elements
            error('Inconsistent matrix sizes: ''configs'' should have ''num_elements'' columns');
        end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case msa_fmt_cfgE_wise  % Estimated CFG-WISE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isfield(x, 'configs') | ~isfield(x, 'perfs')
            error('Expecting data in configuration-wise format');
        end
        
        if msa_cnst_extra_checks
            % check row order of configs
            [sorted, idx] = sortrows(x.configs);
            if any(any(x.configs ~= sorted))
                warning('Configs matrix must be row-sorted. Sorting...');
                x.configs = sorted;
                x.perfs = x.perfs(idx,:);
            end
            
            % check they're all binary
            if any(any(x.configs ~= 0 & x.configs ~= 1))
                error('Expecting a binary configs matrix');
            end            
        end

        % return values
        num_configs = size(x.configs,1);
        [num_perms, perm_cols_check] = subcheckPerms(x,'full');
        num_tasks   = size(x.perfs,2);
        
        % consistency
        if num_configs ~= size(x.perfs,1)
            error('Inconsistent matrix sizes: ''perfs'' and ''configs'' must have equal row counts');
        end
        if size(x.configs,2) ~= x.num_elements
            error('Inconsistent matrix sizes: ''configs'' should have ''num_elements'' columns');
        end
        if ~isnan(perm_cols_check) & perm_cols_check ~= x.num_elements
            error('Inconsistent matrix sizes: ''perms.matrix'' should have ''num_elements'' columns');
        end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case msa_fmt_prm_wise   % PERM-WISE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (2,3,4) existance
        if ~isfield(x, 'max_perturbed') | ~isfield(x, 'perfs') | ~isfield(x, 'perms')
            error('Expecting data in permutation-wise format');
        end
        
        % (2)
        if (x.max_perturbed < 1) | (x.max_perturbed > x.num_elements) | floor(x.max_perturbed) ~= x.max_perturbed
            error('Field ''max_perturbed'' must contain an integer in the range [1,num_elements]');
        end
        
        % return values, + (4) content of 'perms'
        num_configs = nan;
        [num_perms, perm_cols_check] = subcheckPerms(x,'full');
        num_tasks   = 1;

        % (5,6) consistency
        if size(x.perfs,1) ~= num_perms
            error('Inconsistent matrix sizes: ''perms.matrix'' and ''perfs'' must have equal row counts');
        end
        if ~isnan(perm_cols_check) & (perm_cols_check+1 ~= size(x.perfs,2))
            error('Inconsistent matrix sizes: ''perms.matrix'' and ''perfs'' must have compatible column counts');
        end
        
        % (7) if k < n, check that perfs & perms include enough columns
        if x.max_perturbed ~= x.num_elements
            if size(x.perfs,2) < x.max_perturbed+1
                error('Inconsistent matrix sizes: ''perfs'' must have at least ''max_perturbed''+1 columns');
            end
            % for perms, (7) is implied from (6) and the previous check:
            %if ~isnan(perm_cols_check) & perm_cols_check < x.max_perturbed
            %    error('Inconsistent matrix sizes: ''perms.matrix'' must have at least ''max_perturbed'' columns');
            %end
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%
    case msa_fmt_prmSS_only % check perms.sizes & perms.seeds
    %%%%%%%%%%%%%%%%%%%%%%%
        % return values
        num_configs = nan;
        num_perms = subcheckPerms(x,'SSexistance_only');
        num_tasks = 1;
        
    otherwise
        error('Unknown format option');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rows, cols] = subcheckPerms(x, type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(x, 'perms')
    error('Expecting data to have the ''perms'' field');
end

switch type
    case 'full'
        
        alternative = 0;

        % one of the alternatives:
        if ~isfield(x.perms, 'matrix') & ~(isfield(x.perms, 'seeds') & isfield(x.perms, 'sizes'))
            error('Expecting the field ''perms'' to include either subfield ''matrix'' or subfields ''seeds'' and ''sizes''');
        end
        
        % Check the alternative representation:
        if isfield(x.perms, 'seeds') & isfield(x.perms, 'sizes')
            if  length(x.perms.seeds) ~= prod(size(x.perms.seeds)) | ...
                length(x.perms.sizes) ~= prod(size(x.perms.sizes))
                error('Subfields ''seeds'' and ''sizes'' must be vectors');
            end
            if length(x.perms.seeds) ~= length(x.perms.sizes)
                error('Subfields ''seeds'' and ''sizes'' must have equal dimensions');
            end
            rows = sum(x.perms.sizes(:));
            cols = nan;
        end
        
        % Check the prefered representation:
        if isfield(x.perms, 'matrix')
            [rows, cols] = size(x.perms.matrix);
        end
            
    case 'SSexistance_only'
        if ~isfield(x.perms, 'seeds') | ~isfield(x.perms, 'sizes')
            error('Expecting the field ''perms'' to include subfields ''seeds'' and ''sizes''');
        end
        rows = sum(x.perms.sizes);
        cols = nan;
end