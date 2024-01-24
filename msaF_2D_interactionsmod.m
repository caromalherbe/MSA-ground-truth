function [inter_size, CV_col_intact, CV_col_lesioned, CV_both] = msaF_2D_interactionsmod(msaData)
% Computes the 2-dimensional interactions between all pairs of elements,
% given full information of all multi-perturbations performance.
%
% Syntax: [inter_size, CV_col_intact, CV_col_lesioned] = msaF_2D_interactions(msaData)
%
% Input:
%   msaData    - Data in the configuration-wise format, with the full set
%                of 2^n configurations available.
%
% Output: 
%   inter_size - 
%       The symmetric matrix of interactions between all pairs of elements.
%       Notice: inter_size = (CV_col_intact - CV_col_lesioned).
%
%   CV_col_intact -
%       A matrix: CV_col_intact(i,j) is the CV of element i in
%       the subgame where j is always intact.
%
%   CV_col_lesioned -
%       A matrix: CV_col_lesioned(i,j) is the CV of element i in
%       the subgame where j is always lesioned.

% The MSA matlab package, written by David Deutscher, June 2004.

% GENERAL STUFF
error(nargchk(1,1,nargin));
msa_internal_global_consts;
[num_configs, num_perms, num_tasks] = msa_internal_checks(msaData, msa_fmt_cfgF_wise);

% shorthand
N = msaData.num_elements;

% allocate memory
CV_both     = zeros(N);
CV_onlyi    = zeros(N);
CV_onlyj    = zeros(N);

% init sub games 
sub_game_both.num_elements = msaData.num_elements - 1;
sub_game_onlyi.num_elements = msaData.num_elements - 1;
sub_game_onlyj.num_elements = msaData.num_elements - 1;

for i=1:N
    for j=i+1:N
        % measuring the interaction between element i and element j.
        if (msa_cnst_display_level >= 2)
            disp(sprintf('interaction of %d <-> %d',i,j));
        end
        
        i_lesioned  = ~msaData.configs(:,i);
        j_lesioned  = ~msaData.configs(:,j);
        both_equal  = (msaData.configs(:,i) == msaData.configs(:,j));

        % sanity check:
        if sum(both_equal)~=num_configs/2 | sum(j_lesioned)~=num_configs/2 | sum(i_lesioned)~=num_configs/2
            error ('Incorrect configurations set');
        end
        
        % create subgames
        sub_game_onlyi.configs = msaData.configs( j_lesioned , [1:j-1,j+1:end]);
        sub_game_onlyi.perfs = msaData.perfs( j_lesioned , :);
        
        sub_game_onlyj.configs = msaData.configs( i_lesioned , [1:i-1,i+1:end]);
        sub_game_onlyj.perfs = msaData.perfs( i_lesioned , :);

        % when both have equal status, location j is dropped & location i stands 
        % for the new compund player
        sub_game_both.configs  = msaData.configs( both_equal , [1:j-1,j+1:end]);
        sub_game_both.perfs  = msaData.perfs( both_equal , :);

        % compute player's contribution in the subgames
        res = msaF_compute_CVs(sub_game_onlyi);
        CV_onlyi(i,j) = res.sh(i);
        
        res = msaF_compute_CVs(sub_game_onlyj);
        CV_onlyj(i,j) = res.sh(j-1); % always j>i, so after removing i, the j-th player is in the j-1 position
        
        res = msaF_compute_CVs(sub_game_both);
        CV_both(i,j) = res.sh(i);
    end
end

% compute interactions:
% inter_size(i,j) is the interaction between elements i & j, defined as:
%   SH(i, with j intact) - SH(i, with j lesioned)
% OR, equivalently
%   SH(i&j together) - SH(i, with j lesioned) - SH(j, with i intact)
% error I think, SH(j, with i lesioned)
CV_col_lesioned = CV_onlyi + CV_onlyj';
CV_col_intact = (CV_both - CV_onlyj) + (CV_both - CV_onlyi)';
inter_size = CV_col_intact - CV_col_lesioned;
