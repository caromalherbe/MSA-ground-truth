% FOR INTERNAL USE by the MSA matlab package.

% The MSA matlab package, written by David Deutscher, June 2004.

% define these to be common to all functions
global msa_fmt_cfgF_wise;
global msa_fmt_cfgE_wise;
global msa_fmt_prm_wise;
global msa_fmt_prmSS_only;
global msa_cnst_extra_checks;
global msa_cnst_display_level;

% THIS CAN BE CHANGED TO 0/1 TO SUPPRESS/ALLOW EXTRA ASSERTION IN THE CODE
msa_cnst_extra_checks = 1;

% THIS CAN BE CHANGED TO 0/1/2 FOR NO/PARTIAL/FULL DISPLAY
msa_cnst_display_level = 1;

% THESE ARE INTERNAL CONSTANTS
msa_fmt_cfgF_wise   = 1;
msa_fmt_cfgE_wise   = 2;
msa_fmt_prm_wise    = 3;
msa_fmt_prmSS_only  = 4;
