function DSGE_setup
% Use the Dynare installation explicitly selected by the user.
assert(~isempty(which('dynare')), ...
    'Add the matlab directory of Dynare 7.1 to the MATLAB path first.');
dynare_config;
assert(strcmp(dynare_version,'7.1'), ...
    'This package uses the Dynare 7.1 interface. Please select Dynare 7.1.');
end
