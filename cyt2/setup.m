STACK = dbstack('-completenames');
[THISDIR, ~, ~] = fileparts(STACK(1).file);
clear('STACK');

path_with_subdirectories = genpath(THISDIR);
addpath(path_with_subdirectories);

% although the pipeline directory should by now be in the MATLAB path, the next
% line ensures that it comes ahead of all the others
addpath(fullfile(THISDIR, 'pipeline'));

clear('THISDIR');
% -----------------------------------------------------------------------------

global DEBUG_REPRODUCIBILITY;

global DEFAULT_PRNG_SEED
DEFAULT_PRNG_SEED = 1;
