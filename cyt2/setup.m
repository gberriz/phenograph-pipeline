STACK = dbstack('-completenames');
[THISDIR, ~, ~] = fileparts(STACK(1).file);

init_cyt;

addpath(fullfile(THISDIR, 'pipeline'));

% -----------------------------------------------------------------------------

global DEFAULT_PRNG_SEED
DEFAULT_PRNG_SEED = 1;
