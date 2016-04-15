STACK = dbstack('-completenames');
[THISDIR, ~, ~] = fileparts(STACK(1).file);

init_cyt;

addpath(fullfile(THISDIR, 'pipeline'));
