function [] = terse_warning(varargin)
    backtrace = warning('off', 'backtrace');
    verbose = warning('off', 'verbose');
    warning(varargin{:});
    warning('backtrace', backtrace.state);
    warning('verbose', verbose.state);
end

