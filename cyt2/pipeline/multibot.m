function [] = multibot(number_of_runs, inputdir, outputdir, samplesize)
% MULTIBOT execute multiple runs of bot script.
%
% MULTIBOT(NUMBER_OF_RUNS, INPUTDIR, OUTPUTDIR, SAMPLESIZE) executes
% NUMBER_OF_RUNS runs of
%
% BOT(INPUTDIR, fullfile(OUTPUTDIR, SUBDIR), SAMPLESIZE)
%
% where SUBDIR is of the form 'run<NUMERIC_SUFFIX>'.

    setup;
    global DEBUG_REPRODUCIBILITY;
    DEBUG_REPRODUCIBILITY = false;

    global DEV_MODE;
    DEV_MODE = false;

    format_ = indexed_format('run', number_of_runs);
    for run_number = 1:number_of_runs
        dirname = sprintf(format_, run_number);
        outputpath = fullfile(outputdir, dirname);
        maybe_create_dir(outputpath);
        bot(inputdir, outputpath, samplesize, false);
    end

end
