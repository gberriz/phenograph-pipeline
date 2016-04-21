function [] = bot(inputdir, outputdir, samplesize, savesession)
%BOT run PhenoGraph pipeline on a collection of files
%
%    BOT(INPUTDIR, OUTPUTDIR) locates all the *.fcs files under INPUTDIR, uses
%    them as input for a PhenoGraph pipeline, and stores the results under
%    OUTPUTDIR.
%
%    Two sub-directories will be created under OUTPUTDIR, namely 'full' and
%    'truncated'.  The structure of each of these directories will mirror
%    that of INPUTDIR, except that all the files will have extension .tsv
%    instead of .fcs.  (Files INPUTDIR that don't match *.fcs are ignored.)
%
%    BOT(INPUTDIR, OUTPUTDIR, SAMPLESIZE) behaves similarly, but runs PhenoGraph
%    on a size-SAMPLESIZE random sample of the input data.
%
%    BOT(INPUTDIR, OUTPUTDIR, SAMPLESIZE, true), in addition, results in having
%    the current MATLAB session saved in the file matlab_session.mat, under
%    OUTPUTDIR.
%
%    After calling any one of the forms of BOT shown above, it may be called
%    without arguments; in this case, the last explicit set of arguments passed
%    to it are used.

    global BOT_LAST_ARGS

    if nargin == 0
        if wehave(BOT_LAST_ARGS)
            args = BOT_LAST_ARGS;
        else
            args = {};
        end
        bot(args{:});
        return;
    end

    if ~exist('savesession', 'var')
        if ~exist('samplesize', 'var')
            bot(inputdir, outputdir, [], false);
        else
            bot(inputdir, outputdir, samplesize, false);
        end
        return;
    end

    BOT_LAST_ARGS = {inputdir, outputdir, samplesize, savesession};

    tic
    run_pipeline(inputdir, outputdir, samplesize, savesession);
    toc
end

% -----------------------------------------------------------------------------
% input

function files = get_files(inputdir)

    % a bit of a kluge...
    command = sprintf('cd "%s" && find ./ -name "*.fcs" | sort', inputdir);

    output = run_command(command);

    files = cellfun(@(s) regexprep(s, '^\./', ''), split_lines(output), ...
                    'UniformOutput', false);

end

function [data, block_sizes, channel_names] = read_files(files)
    function [fcs_data, fcs_metadata] = read_fcs(path_)
        fprintf('reading %s ...', path_);
        [fcs_data, fcs_metadata] = fca_readfcs(path_);
        fprintf('done\n');
    end
    [datacell, metadatacell] = cellfun(@read_fcs, files, ...
                                       'UniformOutput', false);
    data = cat(1, datacell{:});
    block_sizes = cellfun(@(d) size(d, 1), datacell);
    clear('datacell');

    parameters = cellfun(@(s) s.par, metadatacell, 'UniformOutput', false);

    function names = get_channel_names(parameter)
        raw_names = cellfun(@char, [{parameter.name}; {parameter.name2}], ...
                            'UniformOutput', false);

        function name = get_channel_name(pair)
            names = unique(pair(~strcmp(pair, '')));
            name = strjoin(names, '__');
        end

        names = cellfun(@get_channel_name, num2cell(raw_names, 1), ...
                        'UniformOutput', false);
    end

    channel_names = cellfun(@get_channel_names, parameters, ...
                            'UniformOutput', false);
end

function culled_relative_paths = cull_sources(relative_paths, ...
                                              maximum_number_of_paths)
    number_of_paths = numel(relative_paths);
    if number_of_paths <= maximum_number_of_paths
        culled_relative_paths = relative_paths;
        if number_of_paths < maximum_number_of_paths
            terse_warning('only %d sources available', number_of_paths);
        end
        return;
    end

    function parts = parse(i)
        relative_path = relative_paths{i};
        [~, stem, ~] = fileparts(relative_path);
        function part = check_match(match)
            if wehave(match)
                part = match{1};
            else
                part = stem;
            end
        end
        prefix = check_match(regexpi(stem, '^([a-z]+)', 'match'));
        suffix = check_match(regexpi(stem, '([a-z]+)$', 'match'));
        parts = {prefix, suffix, i};
    end

    rows = arrayfun(@parse, 1:number_of_paths, 'UniformOutput', false);
    table_ = cell2table(cat(1, rows{:}), ...
                        'VariableNames', {'prefix', 'suffix', 'index'});

    groups = findgroups(table_.suffix, table_.prefix);
    width = max(splitapply(@numel, table_.index, groups));

    function padded_indices = pad(indices)
        padded_indices = [to_row(indices) nan(width - numel(indices), 1)];
    end

    indices = to_column(splitapply(@pad, table_.index, groups));
    indices(isnan(indices)) = [];

    culled_relative_paths = relative_paths(indices(1:maximum_number_of_paths));
end

% -----------------------------------------------------------------------------
% output

function [] = save_to_tsv(path_, data_as_table)
    [dirname, ~, ~] = fileparts(path_);
    maybe_create_dir(dirname);
    writetable(data_as_table, path_, 'Delimiter', '\t', 'FileType', 'text');
end

% -----------------------------------------------------------------------------
% sampling

function sample = random_sample(k, n)
    global DEBUG_REPRODUCIBILITY;
    if DEBUG_REPRODUCIBILITY
        global PRNG_SEED;
        global RANDSAMPLE_PRNG_SEED_OFFSET;
        rng(PRNG_SEED + RANDSAMPLE_PRNG_SEED_OFFSET);

        sample = sort(to_column(randsample(1:n, k)));
    else
        sample = partial_fisher_yates(k, n);
    end
end

function sample = partial_fisher_yates(k, n)
    if k > n
        error('k exceeds n');
    end

    map = containers.Map('KeyType', 'uint64', 'ValueType', 'uint64');
    function value = getvalue(key)
        if isKey(map, key)
            value = map(key);
        else
            value = key;
        end
    end

    sample = zeros(k, 1, 'uint64');
    i = 0;
    last = n;
    while i < k
        i = i + 1;
        pick = randi(last);
        sample(i) = getvalue(pick);
        map(pick) = getvalue(last);
        last = last - 1;
    end
end

function mapper = make_mapper(block_sizes, block_values)
    interpolation_table = [0; cumsum(block_sizes)];
    max_index = interpolation_table(end);
    core = @(i) interp1(interpolation_table, ...
                          [NaN block_values], i, 'next');

    function value = mapper_(index)
        if (index < 1) || (max_index < index)
            error('index out of bounds: %d', index);
        end
        value = core(double(index));
    end

    mapper = @mapper_;
end

% -----------------------------------------------------------------------------
% channels, sources, tsne and phenograph tables

function channels_table = make_channels_table(data)
    width = size(data, 2);
    format = sprintf('ch%%0%dd', floor(log10(width)) + 1);
    channel_columns = arrayfun(@(i) sprintf(format, i), 1:width, ...
                               'UniformOutput', false);
    channels_table = array_to_table(data, channel_columns);
end

function sources_table = make_sources_table(block_sizes, data_indices)

    number_of_sources = numel(block_sizes);
    block_indices = 1:number_of_sources;
    % letters = 'A':'Z';
    % block_values = letters(block_indices);

    if wehave(data_indices)
        observation_index_to_source_index = make_mapper(block_sizes, ...
                                                        block_indices);

        source_indices = arrayfun(@(i) observation_index_to_source_index(i), ...
                                  to_column(data_indices));
    else
        blocks = arrayfun(@(i) i * ones(block_sizes(i), 1), ...
                          to_column(block_indices));
        sources_indices = cat(1, blocks{:});
    end
    sources_table = table(categorical(source_indices), ...
                          'VariableNames', {'source'});
end

function tsne_table = run_tsne(data)

    number_of_observations = size(data, 1);
    MAX_TSNE = 1000000;

    if (number_of_observations > MAX_TSNE)
        error('tsne: too many observations (%d)', ...
              number_of_observations);
    end

    % -------------------------------------------------------------------------

    initial_dims = 110; % taken from the cyt code; I don't know why they
                        % chose this value
    perplexity = min(30, (size(data, 1) - 1)/3);

    tsne_mapping = fast_tsne(data, [], initial_dims, perplexity);

    % -------------------------------------------------------------------------
    variable_names = arrayfun(@(i) sprintf('bh_SNE%d', i), ...
                              1:size(tsne_mapping, 2), ...
                              'UniformOutput', false);

    tsne_table = array_to_table(tsne_mapping, variable_names);
end

function phenograph_table = run_phenograph(data)
    [cluster_labels, ~, ~, ~] = ...
        phenograph(data, 20, 'distance', 'mahalanobis');
    phenograph_table = table(categorical(cluster_labels), ...
                             'VariableNames', {'cluster'});
end

% -----------------------------------------------------------------------------
% table utils

function table_ = array_to_table(array, headers)
    table_ = array2table(array, 'VariableNames', headers);
end

% -----------------------------------------------------------------------------
% string manipulation

function new_string = strip_newline(string_)
    new_string = regexprep(string_, '\n$', '');
end

function lines = split_lines(string_)
    lines = to_column(regexp(strip_newline(string_), '\n', 'split'));
end

% -----------------------------------------------------------------------------
% run shell commands

function output = run_command(command)
    [status, output] = system(command);
    if status ~= 0
        error('"%s" failed with status %d:\n%s', command, status, output)
    end
end

% -----------------------------------------------------------------------------
% semantic sugar

function column = to_column(array)
    column = reshape(array, [], 1);
end

function row = to_row(array)
    row = reshape(array, 1, []);
end

% -----------------------------------------------------------------------------
% main routine

function [] = run_pipeline(inputdir, outputdir, samplesize, savesession)

    global DEV_MODE;
    global DEBUG_REPRODUCIBILITY;
    global PRNG_SEED;
    if DEBUG_REPRODUCIBILITY
        if isempty(PRNG_SEED)
            global DEFAULT_PRNG_SEED;
            if wehave(DEFAULT_PRNG_SEED)
                PRNG_SEED = DEFAULT_PRNG_SEED;
            else
                error('no PRNG seed available');
            end
        end
        terse_warning('PRNG_SEED: %d', PRNG_SEED);
    else
        clear('global PRNG_SEED');
    end

    % -------------------------------------------------------------------------

    if DEV_MODE
        global RELATIVE_PATHS;
        if wehave(RELATIVE_PATHS)
            relative_paths = RELATIVE_PATHS;
        else
            relative_paths = cull_sources(get_files(inputdir), 4);
        end
    else
        relative_paths = get_files(inputdir);
    end

    files = cellfun(@(r) fullfile(inputdir, r), relative_paths, ...
                    'UniformOutput', false);

    if DEBUG_REPRODUCIBILITY
        global FCS_FILES;
        FCS_FILES = files;

        global CMPDIR;
        CMPDIR = fileparts(outputdir);
    end

    % -------------------------------------------------------------------------
    [all_data, block_sizes, channel_names] = read_files(files);

    number_of_observations = size(all_data, 1);
    number_of_sources = numel(files);

    if wehave(samplesize)
        if samplesize > number_of_observations
            error('sample size exceeds number of observations');
        end

        indices = random_sample(samplesize, number_of_observations);
        data = all_data(indices, :);
    else
        indices = [];
        data = all_data;
    end

    clear('all_data');

    % -------------------------------------------------------------------------
    phenograph_table = run_phenograph(data);

    channels_table = make_channels_table(data);

    sources_table = make_sources_table(block_sizes, indices);

    tsne_table = run_tsne(data);

    % -------------------------------------------------------------------------
    all_tables = {sources_table, channels_table, tsne_table, phenograph_table};

    save_all(outputdir, relative_paths, all_tables, channel_names, savesession);
end
