function [] = bot(inputdir, outputdir, savesession)
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
%    BOT(INPUTDIR, OUTPUTDIR, true), in addition, results in having the current
%    MATLAB session saved in the file matlab_session.mat, under OUTPUTDIR.



    if ~exist('savesession', 'var')
        bot(inputdir, outputdir, false);
        return;
    end

    tic
    run_pipeline(inputdir, outputdir, savesession);
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

% -----------------------------------------------------------------------------
% output

function [] = maybe_create_dir(dir_)
    [status, message, ~] = mkdir(dir_);
    if status == 0 % sic
        error(message);
    end
end

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

function sources_table = make_sources_table(block_sizes, sample_indices)
    block_indices = 1:numel(block_sizes);
    % letters = 'A':'Z';
    % block_values = letters(block_indices);

    observation_index_to_source_index = make_mapper(block_sizes, block_indices);

    source_indices = arrayfun(@(i) observation_index_to_source_index(i), ...
                            sample_indices);

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
% data normalization

function normalized = normalize_(data)
    repmat_args = {size(data, 1), 1};

    min_ = min(data, [], 1);
    max_ = max(data, [], 1);

    base = repmat(min_, repmat_args{:});
    extent = repmat(max_ - min_, repmat_args{:});

    normalized = (data - base)./extent;
end

% -----------------------------------------------------------------------------
% table utils

function table_ = array_to_table(array, headers)
    table_ = array2table(array, 'VariableNames', headers);
end

function [valid_names, modified] = make_valid_names(names)
    whitespace_free_names = cellfun(@(s) regexprep(s, '\s+', '?'), names, ...
                                    'UniformOutput', false);
    [valid_names, modified] = matlab.lang.makeValidName(whitespace_free_names);

    if ~any(modified); return; end

    % -------------------------------------------------------------------------

    formatter = @(i) sprintf('    ''%s'' -> %s\n', ...
                             names{i}, valid_names{i});

    key_value_pairs = arrayfun(formatter, find(modified), ...
                               'UniformOutput', false);

    message = sprintf(['The following names were converted ' ...
                       'to valid MATLAB identifiers:\n%s'], ...
                      sprintf('%s', key_value_pairs{:}));

    terse_warning(message);
end

function data_as_table = cell_to_table(data_as_cell)
    header_row = data_as_cell(1, :);
    data_rows = data_as_cell(2:end, :);
    data_as_table = ...
        cell2table(data_rows, 'VariableNames', make_valid_names(header_row));
end

% -----------------------------------------------------------------------------
% filename manipulation

function new_path = change_extension(path_, extension)
    [extfree, ~] = splitext(path_);
    new_path = [extfree, extension];
end

function [extfree, ext] = splitext(path_)
    [dirname, basename, ext] = fileparts(path_);
    extfree = fullfile(dirname, basename);
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
% warnings

function [] = terse_warning(varargin)
    backtrace = warning('off', 'backtrace');
    verbose = warning('off', 'verbose');
    warning(varargin{:});
    warning('backtrace', backtrace.state);
    warning('verbose', verbose.state);
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

function [] = run_pipeline(inputdir, outputdir, savesession)

    global DEV_MODE;
    global DEBUG_REPRODUCIBILITY;
    global PRNG_SEED;
    if DEBUG_REPRODUCIBILITY
        if isempty(PRNG_SEED)
            global DEFAULT_PRNG_SEED;
            if isempty(DEFAULT_PRNG_SEED)
                error('no PRNG seed available');
            else
                PRNG_SEED=DEFAULT_PRNG_SEED;
            end
        end
        terse_warning('PRNG_SEED: %d', PRNG_SEED);
    else
        clear('global PRNG_SEED');
    end

    % -------------------------------------------------------------------------

    relative_paths = get_files(inputdir);
    if DEV_MODE
        relative_paths = relative_paths([1 2 21 22]);
        % relative_paths = relative_paths([1 21]);
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

    if DEV_MODE
        % sample_size = 11 * number_of_sources;
        sample_size = 100 * number_of_sources;
    else
        sample_size = 100000;
    end

    if DEBUG_REPRODUCIBILITY
        global SAMPLE_SIZE;
        SAMPLE_SIZE = sample_size;
    end

    if sample_size > number_of_observations
        error('sample size exceeds number of observations');
    end

    sample_indices = random_sample(sample_size, number_of_observations);

    sample = all_data(sample_indices, :);
    clear('all_data');

    % -------------------------------------------------------------------------
    channels_table = make_channels_table(sample);

    sources_table = make_sources_table(block_sizes, sample_indices);

    tsne_table = run_tsne(sample);

    phenograph_table = run_phenograph(sample);

    results_table = [sources_table phenograph_table channels_table tsne_table];

    % -------------------------------------------------------------------------
    channel_columns = channels_table.Properties.VariableNames;

    function means_table = cluster_means(table_)
        stat = 'mean';
        keyvars = {'source', 'cluster'};

        temp_table = grpstats(table_, keyvars, stat);
        temp_table.Properties.RowNames = {};

        % ---------------------------------------------------------------------
        % restore the names of the averaged columns
        for old_name_cell = table_.Properties.VariableNames
            old_name = old_name_cell{:};
            if any(strcmp(keyvars, old_name)); continue; end
            new_name = [stat '_' old_name];
            temp_table.Properties.VariableNames{new_name} = old_name;
        end

        if DEBUG_REPRODUCIBILITY
            means_table = sortrows(temp_table, 'cluster', 'ascend');
            return;
        end

        % ---------------------------------------------------------------------
        % sort table descendingly by the mean signal across all the channels
        temp_table.mean_signal = ...
            mean(table2array(temp_table(:, channel_columns)), 2);

        means_table = sortrows(temp_table, 'mean_signal', 'descend');
        means_table.mean_signal = []; % premature optimization?
    end

    function new_table = replace_channel_columns(table_, data)
        new_table = table_;
        new_table{:, channel_columns} = data;
    end

    means_table = cluster_means(results_table);
    means = table2array(means_table(:, channel_columns));
    truncated_means = max(means, 0);
    truncated_means_table = replace_channel_columns(means_table, ...
                                                    truncated_means);

    % -------------------------------------------------------------------------
    function normalized_table = normalize_table(table_)
        data = table2array(table_(:, channel_columns));
        normalized_data = normalize_(data);
        normalized_table = replace_channel_columns(table_, normalized_data);
        normalized_table.GroupCount = ...
            100 * normalized_table.GroupCount/sum(normalized_table.GroupCount);
        normalized_table.Properties.VariableNames{'GroupCount'} = 'percentage';
    end

    function [] = save_to_csv(path_to_tsv, table_)

        outputpath = change_extension(path_to_tsv, '.csv');

        [dirname, ~, ~] = fileparts(outputpath);
        maybe_create_dir(dirname);

        % ---------------------------------------------------------------------

        subtable = table_(:, [channel_columns 'percentage' 'cluster']);
        subtable.cluster = uint64(subtable.cluster);
        headers = subtable.Properties.VariableNames;
        rows = table2cell(subtable);
        data = cat(1, headers, rows);

        % ---------------------------------------------------------------------

        cell2csv(outputpath, data);

    end

    function renamed_table = rename_channel_columns(table, channel_names_i)
        old_names = channels_table.Properties.VariableNames;
        new_names = make_valid_names(channel_names_i);
        renamed_table = table;
        for i = 1:numel(new_names)
            renamed_table.Properties.VariableNames{old_names{i}} = new_names{i};
        end
    end

    function [] = save_normalized_tables(subdir, table_)

        subpath = fullfile(outputdir, subdir);

        for i = categorical(1:number_of_sources)

            subtable = table_(table_.source == i, :);
            subtable.source = [];
            normalized_subtable = normalize_table(subtable);

            % -----------------------------------------------------------------
            relative_path = change_extension(relative_paths{i}, '.tsv');
            outputpath = fullfile(subpath, relative_path);

            % -----------------------------------------------------------------
            if DEBUG_REPRODUCIBILITY
                subtable_to_save = normalized_subtable;
                wanted_columns = [channel_columns 'percentage' 'cluster'];
            else
                subtable_to_save = ...
                    rename_channel_columns(normalized_subtable, ...
                                           channel_names{i});
                wanted_columns = ['cluster' 'percentage' channel_columns ...
                                  tsne_table.Properties.VariableNames];
            end

            save_to_tsv(outputpath, subtable_to_save(:, wanted_columns));
        end

    end

    maybe_create_dir(outputdir);

    save_normalized_tables('full',                means_table);
    save_normalized_tables('truncated', truncated_means_table);

    % -------------------------------------------------------------------------
    if savesession
        sessionfile = fullfile(outputdir, 'matlab_session.mat');
        save(sessionfile, '-v7.3');
    end
end
