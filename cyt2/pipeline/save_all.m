function [] = save_all(outputdir, relative_paths, all_tables, channel_names, ...
                       savesession)

    % -------------------------------------------------------------------------

    [sources_table, channels_table, tsne_table, phenograph_table] = ...
        deal(all_tables{:});

    % -------------------------------------------------------------------------

    number_of_sources = numel(relative_paths);

    % -------------------------------------------------------------------------

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

        global DEBUG_REPRODUCIBILITY;
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
        normalized_data = interpolate_to_unit_interval(data);
        normalized_table = replace_channel_columns(table_, normalized_data);
        normalized_table.GroupCount = ...
            100 * normalized_table.GroupCount/sum(normalized_table.GroupCount);
        normalized_table.Properties.VariableNames{'GroupCount'} = 'percentage';
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
            global DEBUG_REPRODUCIBILITY;
            if bool(DEBUG_REPRODUCIBILITY) || isempty(channel_names)
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
% table utils

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
