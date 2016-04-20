function [] = save_to_tsv(path_, data_as_table)
    [dirname, ~, ~] = fileparts(path_);
    maybe_create_dir(dirname);
    writetable(data_as_table, path_, 'Delimiter', '\t', 'FileType', 'text');
end
