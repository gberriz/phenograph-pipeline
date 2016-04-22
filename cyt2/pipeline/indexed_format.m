function format_ = indexed_format(prefix, max_value)
    format_ = sprintf('%s%%0%dd', prefix, floor(log10(max_value)) + 1);
end

