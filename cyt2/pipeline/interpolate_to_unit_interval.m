function normalized = interpolate_to_unit_interval(data)
    repmat_args = {size(data, 1), 1};

    min_ = min(data, [], 1);
    max_ = max(data, [], 1);

    base = repmat(min_, repmat_args{:});
    extent = repmat(max_ - min_, repmat_args{:});

    normalized = (data - base)./extent;
end
