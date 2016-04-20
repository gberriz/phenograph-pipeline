function [] = maybe_create_dir(dir_)
    [status, message, ~] = mkdir(dir_);
    if status == 0 % sic
        error(message);
    end
end
