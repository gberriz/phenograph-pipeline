function mapped_data = fast_tsne(data, no_dims, initial_dims, perplexity, theta)
%FAST_TSNE Runs the C++ implementation of Barnes-Hut t-SNE
%
%   mapped_data = fast_tsne(data, no_dims, initial_dims, perplexity, theta)
%
% Runs the C++ implementation of Barnes-Hut-SNE. The high-dimensional
% datapoints are specified in the NxD matrix data. The dimensionality of the
% datapoints is reduced to initial_dims dimensions using PCA (default = 50)
% before t-SNE is performed. Next, t-SNE reduces the points to no_dims
% dimensions. The perplexity of the input similarities may be specified
% through the perplexity variable (default = 30). The variable theta sets
% the trade-off parameter between speed and accuracy: theta = 0 corresponds
% to standard, slow t-SNE, while theta = 1 makes very crude approximations.
% Appropriate values for theta are between 0.1 and 0.7 (default = 0.5).
% The function returns the two-dimensional data points in mapped_data.
%
% NOTE: The function is designed to run on large (N > 5000) data sets. It
% may give poor performance on very small data sets (it is better to use a
% standard t-SNE implementation on such data).


% Copyright (c) 2014, Laurens van der Maaten (Delft University of Technology)
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
% 3. All advertising materials mentioning features or use of this software
%    must display the following acknowledgement:
%    This product includes software developed by the Delft University of Technology.
% 4. Neither the name of the Delft University of Technology nor the names of
%    its contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY LAURENS VAN DER MAATEN ''AS IS'' AND ANY EXPRESS
% OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
% EVENT SHALL LAURENS VAN DER MAATEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
% BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
% IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
% OF SUCH DAMAGE.


    if ~exist('no_dims', 'var') || isempty(no_dims)
        no_dims = 2;
    end
    if ~exist('initial_dims', 'var') || isempty(initial_dims)
        initial_dims = 50;
    end
    if ~exist('perplexity', 'var') || isempty(perplexity)
        perplexity = 30;
    end
    if ~exist('theta', 'var') || isempty(theta)
        theta = 0.5;
    end

    % Perform the initial dimensionality reduction using PCA

    % X = double(X);
    % X = bsxfun(@minus, X, mean(X, 1));

    % covX = X' * X;
    % [M, lambda] = eig(covX);
    % [~, ind] = sort(diag(lambda), 'descend');
    % if initial_dims > size(M, 2)
    %     initial_dims = size(M, 2);
    % end
    % M = M(:,ind(1:initial_dims));

    % X = X * M;
    % clear covX M lambda

    reduced_data = reduce_data(double(data), initial_dims);

    % Run the fast diffusion SNE implementation
    global DEBUG_REPRODUCIBILITY;
    if DEBUG_REPRODUCIBILITY
        global PRNG_SEED;
        global TSNE_PRNG_SEED_OFFSET;
        rand_seed = PRNG_SEED + TSNE_PRNG_SEED_OFFSET;
        write_data(reduced_data, no_dims, theta, perplexity, rand_seed);
    else
        write_data(reduced_data, no_dims, theta, perplexity);
    end

    bh_tsne_executable = find_bh_tsne();

    tic
    status = system(bh_tsne_executable);
    toc
    if status ~= 0
        error('execution of %s failed', bh_tsne_executable);
    end

    [mapped_data, ~, ~] = read_data;

    % landmarks = landmarks + 1;              % correct for Matlab indexing
    % delete('data.dat');
    delete('result.dat');
end


function reduced_data = reduce_data(data, dims)
    centered_data = bsxfun(@minus, data, mean(data, 1));

    covariance_matrix = centered_data' * centered_data;
    [change_of_basis_matrix, eigenvalues] = eig(covariance_matrix);
    new_dims = min(dims, size(change_of_basis_matrix, 2));

    [~, indices] = sort(diag(eigenvalues), 'descend');

    principal_components_matrix = change_of_basis_matrix(:, indices(1:new_dims));
    reduced_data = centered_data * principal_components_matrix;

    global INITIALDIMS;
    INITIALDIMS = new_dims;
    terse_warning('INITIALDIMS = %d', INITIALDIMS);
end


function [] = terse_warning(varargin)
    backtrace = warning('off', 'backtrace');
    verbose = warning('off', 'verbose');
    warning(varargin{:});
    warning('backtrace', backtrace.state);
    warning('verbose', verbose.state);
end


function bh_tsne_path = find_bh_tsne()
    stack = dbstack('-completenames');

    [bh_tsne_dir, ~, ~] = fileparts(stack(1).file);

    make_path = @(basename) fullfile(bh_tsne_dir, basename);

    default_basename = 'bh_tsne';

    candidate_basenames = { default_basename, ...
                            [default_basename '_' computer('arch')] };

    for i = 1:numel(candidate_basenames)

        candidate = make_path(candidate_basenames{i});
        terse_warning('Looking for bh_tsne executable: trying %s...', ...
                      candidate);

        if exist(candidate, 'file') == 2
            bh_tsne_path = candidate;
            terse_warning('Success');
            return;
        end
    end

    error('Found no recognizable bh_tsne executable in %s', bh_tsne_dir);
end


% Writes the datafile for the fast t-SNE implementation
function write_data(data, no_dims, theta, perplexity, rand_seed)

    [number_of_columns, number_of_rows] = size(data);

    outstream = fopen('data.dat', 'wb');

    fwrite(outstream, number_of_columns, 'integer*4');
    fwrite(outstream, number_of_rows, 'integer*4');
    fwrite(outstream, theta, 'double');
    fwrite(outstream, perplexity, 'double');
    fwrite(outstream, no_dims, 'integer*4');
    fwrite(outstream, data', 'double');

    if exist('rand_seed', 'var') && ~isempty(rand_seed)
        fwrite(outstream, rand_seed, 'integer*4');
    end

    fclose(outstream);
end


% Reads the result file from the fast t-SNE implementation
function [data, landmarks, costs] = read_data
    instream = fopen('result.dat', 'rb');

    number_of_columns = fread(instream, 1, 'integer*4');
    number_of_rows = fread(instream, 1, 'integer*4');

    number_of_elements = number_of_columns * number_of_rows;
    shape = [number_of_rows number_of_columns];

    data = reshape(fread(instream, number_of_elements, 'double'), shape)';

    landmarks = fread(instream, number_of_columns, 'integer*4');
    costs = fread(instream, number_of_columns, 'double'); % this vector contains
                                                          % only zeros
    fclose(instream);
end
