function [c,Q,bestpartition,bestpartitionhierarchy] = LouvainfromBin( filename, numiters )
%the path to louvain in each cyt

    [curr_path, ~, ~] = fileparts(mfilename('fullpath'));
    curr_path = [curr_path filesep];
    ps = [curr_path 'Louvain' filesep];

    filename = strrep( filename, '.bin', '' );
    % begin
    fprintf(1, 'MATLAB: calling convert:\n');
    command = [ps 'convert -i ' filename '.bin -o ' filename '_graph.bin -w ' filename '_graph.weights' ];
    run_command( command );

    global DEBUG_REPRODUCIBILITY;
    if DEBUG_REPRODUCIBILITY
        global PRNG_SEED;
        global PHENOGRAPH_PRNG_SEED_OFFSET;
        rand_seed = PRNG_SEED + PHENOGRAPH_PRNG_SEED_OFFSET;
        rng(rand_seed);
        have_seeded_prng = true;
    else
        have_seeded_prng = false;
    end

    seedopt = '';

    % run community detection
    for iter = 1:numiters

        if have_seeded_prng
            seedopt = sprintf('-s %u', randi(2^32) - 1);
        end

        fprintf(1,'MATLAB: running community detection, ITERATION %i\n', iter );
        command = sprintf(['%scommunity "%s_graph.bin" -l -1 -v ' ...
                           '-w "%s_graph.weights" %s > "%s.tree"'], ...
                          ps, filename, filename, seedopt, filename);

        r = run_command( command );
        fprintf(1, '%s\n', r);

        try
            % find each iteration's modularity
            q = find_modularity( r );
            fprintf( 1, 'MATLAB: modularity scores:\n' );
        catch
            cleanup(filename)
            error('Unable to find modularity score in the stderr: %s.\nCheck that the correct path to the Louvain code is specified.', r);
        end


        % find number of lvevls
        command = [ps 'hierarchy ' filename '.tree' ];
        r = run_command( command );
        fprintf(1, '\n' );

        r = strtok(r, 10);
        r = regexprep( r, 'Number of levels: ', '' );
        num_levels = str2double( r )-1;

        fprintf( 1, 'MATLAB: max level is %d\n', num_levels );

        % import levels
        for level = 1:num_levels
            fprintf( 1, 'MATLAB: importing level %d\n', level );
            command = [ps 'hierarchy ' filename '.tree -l ' num2str( level ) ' > ' filename '.tmp' ];
            run_command( command );
            hierarchy_output = load( [filename '.tmp'] );
            c{iter,level} = hierarchy_output(:,2) + 1;
            Q{iter,level} = q(level);
        end

    end
    % find best partition
    maxmod = 0;
    for i = 1:numel(Q)
        if Q{i} > maxmod
            maxmod = Q{i};
            [I,J] = ind2sub( size(Q), i );
        end
    end
    bestpartition = c{I,J};
    bestpartitionhierarchy = c(I,:);
    % delete temporary files
    cleanup(filename)
end

%-------------------------
function Q = find_modularity( r )
% Q = find_modularity( r )
% convert the text output into modularity score of each iteration
    signature = '  modularity increased from %f to %f';
    idx = 0;
    while( ~isempty( r ) )
        % read a line and match it to the signature
        [token, r] = strtok( r, char( 10 ) );
        a = sscanf( token, signature );

        if( ~isempty( a ) )
            % signature matched copy new modularity
            idx = idx + 1;
            Q( idx ) = a( 2 );
        end
    end
end

%-------------------------
function cleanup(filename)
    files{1} = dir([filename '*.tmp']);
    files{2} = dir([filename '*.tree']);
    files{3} = dir([filename '*.weights']);
    files{4} = dir([filename '*.bin']);
    for i = 1:4
        for j = 1:length(files{i})
            delete( files{i}(j).name )
        end
    end
end

function output = run_command(command)
    fprintf(1, '%s\n', command);
    [status, output] = system(command);
    if status ~= 0
        error('"%s" failed with status %d:\n%s', command, status, output)
    end
end
