function [graph_traj, lnn] = wanderlust( data, distance, k, l, num_graphs, s, num_landmarks, verbose, partial_order)
% [graph_traj, lnn] = wanderlust( data, distance, k, l, num_graphs, s, num_landmarks, power_k, verbose )
%
% run the wonderlust algorithm on data. generate num_graphs klNN graphs; starting point is s, choose num_landmarks
% random landmarks. weigh landmark-node pairs by distance to the power of power_k.
%
% (use distance as distance metric)
%
% (alternatively data can be lnn graph)
%
% TODO fix initial comment
snn=false;
% build lNN graph
if( issparse( data ) )
	if( verbose )
		disp 'using prebuilt lNN graph';
	end
	lnn = data;
else
	if( verbose )
		disp 'building lNN graph';
    end

    tic
	lnn = parfor_spdists_knngraph( data, l, 'distance', distance, 'chunk_size', 2000, 'SNN', true, 'verbose', true );
    
%     if (snn)
%         lnn = knn2jaccard(lnn);
%     end
    sprintf('lnn computed: %gs', toc/60);
end

% generate klNN graphs and iteratively refine a trajectory in each
for graph_iter = 1:num_graphs
	if( verbose )
		fprintf( 1, 'iter #%d:\n', graph_iter );
	end

	iter_t = tic;

	% randomly generate a klNN graph
	if( verbose )
		fprintf( 1, 'entering knn graph: ' );
	end

	klnn = spdists_klnn( lnn, k, verbose );
 	klnn = spdists_undirected( klnn );
 
	if( verbose )
		fprintf( 1, ' done (%3.2fs)\n', toc( iter_t ) );
		fprintf( 1, 'entering trajectory landmarks: ' );
	end

	% run traj. landmarks
	[ traj, dist, iter_l ] = trajectory_landmarks( klnn, s, num_landmarks, partial_order, verbose, data, distance );

	if( verbose )
		fprintf( 1, ' done (%3.2fs)...\n', toc( iter_t ) );
    end
 
	% calculate weighed trajectory
    sdv = mean ( std ( dist) )*3;
	dist_power = dist;
%     dist_power(:, :) = 1;
%     dist_power = ((2*pi*(sdv^2))^(-1)) * exp( -.5 * (dist_power / sdv).^2);
	dist_power_norm = dist_power ./ repmat( sum( dist_power ), size( dist_power, 1 ), 1 );
     dist_power_norm = 1-dist_power_norm;
%     try
%     str_input = input(sprintf('plot the first 2 features?'),'s');
% 
%     if strcmpi(str_input, 'y') || strcmpi(str_input, 'yes') || strcmpi(str_input, 'yea')
%         [~, paths, ~] = graphshortestpath( klnn, iter_l(1), 'METHOD','Dijkstra');
%         col1=1;
%         col2 =2;
%         set(0,'DefaultLineMarkerSize',2);
%         myplotclr(data(:, col1), data(:, col2), dist_power_norm_s(randperm(length(dist_power_norm_s))), dist_power_norm_s, 'o', colormap,[min(dist_power_norm_s), max(dist_power_norm_s)], false);
%         hold on;
%         set(0,'DefaultLineMarkerSize',18);
%         myplotclr(data(iter_l, col1), data(iter_l, col2), ones(length(iter_l), 1)*(max(dist_power_norm_s)), ones(length(iter_l), 1)*(max(dist_power_norm_s) ), 'X', colormap,[min(dist_power_norm_s), max(dist_power_norm_s)], false);
%         pathi = 54734;
%         pathi = iter_l(1);
%         hold on;
%         plot(data(paths{pathi}, 1),data(paths{pathi}, 2), '-', 'Color', 'blue', 'MarkerSize', 18);
%     end
%     catch
%         fprintf( 1, 'thrown out');
%     end
%     
	t( 1, : ) = sum( traj .* dist_power_norm );
    if any(isnan(t))
            strPrompt = sprintf('there are nans in t, should I stop?');
            str_input = input(strPrompt,'s');
            if strcmpi(str_input, 'y') || strcmpi(str_input, 'yes') || strcmpi(str_input, 'yea')
                data(find(traj(1, :) == inf), :)
                return;
            end 
    end
    
	% iteratively realign trajectory (because landmarks moved)
	converged = 0; user_break = 0; realign_iter = 1;
	while( ~converged && ~user_break)
		realign_iter = realign_iter + 1;

		traj = dist;
		for idx = 1:size( dist, 1 )
			% find position of landmark in previous iteration
			idx_val = t( realign_iter - 1, iter_l( idx ) );
			% convert all cells before starting point to the negative
			before_indices = find( t( realign_iter - 1, : ) < idx_val );
			traj( idx, before_indices ) = -dist( idx, before_indices );
			% set zero to position of starting point
			traj( idx, : ) = traj( idx, : ) + idx_val;
		end

		% calculate weighed trajectory
		t( realign_iter, : ) = sum( traj .* dist_power_norm );

		% check for convergence
		converged = corr( t( realign_iter, : )', t( realign_iter - 1, : )' ) > 0.99;
        
        if (mod(realign_iter,1000)==0)
            mean ( std ( dist) )
            t(realign_iter-4:realign_iter, :);
            last_corr = corr( t( realign_iter, : )', t( realign_iter - 1, : )' );
            strPrompt = sprintf('iter=%g, corr=%g. should I stop (y\n)?', realign_iter, last_corr);
            str_input = input(strPrompt,'s');
            if strcmpi(str_input, 'y') || strcmpi(str_input, 'yes') || strcmpi(str_input, 'yea')
                user_break = true;
            end
        end
	end
	fprintf( 1, '%d realignment iterations, ', realign_iter );

	% save final trajectory for this graph
	graph_traj( graph_iter, : ) = t( realign_iter, : );

	if( verbose )
		toc( iter_t );

		fprintf( 1, '\n' );
	end
end
end


	function spdists = spdists_klnn( spdists, k, verbose )
	% spdists = spdists_klnn( spdists, k, verbose )
	%
	% given a lNN graph spdists, choose k neighbors randomly out of l for each node

	remove_edges = [];

	for idx = 1:length( spdists )
		% remove l-k neighbors at random
		neighs = find( spdists( :, idx ) );
		l = length( neighs ); % count number of neighbors
		remove_indices = neighs( randsample( length( neighs ), l - k ) );
		idx_remove_edges = sub2ind( size( spdists ), remove_indices, ones( l - k, 1 ) * idx );
		remove_edges = [ remove_edges; idx_remove_edges ];

		if( verbose )
			if( mod( idx, 50000 ) == 0 )
				fprintf( 1, '%3.2f%%', idx / length( spdists ) * 100 );
			elseif( mod( idx, 10000 ) == 0 )
				fprintf( 1, '.' );
			end
		end
	end

	spdists( remove_edges ) = 0;
    end

	function [ traj, dist, l ] = trajectory_landmarks( spdists, s, n, partial_order, verbose,data, distance )
	% [ traj, dist, l ] = trajectory_landmarks( spdists, s, n, verbose )
	%
	% calculate the trajectory score of each point in spdists.
	%
	% s: list of indices of possible starting points. one of these points will be used to generate a reference
	% trajectory; the landmark shortest paths will be aligned to this reference.
	% n: list of landmark indices to use; or, alternatively, the number of landmarks to choose randomly from all
	% points.
	%
	% traj is a |n|x|spdists| matrix, where row i is the aligned shortest path from landmark i to each other point.
	% dist is a |n|x|spdists| matrix, where row i is the shortest path from landmark i to each other point. l is
	% the list of landmarks, l(1) is the starting point.

	if( length( s ) > 1 )
		% if given a list of possible starting points, choose one
		s = randsample( s, 1 );
	end

	if( length( n ) == 1 )
 		% if not given landmarks list, decide on random landmarks
 		n = randsample( length( spdists ), n - 1 - length(partial_order) );
        
        % flock centers once
        [IDX, ~] = knnsearch(data, data(n, :), 'distance', distance, 'K', 10);     
        for i=1:numel(n)
            n(i) = knnsearch(data, median(data(IDX(i, :), :)), 'distance', distance); 
        end
        
%         % shift random list to centroids
%         [~, C] = kmeans(data, n - 1 - length(partial_order));
%         n = knnsearch(data, C, 'distance', distance);     
    end

    partial_order = [s;partial_order(:)]; % partial_order includes start point
	l = [ partial_order; n ]; % add extra landmarks if user specified
    
    % prune out weak edges 
%     res = 256;
% 	[band, density, x, y] = kde2d(data, res);
    
        % maybe remove edges jumping over low density regions?
        % BW1 = edge(density,'canny', .001);
        % imshow(BW1)
        % set(gca,'YDir','normal')
        
        %remove bad edges if the path to a landmark has a suspicious jump
%     mapX = round( interp1(diag(x), 1:res, data(:, 1), 'linear', 'extrap') );
%     mapY = round( interp1(diag(y), 1:res, data(:, 2), 'linear', 'extrap') );
%     recalc = true;
%     cou = 0;
%     while recalc && cou < 20
%         cou = cou+1;
%         [dist( 1, : ), paths, ~] = graphshortestpath( spdists, s);%, 'directed', false );
%         recalc = false;
%         for pathidx=2:numel(l) %iterate over path to each landmark
%             l_path = paths{l(pathidx)};
%             path_jumpsX = abs(mapX(l_path(2:end))-mapX(l_path(1:end-1)));
%             path_jumpsY = abs(mapY(l_path(2:end))-mapY(l_path(1:end-1)));
%             path_jumps = path_jumpsX + path_jumpsY;
%             bad_nodes = find(path_jumps > (mean(path_jumps) + 4*std(path_jumps)));
%             if any(bad_nodes)
%                 disp(sprintf('removing %g bad connections\n', numel(bad_nodes)));
%                 spdists(sub2ind(size(spdists), l_path(bad_nodes), l_path(bad_nodes+1))) = 0;
%                 spdists(sub2ind(size(spdists), l_path(bad_nodes+1), l_path(bad_nodes))) = 0;
%                 recalc = true;
%             end
%         end
%     end
    

	% calculate all shortest paths
	for idx = 1:length( l )
		[dist( idx, : ), paths, ~] = graphshortestpath( spdists, l( idx ),'METHOD','Dijkstra');%, 'directed', false );
		unreachable = find(dist(idx,:)==inf);
        dist(idx,unreachable) = max(dist(idx, dist(idx,:)<inf));
		if( verbose )
			fprintf( 1, '.' );
        end
	end
    
    % adjust paths according to partial order by redirecting
    nPartialOrder = length(partial_order);
    for radius = 1:nPartialOrder 
        for landmark_row = 1:nPartialOrder
            if (landmark_row + radius <= nPartialOrder)
                a = landmark_row;
                b = landmark_row + (radius-1);
                c = landmark_row + radius;
                dist(a, partial_order(c)) = dist(a, partial_order(b)) + dist(b, partial_order(c));
            end
            if (landmark_row - radius >= 1)
                a = landmark_row;
                b = landmark_row - (radius-1);
                c = landmark_row - radius;
                dist(a, partial_order(c)) = dist(a, partial_order(b)) + dist(b, partial_order(c));
            end
        end
    end

	% align to dist_1
	traj = dist;
    for idx = 2:length(partial_order)
        [~, closest_landmark_row] = min(dist); %closest landmark will determine directionality
        traj(idx, closest_landmark_row < idx) = -dist(idx, closest_landmark_row < idx);
        traj( idx, : ) = traj( idx, : ) + dist( 1, l( idx ) );
    end
    
    if length( l ) > length(partial_order)
        for idx = length(partial_order)+1:length( l )
            % find position of landmark in dist_1
            idx_val = dist( 1, l( idx ) );
            % convert all cells before starting point to the negative
            before_indices = find( dist( 1, : ) < idx_val );
            traj( idx, before_indices ) = -dist( idx, before_indices );
            % set zero to position of starting point
            traj( idx, : ) = traj( idx, : ) + idx_val;
        end
    end
end


