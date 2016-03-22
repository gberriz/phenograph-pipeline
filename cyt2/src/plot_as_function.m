function plot_as_function(x, Y, varargin)
% plots Y (NxM) as a function of vector x (N length).  
% an average of points in Y is computed for 'num_locs' along vector x.
%
% 'labels': a string cell array to enter for the legend
%
% specify 'avg_type': default 'gaussian_var'
%
% 'sliding': a fixed averaged is computed within a sliding window. takes
% into account points within two windows to each side. so let's say 100
% locations would mean each location is averaging points from two windows to
% the right and to the left
%
% 'linear': same as sliding the each point's contibution dimishes linearly
% 
% 'squared': same as sliding the each point's contibution dimishes
% polynomialy
%
% 'gaussian': all points contibute to the weighted average but beyond 2
% standard deviations the contibution is insignificant.
%
% 'gaussian_var': same as gaussian but the standard deviation is variable
% as a function of the normalized density of the points along x. so in
% saprser locations on x, we use a larger standard deviation to let farther
% points contribute.
%
%
% 'show_error': true or false, shows error bars along the curve, default
% 'false'
% 
% Michelle Tadmor, Columbia University, 2013

clear persistent_ksdensity;
legend_flag = false;
avg_type = 'gaussian_var';
num_locs = 100;
show_error = false;
normalize = true;
rank = false;
make_distribution_movie = false;
smoothness_factor = 0.5;
svGolay = true;
control_density = false;
matPatchColors = [0.75 0.75 0.75; 0.6 0.6 0.6; 0 0 1];
branch = zeros(1, numel(x));

for i=1:length(varargin)-1
    if(strcmp(varargin{i},'num_locs'))
        num_locs = varargin{i+1};
    elseif(strcmp(varargin{i},'avg_type'))
        avg_type = varargin{i+1};
    elseif(strcmp(varargin{i},'labels'))
        legend_flag = true;
        labels = varargin{i+1};
    elseif(strcmp(varargin{i},'show_error'))
        show_error = varargin{i+1};  
    elseif(strcmp(varargin{i},'normalize'))
        normalize = varargin{i+1};
    elseif(strcmp(varargin{i},'rank'))
        rank = varargin{i+1};
        if rank
            control_density = false;
        end
    elseif(strcmp(varargin{i},'svGolay'))
        svGolay = varargin{i+1};
    elseif(strcmp(varargin{i},'smooth'))
        smoothness_factor = varargin{i+1};
    elseif(strcmp(varargin{i},'branch'))
        branch = varargin{i+1};
    elseif(strcmp(varargin{i},'movie'))
        make_distribution_movie = true;
        movie_name = varargin{i+1};
    end
end

if (rank)
    x = tiedrank(x);
end    

weights = zeros(num_locs, length(x));
Y_vals = zeros(num_locs, size(Y, 2));
Y_errs = zeros(num_locs, size(Y, 2));
Yn = size(Y, 1);

% compute a weight for each value (data point), at each plot location
for i=1:num_locs
	weights(i, :) = compute_weights(x, (i/num_locs)*range(x)+min(x), avg_type, smoothness_factor);
end

%clean branch
if any(branch)
    pre_cleaning = weights;

    weights_by_branch = sum(weights(:, branch==0), 2);
    weights_by_branch(:, end+1) = sum(weights(:, branch~=0), 2)
    branching_loc = find(weights_by_branch(:, 1)<weights_by_branch(:, 2));
    weights(branching_loc(1):num_locs, branch==0) = 0;
    weights(1:branching_loc(1)-1, branch~=0) = 0;
   
end

real_weights = weights;
for bri=unique(branch)'
    if any(branch)
        weights = real_weights;
        weights(:, branch~=bri) = 0;
        visible_locs = find(sum(weights, 2));
        if (bri ~= 0) % if on a branch
            trans_length = ceil(num_locs/10);
            for transi = 1:(trans_length-1)
                weights(visible_locs(transi), branch~=bri) = (1-transi/trans_length)*pre_cleaning(visible_locs(transi), branch~=bri);
            end
        else
            weights(visible_locs(end)+1, :) = pre_cleaning(visible_locs(end)+1, :);
        end
    end
% Compute weighted averages at each location
X = linspace(min(x), max(x), num_locs);
Y_vals = weights*Y./repmat(sum(weights, 2), 1, size(Y, 2));

if control_density
    dens = sum(weights, 2)';
    dens = dens./max(dens);
    dens = floor(dens*10);
    X = linspace(min(x), max(x), num_locs+sum(dens));
    Yrow = 0;
	Y_vals_stretched = zeros(0, size(Y_vals, 2));
    for densi=1:numel(dens)
        Yrow = Yrow+1;
        copies = 0;
        while copies<=dens(densi)
            Y_vals_stretched(end+1, :) = Y_vals(Yrow, :);
            copies = copies+1;
        end
    end
    Y_vals = Y_vals_stretched;
	
    % smooth the streching
    for col=1:size(Y_vals, 2)
        Y_vals(:, col) = smooth(X, Y_vals(:, col),sqrt(num_locs*2), 'sgolay');
    end          

end

if (svGolay)  
    for col=1:size(Y_vals, 2)
        Y_vals(:, col) = smooth(X, Y_vals(:, col),sqrt(num_locs*2), 'sgolay');
    end          
end

Y_vals_orig = Y_vals;
if (normalize)
    % we want to normalize while accounting for branches
    y_vals_all = [];
    for bri=unique(branch)'
        weights = real_weights;
        weights(:, branch~=bri) = 0;
        % Compute weighted averages at each location
        y_vals_all = [y_vals_all; weights*Y./repmat(sum(weights, 2), 1, size(Y, 2))];
    end
    Y_vals = Y_vals-repmat(prctile(y_vals_all, 1, 1), size(Y_vals,1),1);
    y_vals_all = y_vals_all-repmat(prctile(y_vals_all, 1, 1), size(y_vals_all,1),1);
    Y_vals = Y_vals./repmat(prctile((y_vals_all), 97, 1),size(Y_vals,1),1);
end

matColors = distinguishable_colors(size(Y, 2));
set(gca, 'ColorOrder', matColors);
set(0, 'DefaultAxesColorOrder', matColors);

plot(X, Y_vals(:, 1), 'Color', matColors(1, :));    
if (size(Y, 2)> 1)
    for col=2:size(Y, 2)
        hold on;
        plot(X, Y_vals(:, col), 'Color', matColors(col, :));        
    end
end
if (show_error)

    for i=1:num_locs       
        %symmetrical for error bars
        Y_errs(i, :) = weights(i, :)*(((Y-repmat(Y_vals_orig(i, :),Yn,1)-repmat(prctile(Y_vals, 1, 1), size(Y,1),1))./...
            repmat(prctile((Y_vals_orig), 99, 1),size(Y,1),1))).^2/sum(weights(i, :));        
    end
    
    for yi=1:size(Y, 2)
        % plot light grey background first for variance
        fill([X, fliplr(X)], [(Y_vals(:, yi)-Y_errs(:, yi))', fliplr((Y_vals(:, yi)+Y_errs(:, yi))')],...
            matColors(yi,:),'linestyle','none','facealpha',.5);
    end
end

if (any(branch))
    hold on;
end
end
if ~(rank || control_density)    
    try
    dens = sum(weights, 2)';
    ca = axis;
    hold on;
    imagesc(X, ca(3):0.04:ca(3)+.05, dens, [0, max(dens)]);
    colorbar;
    axis(ca);
    catch e
    disp(getReport(e,'extended'));
    end
end

if (legend_flag)
    legend(remove_repeating_strings(labels), 'Interpreter', 'none');
end

if (make_distribution_movie && askuser('Are you sure you''d like to make a movie?'));
    
    distances_per_window = cell(num_locs);
    for i=1:num_locs
        distances_per_window{i} = pdist2(Y(weights(i, :) > .4, :), Y(weights(i, :) > .4, :), 'cosine');
    end

    current_figure = gcf;
    try 
    fid = figure;
    hold on;
    [f,xi] = ksdensity(distances_per_window{1}(:),'npoints',200);
    hLine_new = plot(xi, f);
    title(sprintf('%.2f: %g points', (1/num_locs), sum(weights(1, :) > .4)));

    % Set up the movie.
    writerObj = VideoWriter(sprintf('%s.avi', movie_name)); % Name it.
    writerObj.FrameRate = 10;           % How many frames per second.

    open(writerObj); 

    for i=2:num_locs
        hLine = hLine_new;

        figure(fid); % Makes sure you use your desired frame.
        if (isempty(distances_per_window{i}(:)))
            xi = 0;
            f = 0;
        else
            [f,xi] = ksdensity(distances_per_window{i}(:),'npoints',250);
        end
        hLine_new = plot(xi, f);
        title(sprintf('%.2f: %g points', (i/num_locs), sum(weights(i, :) > .4)));
        delete(hLine);


        %if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
            frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
            writeVideo(writerObj, frame);
        %end

    end
    hold off
    close(writerObj); % Saves the movie.
    
    catch e
        % Return to base figure
        set(0,'CurrentFigure', current_figure);
        disp (getReport(e));
    end
end

end

function weights = compute_weights(points, loc, type, factor)
    
    range = quantile(points, .98) - quantile(points, .02);
    min_std_dev = factor*.18*range; % minimum std_dev for dense regions
    max_std_dev = .19*range; % max std_dev for sparse regions  
    linear_slope = 10/range;
    
    if strcmpi('sliding', type) %set '1's on the indices in the windows 
        weights = (points < (loc + 2*min_std_dev)) & ...
                  (points > (loc - 2*min_std_dev));
    
    elseif strcmpi('linear', type)
        weights = 1 - linear_slope*(abs(points - loc));
        weights(weights<0) = 0;
        weights = weights/sum(weights);       
    
    elseif strcmpi('squared', type)
        weights = 1 - ((linear_slope*(points - loc)).^2);
        weights(weights<0) = 0;
        weights = weights/sum(weights);
             
    elseif strcmpi('gaussian_var', type)
        [f, xi] = persistent_ksdensity(points);
    	d = f(minind(abs(xi-loc)))/max(f);

        std_dev = d*(min_std_dev)+(1-d)*max_std_dev;

        weights = ((2*pi*(std_dev^2))^(-1))*exp(-.5*((points - loc)/std_dev).^2);
    
    else % default is strcmpi('gaussian', type)
        weights = ((2*pi*(min_std_dev)^2)^(-1))*exp(-.5*((points - loc)/min_std_dev).^2);
    end
% 	weights = mynormalize(weights(:), 100);
end

function ind = minind(x)
    [~, ind] = min(x);
end


function [f, xi] = persistent_ksdensity(points)

    persistent f_p;
    persistent xi_p;    

    if isempty(f_p)
        [f_p, xi_p] = ksdensity(points);
    end
    f = f_p;
    xi = xi_p;
end