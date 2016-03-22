function scatter_by_point(x, y, gcolors, dot_size)
        
    groups = unique(gcolors);
    nGroups = numel(groups);
    
    map = distinguishable_colors(nGroups);
    
    lmap = length(map); 
    lgcolors = length(unique(gcolors));
    len_add = lgcolors - lmap;
    if (len_add>0)
        newMap=map(randsample(1:lmap,lgcolors-lmap,true),:);
        map(lmap+1:lgcolors,:)=newMap;
    end
     
     
    color_per_group = colormap(map);
    
    % plot each color in a seperate group (otherwise the legend won't work)
    for i=1:nGroups
        curr_group = groups(i);
%         if (curr_group~=0)
            curr_color_inds = (curr_group==gcolors);

            scatter(x(curr_color_inds),y(curr_color_inds),...
                    dot_size(curr_color_inds),...
                    color_per_group(i, :), 'fill');

             hold on; 
%        end
    end
    
    %Add contour to plot
    %plotContour(x, y);
    
	% find axis limits
    x_range  = max(x)-min(x);
    x_buffer = x_range*.03;
    x_lim    = [min(x)-x_buffer, max(x)+x_buffer];
    
    y_range  = max(y)-min(y);
    y_buffer = y_range*.03;
    y_lim    = [min(y)-y_buffer, max(y)+y_buffer];

    xlim(x_lim);
    ylim(y_lim);
end


%Add contour to plot
function plotContour(x1, x2)
    [~, density, x, y] = kde2d([x1, x2], 64);
    hold on;
%         cmap = jet;
%         cmap(1, :) = [1, 1, 1];
%         colormap(cmap);
    cmap=[0,0,0];
    colormap (cmap);
   
    contour(x, y, density, 8 );
    
    hold off;

end
