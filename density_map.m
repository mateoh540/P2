function d = density_map(x,y,t)
%{
Let "n" be the length of the grid.

--- Inputs ---
x, y: n by n meshgrid of x and y coordinates
   t: time  

--- Outputs ---
d: n by n meshgrid of densities.  All densities must be on [0,inf)
%}

    % Calculate the grid size
    n = numel(x);

    % Initialize density map
    d = zeros(size(x));

    % Define the forest fire points, this will govern the strength of the
    % fire at that point
    fire_points = [x(:), y(:)]; 

    % Compute distances and update density map
    for i = 1:n
        for j = 1:n
            % Compute distances from each point on the grid to fire_points
            distances = sqrt((fire_points(:,1) - x(i, j)).^2 + (fire_points(:,2) - y(i, j)).^2);

            % Compute density as inversely proportional to distance
            d(i, j) = sum(1 ./ distances); % Aggregate inverse distances
        end
    end

    % Ensure densities are on [0, inf)
    d = max(d, 0); % Set negative densities to 0

    

end
