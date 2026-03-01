function violin(data, varargin)
    % violin creates a half-violin plot with a "rain" scatter plot
    %
    % Usage:
    %   half_violin_with_rain(data)
    %   half_violin_with_rain(data, 'param', value, ...)
    %
    % Inputs:
    %   data - a vector of data points
    %
    % Optional Parameters:
    %   'x'        - x location of the violin (default: 1)
    %   'color'    - color of the violin and points (default: [0.5 0.5 0.5])
    %   'width'    - width of the violin (default: 0.3)
    %   'side'     - side of the violin to plot ('left' or 'right', default: 'right')
    %   'pointSize'- size of the points in the rain plot (default: 10)
    %
    % Example:
    %   data = randn(100,1);
    %   violin(data, 'x', 1, 'color', [0.1 0.2 0.5], 'width', 0.5, 'side', 'right');

        % Default parameters
    p = inputParser;
    addParameter(p, 'x', 1);
    addParameter(p, 'color', [0.5 0.5 0.5]);
    addParameter(p, 'width', 0.3);
    addParameter(p, 'side', 'right');
    addParameter(p, 'pointSize', 10);
    addParameter(p, 'jitter', 0.2);
    parse(p, varargin{:});

    x = p.Results.x;
    color = p.Results.color;
    width = p.Results.width;
    side = p.Results.side;
    pointSize = p.Results.pointSize;
    jitter = p.Results.jitter;

    % Calculate the kernel density estimate
    [f, xi] = ksdensity(data);
    f = f / max(f) * width; % Normalize and scale by width

    % Create the half-violin plot
    hold on;
    if strcmp(side, 'right')
        fill([x + f, x * ones(size(f))], [xi, fliplr(xi)], color, 'EdgeColor', 'none');
        %plot(x + f, xi, 'k-', 'LineWidth', 1);
        scatter(x - width/2 + (rand(size(data)) - 0.5) * jitter, data, pointSize, 'filled', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none');
    elseif strcmp(side, 'left')
        fill([x - f, x * ones(size(f))], [xi, fliplr(xi)], color, 'EdgeColor', 'none');
        %plot(x - f, xi, 'k-', 'LineWidth', 1);
        scatter(x + width/2 + (rand(size(data)) - 0.5) * jitter, data, pointSize, 'filled', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none');
    end
    %plot([x x], [min(data), max(data)], 'k-', 'LineWidth', 1);
    hold off;
end