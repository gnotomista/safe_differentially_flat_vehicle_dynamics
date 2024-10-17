function h_p = plot_bicycle_rearaxle(x, l_f, l_r, h_p, varargin)
l = l_f + l_r;
NO_FIGURE = false;
if ~isempty(varargin)
    if strcmp(varargin{1}, 'dynamic')
        x_tmp = zeros(5, 1);
        x_tmp(1:2) = x(1:2);
        x_tmp(3) = x(5);
        x_tmp(5) = x(7);
        x = x_tmp;
    end
    if any(cellfun(@(x) strcmp(x, 'no-figure'), varargin))
        NO_FIGURE = true;
    end
end
if isempty(h_p)
    if ~NO_FIGURE
        figure('units', 'normalized', 'position', [0 0 1 1])
        hold on
        grid on
        axis equal
        axis([-10, 10, -10, 10])
    end
    h_p{1} = line(x(1) + cos(x(3)) * [0, l], x(2) + sin(x(3)) * [0, l], 'Color', [0.5 0.5 0.5], 'LineWidth', 4);
    h_p{2} = line(x(1) + l * cos(x(3)) + l / 4 * cos(x(3) + x(5)) * [-1, 1], ...
        x(2) + l * sin(x(3)) + l / 4 * sin(x(3) + x(5)) * [-1, 1], 'Color', 'k', 'LineWidth', 4);
else
    h_p{1}.XData = x(1) + cos(x(3)) * [0, l];
    h_p{1}.YData = x(2) + sin(x(3)) * [0, l];
    h_p{2}.XData = x(1) + l * cos(x(3)) + l / 4 * cos(x(3) + x(5)) * [-1, 1];
    h_p{2}.YData = x(2) + l * sin(x(3)) + l / 4 * sin(x(3) + x(5)) * [-1, 1];
end
end