function [fig, ax] = EulerRiemannExactPlot( ...
    rho_l, u_l, p_l, rho_r, u_r, p_r, x_l, x_r, t, x_c )

    if nargin < 10
        x_c = 0.0;
    end

    xlist = linspace(x_l, x_r, 1000);

    gamma = 1.4;
    [rho, u, p, more_info] = EulerRiemannExactSolver( ...
        rho_l, u_l, p_l, rho_r, u_r, p_r, gamma, xlist, x_c, t);

    % Calculate internal energy
    e = zeros(size(p));
    for w = 1:length(e)
        e(w) = p(w) / ((gamma - 1) * rho(w));
    end

    primitive = {rho, u, p, e};
    names = {"Density", "Velocity", "Pressure", "Internal Energy"};

    fig = figure;
    for w = 1:4
        subplot(2, 2, w);
        q = primitive{w};
        plot(xlist, q, 'LineWidth', 2);
        title(names{w});

        qmax = max(q);
        qmin = min(q);
        qdiff = qmax - qmin;
        ylim([qmin - 0.1 * qdiff, qmax + 0.1 * qdiff]);
    end

    % Plot the wave lines
    figure;
    ax = axes;
    hold on;
    tlist = linspace(0, t, 500);

    % Plot left shock or rarefaction
    if strcmp(more_info.left_type, 'Shock')
        w_1_line = x_c + more_info.w_1_l * tlist;
        plot(w_1_line, tlist, 'r', 'DisplayName', 'Left shock');
    else
        w_1_l_line = x_c + more_info.w_1_l * tlist;
        w_1_r_line = x_c + more_info.w_1_r * tlist;
        fill_betweenx(tlist, w_1_l_line, w_1_r_line, 'r', 0.2, 'Left Rarefaction');
    end

    % Contact discontinuity in the middle
    w_2_line = x_c + more_info.w_2 * tlist;
    plot(w_2_line, tlist, 'g', 'DisplayName', 'Contract discontinuity');

    % Plot right shock or rarefaction
    if strcmp(more_info.right_type, 'Shock')
        w_3_line = x_c + more_info.w_3_l * tlist;
        plot(w_3_line, tlist, 'b', 'DisplayName', 'Right shock');
    else
        w_3_l_line = x_c + more_info.w_3_l * tlist;
        w_3_r_line = x_c + more_info.w_3_r * tlist;
        fill_betweenx(tlist, w_3_l_line, w_3_r_line, 'b', 0.2, 'Right Rarefaction');
    end

    % Adjust axes and labels
    ax.XAxisLocation = 'bottom';
    ax.YAxisLocation = 'left';
    ax.XTick = [x_l, x_c, x_r];
    ax.YTick = t;
    ax.YTickLabel = {'t'};
    xlim([x_l, x_r]);
    ylim([0, t]);
    legend;
    hold off;

    subtitle(sprintf('Riemann Problem (t = %.4f)', t));
end

function fill_betweenx(tlist, y1, y2, color, alpha, label)
    % Fill between two curves, adjusted for x and y
    patch([y1, fliplr(y2)], [tlist, fliplr(tlist)], color, ...
        'FaceAlpha', alpha, 'EdgeColor', 'none', 'DisplayName', label);
end
