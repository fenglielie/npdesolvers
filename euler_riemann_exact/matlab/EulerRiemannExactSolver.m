function [rho_out, u_out, p_out, more_info] = EulerRiemannExactSolver( ...
        rho_l, u_l, p_l, rho_r, u_r, p_r, gamma, xlist, x_c, t)

    if nargin < 7
        gamma = 1.4;
    end
    if nargin < 8
        xlist = linspace(-1.0, 1.0, 1000);
    end
    if nargin < 9
        x_c = 0.0;
    end
    if nargin < 10
        t = 0.8 * max(abs(xlist - x_c)) / max(abs([u_l, u_r]));
    end

    % Calculate sound speed in left and right states
    c_l = sqrt(gamma * p_l / rho_l);
    c_r = sqrt(gamma * p_r / rho_r);

    % Calculate constants alpha and beta
    alpha = (gamma + 1.0) / (gamma - 1.0);
    beta = (gamma - 1.0) / (2.0 * gamma);

    % Check for cavitation (solution for cavitation not supported)
    if u_l - u_r + 2 * (c_l + c_r) / (gamma - 1.0) < 0
        disp('Cavitation detected! Exiting.');
        rho_out = []; u_out = []; p_out = []; more_info = [];
        return;
    end

    % Define integral curves and Hugoniot locus for 1-wave and 3-wave

    integral_curve_1 = @(p) u_l + 2 * c_l / (gamma - 1.0) * (1.0 - (p / p_l) ^ beta);
    integral_curve_3 = @(p) u_r - 2 * c_r / (gamma - 1.0) * (1.0 - (p / p_r) ^ beta);

    hugoniot_locus_1 = @(p) u_l + 2 * c_l / sqrt(2 * gamma * (gamma - 1.0)) * ...
        ((1 - p / p_l) / sqrt(1 + alpha * p / p_l));
    hugoniot_locus_3 = @(p) u_r - 2 * c_r / sqrt(2 * gamma * (gamma - 1.0)) * ...
        ((1 - p / p_r) / sqrt(1 + alpha * p / p_r));

    phi_l = @(p) (p >= p_l) * hugoniot_locus_1(p) + (p < p_l) * integral_curve_1(p);
    phi_r = @(p) (p >= p_r) * hugoniot_locus_3(p) + (p < p_r) * integral_curve_3(p);

    % Construct the intersection equation in the (p-v) plane
    func = @(p) phi_l(p) - phi_r(p);

    % Initial guess p0
    p0_PV = (p_l + p_r) / 2.0 - 1 / 8 * (u_r - u_l) * (rho_l + rho_r) * (c_l + c_r);
    p0 = max(p0_PV, 1e-8);

    % Solve for the intersection point to get the intermediate state p and u
    options = optimset('TolFun', 1.0e-12,'Display', 'off');
    [p_s_tmp, ~, exitflag] = fsolve(func, p0, options);
    p_s = p_s_tmp;
    u_s = 0.5 * (phi_l(p_s) + phi_r(p_s));

    % Warning if the solution fails to converge
    if exitflag <= 0
        disp('Warning: fsolve did not converge.');
    end

    % Calculate the density on the left and right side of the contact discontinuity

    if p_s <= p_l
        rho_s_l = (p_s / p_l) ^ (1.0 / gamma) * rho_l;  % Rarefaction wave
    else
        rho_s_l = ((1.0 + alpha * p_s / p_l) / ((p_s / p_l) + alpha)) * rho_l;  % Shock wave
    end

    if p_s <= p_r
        rho_s_r = (p_s / p_r) ^ (1.0 / gamma) * rho_r;  % Rarefaction wave
    else
        rho_s_r = ((1.0 + alpha * p_s / p_r) / ((p_s / p_r) + alpha)) * rho_r;  % Shock wave
    end

    % Calculate sound speed in the intermediate state
    c_s_l = sqrt(gamma * p_s / rho_s_l);
    c_s_r = sqrt(gamma * p_s / rho_s_r);

    % Calculate wave speeds

    % 1-wave
    if p_s > p_l  % Shock wave
        w_1_l = (rho_l * u_l - rho_s_l * u_s) / (rho_l - rho_s_l);
        w_1_r = w_1_l;
    else  % Rarefaction wave
        w_1_l = u_l - c_l;
        w_1_r = u_s - c_s_l;
    end

    % 2-wave
    w_2 = u_s;

    % 3-wave
    if p_s > p_r  % Shock wave
        w_3_l = (rho_r * u_r - rho_s_r * u_s) / (rho_r - rho_s_r);
        w_3_r = w_3_l;
    else  % Rarefaction wave
        w_3_l = u_s + c_s_r;
        w_3_r = u_r + c_r;
    end

    w_max = max(abs([w_1_l, w_1_r, w_2, w_3_l, w_3_r]));

    % Warning if the wave is out of range
    if t * w_max > max(abs(xlist - x_c))
        disp('Warning: wave is out of range.');
    end

    % Solve for the state inside the rarefaction wave
    xi = (xlist - x_c) / t;
    u_1_fan = ((gamma - 1.0) * u_l + 2 * (c_l + xi)) / (gamma + 1.0);
    u_3_fan = ((gamma - 1.0) * u_r - 2 * (c_r - xi)) / (gamma + 1.0);
    rho_1_fan = (rho_l ^ gamma * (u_1_fan - xi) .^ 2 / (gamma * p_l)) .^ (1.0 / (gamma - 1.0));
    rho_3_fan = (rho_r ^ gamma * (xi - u_3_fan) .^ 2 / (gamma * p_r)) .^ (1.0 / (gamma - 1.0));
    p_1_fan = p_l * (rho_1_fan / rho_l) .^ gamma;
    p_3_fan = p_r * (rho_3_fan / rho_r) .^ gamma;

    % Calculate return values

    rho_out = zeros(size(xlist));
    u_out = zeros(size(xlist));
    p_out = zeros(size(xlist));

    for i = 1:length(xi)
        if xi(i) <= w_1_l  % Left of the 1-wave
            rho_out(i) = rho_l;
            u_out(i) = u_l;
            p_out(i) = p_l;
        elseif xi(i) <= w_1_r  % Inside the 1-wave (if it's a rarefaction wave)
            rho_out(i) = rho_1_fan(i);
            u_out(i) = u_1_fan(i);
            p_out(i) = p_1_fan(i);
        elseif xi(i) <= w_2  % Between the 1-wave and the 2-wave
            rho_out(i) = rho_s_l;
            u_out(i) = u_s;
            p_out(i) = p_s;
        elseif xi(i) <= w_3_l  % Between the 2-wave and the 3-wave
            rho_out(i) = rho_s_r;
            u_out(i) = u_s;
            p_out(i) = p_s;
        elseif xi(i) <= w_3_r  % Inside the 3-wave (if it's a rarefaction wave)
            rho_out(i) = rho_3_fan(i);
            u_out(i) = u_3_fan(i);
            p_out(i) = p_3_fan(i);
        else  % Right of the 3-wave
            rho_out(i) = rho_r;
            u_out(i) = u_r;
            p_out(i) = p_r;
        end
    end

    % Additional info
    more_info = struct();
    more_info.p_s = p_s;
    more_info.u_s = u_s;
    more_info.rho_s_l = rho_s_l;
    more_info.rho_s_r = rho_s_r;
    more_info.w_1_l = w_1_l;
    more_info.w_1_r = w_1_r;
    more_info.w_2 = w_2;
    more_info.w_3_l = w_3_l;
    more_info.w_3_r = w_3_r;

    if p_s > p_l
        left_type = 'Shock';
    else
        left_type = 'Rarefaction';
    end

    center_type = 'Contract discontinuity';

    if p_s > p_r
        right_type = 'Shock';
    else
        right_type = 'Rarefaction';
    end

    more_info.left_type = left_type;
    more_info.center_type = center_type;
    more_info.right_type = right_type;

    more_info.type_msg = sprintf('[L] %s\n[C] %s\n[R] %s', left_type, center_type, right_type);
    more_info.key_msg = sprintf('%s=%f\t%s=%f\t%s=%f\t%s=%f', 'p_s', p_s, 'u_s', u_s, 'rho_s_l', rho_s_l, 'rho_s_r', rho_s_r);
end
