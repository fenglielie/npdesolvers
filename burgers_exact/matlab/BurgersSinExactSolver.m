% Exact solution of Burgers' equation with 2pi periodic boundary condition
% u_t + u*u_x = 0,
% u0(x) = alpha + beta sin(x) in [-pi, pi]
function u = BurgersSinExactSolver(x,t,alpha,beta)
    %  u(x,0) = alpha + beta*sin(x)
    u = alpha + beta * BurgersSinNewton(x - alpha*t, beta*t);
end


function u = BurgersSinNewton(x, t)
    % u(x,0) = sin(x)
    x = mod(x + pi, 2*pi) - pi; % x in [-pi,pi]

    iter = 1;
    u = x/(pi/2+t); % u0
    while iter < 1e5
        du = (u-sin(x-u*t))./(1 + cos(x - u*t)*t);
        u = u - du;
        if max(abs(du)) < 1e-10
            break
        end
        iter = iter + 1;
    end
end
