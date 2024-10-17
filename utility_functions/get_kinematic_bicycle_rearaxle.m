function [f, g, dyn] = get_kinematic_bicycle_rearaxle(l_f, l_r)
    function f_v = f_h(x)
        f_v = [x(4) * cos(x(3));
            x(4) * sin(x(3));
            x(4) / (l_f + l_r) * tan(x(5));
            0;
            0];
    end
    function g_v = g_h(x)
        g_v = [zeros(3,2); eye(2)];
    end
    function x_dot = kinematic_bicycle_h(x, u)
        x_dot = f(x) + g(x) * u;
    end
f = @f_h;
g = @g_h;
dyn = @kinematic_bicycle_h;
end