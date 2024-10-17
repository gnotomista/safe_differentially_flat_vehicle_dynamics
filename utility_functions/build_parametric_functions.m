function [h, h_dot, h_dot_dot, h_p, h_p_dot, h_pp, zeta] = build_parametric_functions(o, D, GAMMA, GAMMA_P, GAMMA_PP, ALPHA, phi)

h = @(x1, x2, x3) norm(x1 - o) ^ 2 - D ^ 2;
h_dot = @(x1, x2, x3) 2 * (x1 - o)' * x2;
h_dot_dot = @(x1, x2, x3) 2 * (x1 - o)' * x3 + 2 * norm(x2) ^ 2;

h_p = @(x1, x2, x3) h_dot(x1, x2, x3) + GAMMA * h(x1, x2, x3);
h_p_dot = @(x1, x2, x3) h_dot_dot(x1, x2, x3) + GAMMA * h_dot(x1, x2, x3);

h_pp = @(x1, x2, x3) h_p_dot(x1, x2, x3) + GAMMA_P * h_p(x1, x2, x3);

zeta = @(x1, x2, x3, r) (...
    2 * (x1 - o)' * ALPHA * phi(x1, x2, x3, r) + 6 * x2' * x3 + ...
    2 * (GAMMA_P + GAMMA_PP) * (x1 - o)' * x3 + ...
    2 * (GAMMA_P + GAMMA_PP) * norm(x2) ^ 2 + ...
    2 * (GAMMA*GAMMA_P + GAMMA_P*GAMMA_PP) * (x1 - o)' * x2 + ...
    GAMMA*GAMMA_P*GAMMA_PP * (norm(x1 - o) ^ 2 - D ^ 2) ...
    ) / 2;

end