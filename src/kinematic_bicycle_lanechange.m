clc
clear
close all

WITH_CBFS = false;

% simulation and control parameters
DT = 0.01;
T = 0.6;
l_f = 1.4;
l_r = 1.6;
SIGMA = 2.0;
L = 4;

[f, g, full_dyn] = get_kinematic_bicycle_rearaxle(l_f, l_r);

% initialize system state
x0 = [SIGMA; 0; pi/2; 10; 0];
t0 = 0;
x = x0;
% initial state of chain of integrators (dynamics of flat outputs)
x_df = [x0(1:2);
    x0(4) * cos(x0(3));
    x0(4) * sin(x0(3))];
u_df = zeros(2,1);

% desired final pose
xf = [-SIGMA; L; pi/2; 10; 0];
tf = 1;

% plan trajectory
t_straight = 6;
v_straight = 5;
y_offset = 5;
A = [1 t0 t0^2 t0^3 0 0 0 0;
    0 0 0 0 1 t0 t0^2 t0^3;
    1 tf tf^2 tf^3 0 0 0 0;
    0 0 0 0 1 tf tf^2 tf^3;
    0 1 2*t0 3*t0^2 0 0 0 0;
    0 0 0 0 0 1 2*t0 3*t0^2;
    0 1 2*tf 3*tf^2 0 0 0 0;
    0 0 0 0 0 1 2*tf 3*tf^2];
c = [x0(1); x0(2); xf(1); xf(2); x0(4)*cos(x0(3)); x0(4)*sin(x0(3)); xf(4)*cos(xf(3)); xf(4)*sin(xf(3))];
param = A \ c;
a0 = param(1);
a1 = param(2);
a2 = param(3);
a3 = param(4);
b0 = param(5);
b1 = param(6);
b2 = param(7);
b3 = param(8);

% tracking and safety parameters
ALPHA = 5;
x1 = x_df(1:2);
x2 = x_df(3:4);
if x1(1) < -SIGMA
    error('unsafe start')
end
if x2(1) < 0
    GAMMA = -1 * (x1(1) + SIGMA) / x2(1);
else
    GAMMA = 1e1;
end
if u_df(1) + GAMMA * x2(1) < 0
    GAMMA_P = -1 * (x2(1) + GAMMA * (x1(1) + SIGMA)) / (u_df(1) + GAMMA * x2(1));
else
    GAMMA_P = 3e0;
end
GAMMA_PP = 1e2;

% log
logger.r = [];
logger.x = [];
logger.u = [];

% plot
h_p = plot_bicycle_rearaxle(x, l_f, l_r, []);
handle_x = scatter(x1' * [1; 0], x1' * [0; 1], 500, 'r.');
handle_x_traj = plot(x1' * [1; 0], x1' * [0; 1], 'r', 'LineWidth', 2);
handle_r = scatter(nan, nan, 100, 'bo');
handle_r_traj = plot(nan, nan, 'b', 'LineWidth', 2);
axis([-5*SIGMA, 5*SIGMA, 0, 100])

% main simulation loop
for t = t0 : DT : 2 * t_straight + tf

    % compute reference
    if t < t_straight
        y1 = SIGMA;
        y2 = y_offset + v_straight*t;
        y1_dot = 0;
        y2_dot = v_straight;
        y1_dot_dot = 0;
        y2_dot_dot = 0;
        y1_dot_dot_dot = 0;
        y2_dot_dot_dot = 0;
    elseif t < t_straight + tf
        y1 = a0 + a1*(t-t_straight) + a2*(t-t_straight)^2 + a3*(t-t_straight)^3;
        y2 = y_offset + v_straight*t_straight + b0 + b1*(t-t_straight) + b2*(t-t_straight)^2 + b3*(t-t_straight)^3;
        y1_dot = a1 + 2*a2*(t-t_straight) + 3*a3*(t-t_straight)^2;
        y2_dot = b1 + 2*b2*(t-t_straight) + 3*b3*(t-t_straight)^2;
        y1_dot_dot = 2*a2 + 6*a3*(t-t_straight);
        y2_dot_dot = 2*b2 + 6*b3*(t-t_straight);
        y1_dot_dot_dot = 6*a3;
        y2_dot_dot_dot = 6*b3;
    else
        y1 = -SIGMA;
        y2 = y_offset + v_straight*t_straight ...
            + (b0 + b1*tf + b2*tf^2 + b3*tf^3) ...
            + b1 * (t-t_straight-tf);
        y1_dot = 0;
        y2_dot = b1;
        y1_dot_dot = 0;
        y2_dot_dot = 0;
        y1_dot_dot_dot = 0;
        y2_dot_dot_dot = 0;
    end

    % get states
    x1 = x(1:2);
    x2 = x(4) * [cos(x(3)); sin(x(3))];

    % compute tracking controller
    phi = 2 / T ^ 2 * ([y1;y2] - x1) - 2 / T * x2 - u_df;

    % compute safety controller (closed form)
    zeta = ALPHA * phi(1) ...
        + (GAMMA + GAMMA_P + GAMMA_PP) * u_df(1) ...
        + (GAMMA*GAMMA_P + GAMMA*GAMMA_PP + GAMMA_P*GAMMA_PP) * x2(1) ...
        + GAMMA*GAMMA_P*GAMMA_PP * (x1(1) + SIGMA);
    w1 = max(0, -zeta);
    w = [w1; 0];

    % simulate flat output dynamics
    v_df = ALPHA * phi ...
        + WITH_CBFS * w;
    x_dot_df = [zeros(2), eye(2);
        zeros(2), zeros(2)] * x_df + [zeros(2); eye(2)] * u_df;
    u_dot_df = v_df;
    x_df = x_df + x_dot_df * DT;
    u_df = u_df + u_dot_df * DT;

    % compute system states and inputs from flat outputs
    x1 = x_df(1:2);
    x2 = x_df(3:4);
    if norm(x2) < eps
        a = 0;
        delta_dot = 0;
    else
        a = (x_df(3)*u_df(1) + x_df(4)*u_df(2)) / sqrt(x_df(3)^2+x_df(4)^2);
        delta_dot = (l_f + l_r) / norm(x2)^6 * ((u_dot_df(2) * x2(1) - u_dot_df(1) * x2(2))*norm(x2)^3 - (u_df(2) * x2(1) - u_df(1) * x2(2))*3*norm(x2)*(x2'*u_df));
    end

    % simulate system dynamics
    u = [a; delta_dot];
    x_dot = f(x) + g(x) * u;
    x = x + x_dot * DT;

    % log
    logger.r(:, end+1) = [y1;y2];
    logger.x(:, end+1) = x;
    logger.u(:, end+1) = u;

    % plot
    h_p = plot_bicycle_rearaxle(x, l_f, l_r, h_p);
    handle_x.XData = x(1);
    handle_x.YData = x(2);
    handle_x_traj.XData = [handle_x_traj.XData, x(1)];
    handle_x_traj.YData = [handle_x_traj.YData, x(2)];
    handle_r.XData = y1;
    handle_r.YData = y2;
    handle_r_traj.XData = [handle_r_traj.XData, y1];
    handle_r_traj.YData = [handle_r_traj.YData, y2];
    drawnow limitrate

end
