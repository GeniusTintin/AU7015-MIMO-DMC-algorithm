%% SISO DMC simulation 
% -- Author: Dingran Yuan
% -- Date: 29/03/2022

clearvars
clc;

%% Model for simulation
num = 1;
den = [5,1];
sys = tf(num, den);
step(sys, 0:1:100)
grid on
%% Initialise parameters
Ts = 1; % Sampling time
P = 10; % Prediction horizon
M = 5; % Control horizon
N = 50; % Model length
w = 2; % Final reference value
y0 = 0.1; % Initial output

% Model parameter derived from step response
[a_step, t_step] = step(sys, 0:Ts:N - 1);

A = zeros(P, M); 
S = zeros(N, N);
Q = 5 * eye(P);
R = 1 * eye(M);
H = ones(N, 1);
YN = zeros(N, 1);

y_out = 0; % output computed by Runge-Kutta method

% Construct A matrix
for i = 1: P
    for j = 1: M
        if(i-j+1 > 0)
            A(i, j) = a_step(i - j + 1);
        end
    end
end

% Construction for Shifting matrix S
S(1:end-1, 2:end) = eye(N-1, N-1);
S(N,N) = 1;

% Construction c
c = zeros(M, 1);
c(1) = 1;

% Calculate dT offline
dT = c' / (A'*Q*A + R) * A' * Q;

%% Start DMC
t_sim = 100.0; % Simulate for 20s
mpc_time = 0.0;

x_array = [y0]; % Store output for plotting

u = 0; % Control signal
Ycurr = y0;
figure
while(mpc_time < t_sim)
    
    wP = w * ones(P, 1);
    YP = YN(1:P); % Crop from YN
    du = dT * (wP - YP); % Calculate the control input
    
    %// control of the YN
    YN1 = YN + a_step * du;

    u = u + du;
    F = @(t, x) -5*x + u;
    [t, Ytmp] = ode45(F, [0, Ts], Ycurr);
    Ycurr = Ytmp(end); % Compute the model output

    e1 = Ycurr - YN1(1);

    Ycor = YN1 + H * e1;
    x_array = [x_array, Ycor(1)]; % save the current value to the array
    YN = S * Ycor; % Shift the predicted model
    
    mpc_time = mpc_time + Ts;

    %% visualise by step
    plot(0:Ts:mpc_time, x_array);
    grid on
    grid minor
    syms s 
    n = sym(num);
    d = sym(den);
    ns = poly2sym(n,s);
    ds = poly2sym(d,s);
    tfsym = ns/ds;
    tftitle = latex(tfsym);
    title(sprintf('SISO MPC for system: $$ %s $$', tftitle), 'Interpreter','latex')
    xlabel('{t} (sec)');
    ylabel('Output value');
    ylim([-3,4])
    hold on
    yline(w,'r--','Reference')
    hold on
    plot(mpc_time:Ts:mpc_time + (P-1)*Ts,YP, 'g.') % prediction
    hold off
    pause(0.1)
end

