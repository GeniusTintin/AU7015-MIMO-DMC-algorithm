%% MIMO DMC simulation 
% -- Author: Dingran Yuan
% -- Date: 29/03/2022

clearvars
clc;

%% Model simulation (MIMO state space model)
A_sys = [-2,-1;6,-4];
B_sys = [4,-2;6,3];
C_sys = [1,0;0,1];
D_sys = [-1,-2;-1,1];

sys = ss(A_sys,B_sys,C_sys,D_sys);
step(sys, 0:0.1: 50*0.1)
grid on
grid minor

%% Initialise parameters
Ts = 0.1; % Sampling time
P = 20; % Prediction horizon
M = 5; % Control horizon
N = 50; % Model length

[a_step, t_step] = step(sys, 0:Ts:50 - Ts); % Step response of the system

nIn = size(a_step, 3); % number of inputs
nOut = size(a_step, 2); % number of outputs

a11 = a_step(1:50,1,1);
a12 = a_step(1:50,1,2); % Step response from in(2) to out(1)
a21 = a_step(1:50,2,1); % Step response from in(1) to out(2)
a22 = a_step(1:50,2,2);
% plot(t_step(1:50),a_step(1:50,2,1))

A11 = zeros(P, M); 
A12 = zeros(P, M); 
A21 = zeros(P, M); 
A22 = zeros(P, M); 

% coloum matrix for YN & YP 
YN = zeros(nOut*N, 1);
YP = zeros(nOut*P, 1);

% Elements used to compute D matrix
L = zeros(nIn, nIn*M);
S = zeros(nOut*N, nOut*N);
Q = zeros(nOut*P, nOut*P);
R = zeros(nIn*M, nIn*M);
H = zeros(nOut*N, nOut);

% Construct A matrices
for i = 1: P
    for j = 1: M
        if(i-j+1 > 0)
            A11(i, j) = a11(i - j + 1);
            A12(i, j) = a12(i - j + 1);
            A21(i, j) = a21(i - j + 1);
            A22(i, j) = a22(i - j + 1);
        end
    end
end
A = [A11, A12; A21, A22];  
AN = [a11, a12; a21, a22];

% weight for system error and input (here we penalise the error more)
Q_arg = [3 4]; % weight for system error
R_arg = [5 5]; % weight for control input
Out_ref = [5 10]; % reference of the output
wr = zeros(nOut * P, 1);

% Construct wr & Q & S & R & L
for i = 1:nOut
    % Construct wr
    wr((i-1)*P+1: i*P) = Out_ref(i) * ones(P, 1);

    % Construct Q
    Q((i-1)*P+1:i*P, (i-1)*P+1:i*P) = Q_arg(i)*eye(P);
    
    % Construct S
    S((i-1)*N+1:i*N-1, (i-1)*N+2:i*N) = eye(N-1);
    S(i*N, i*N) = 1;

    % Construct H
    H((i-1)*N + 1:i*N, i) = ones(N, 1);
end
for i = 1:nIn
    % Construct R
    R((i-1)*M+1:i*M, (i-1)*M+1:i*M) = R_arg(i)*eye(M);

    % Construct L
    L(i,(i-1)*M + 1) = 1;
end

% Offline computation of matrix D
D = L / (A'*Q*A+R)*A'*Q; 

%% MIMO DMC loop

t_sim = 15.0; % Simulation time
mpc_time = 0.0;
u = [0;0]; % Control signal
X_curr = [0;0]; % initial state

Y_array = C_sys*X_curr + D_sys*u; % Store output for plotting
Y1_min = Y_array(1); Y1_max = Y_array(1); 
Y2_min = Y_array(2); Y2_max = Y_array(2); 
u_array = [0;0];
figure

while(mpc_time < t_sim)
    for i = 1:nOut
        YP((i-1)*P+1:i*P) = YN((i-1)*N+1:(i-1)*N+P);
    end
    du = D * (wr - YP);

    YN1 = YN + AN * du;
    Y1 = [YN1(1);YN1(N+1)]; % The current estimation after control input
    
    Y1_min = min(Y1(1), Y1_min); Y1_max = max(Y1(1), Y1_max);
    Y2_min = min(Y1(2), Y2_min); Y2_max = max(Y1(2), Y2_max);

    u = u + du;
    [X_curr, Y_curr] = RKsolver(sys, u, X_curr, Ts);
    e1 = X_curr - Y1; % calculate the error with the real output
    Y_cor = YN1 + H*e1;
    YN = S * Y_cor;

    Y_array = [Y_array, [YN(1);YN(N+1)]];
    u_array = [u_array, u];
    mpc_time = mpc_time + Ts;

    %% visualise by step
    subplot(2,2,2)
    plot(linspace(0, mpc_time, size(Y_array,2)), Y_array(1,:), 'b-');
    ylim([Y1_min - 0.4*(Y1_max - Y1_min), Y1_max + 0.4*(Y1_max - Y1_min)])
    grid on
    grid minor
    yline(Out_ref(1), 'r--', 'Reference')
    xlabel('{t} (sec)');
    ylabel('Output value {y_1}');
    title('Output {y_1} overtime')

    subplot(2,2,1)
    plot(linspace(0, mpc_time, size(Y_array,2)), u_array(1,:),"Color",'#D95319');
    grid on
    grid minor
    xlabel('{t} (sec)');
    ylabel('Input {u_1}');
    title('Input {u_1} overtime')

    subplot(2,2,4)
    plot(linspace(0, mpc_time, size(Y_array,2)), Y_array(2,:), 'b-');
    ylim([Y2_min - 0.4*(Y2_max - Y2_min), Y2_max + 0.4*(Y2_max - Y2_min)])
    grid on
    grid minor
    xlabel('{t} (sec)');
    ylabel('Output value {y_2}');
    yline(Out_ref(2),'r--','Reference')
    title('Output {y_2} overtime')

    subplot(2,2,3)
    plot(linspace(0, mpc_time, size(Y_array,2)), u_array(2,:), "Color",'#D95319');
    grid on
    grid minor
    xlabel('{t} (sec)');
    ylabel('Input {u_2}');
    title('Input {u_2} overtime')

    drawnow
end

%% Runge-Kutta method for simulation
function [x1, y1] = RKsolver(sys, u, x0, Ts)
%--
%-> u: the input of the system
%-> sys: the state space model of the system
%<- Fu: the output differential equition 
    A = sys.A; B = sys.B;
    C = sys.C; D = sys.D;
    F = @(t, x) A*x + B*u;
    [~, X_tmp] = ode45(F, [0, Ts], x0);
%     plot(t, X_tmp)
    x1 = X_tmp(end,:)';
    y1 = C*x1 + D*u;
end

%% Appendix codes
% Construct the output reference matrix
% for i = 1:nOut
%     wr((i-1)*P+1: i*P) = Out_ref(i) * ones(P, 1);
% end

% Construct L matrix
% for i = 1:nIn
%     L(i,(i-1)*M + 1) = 1;
% end

% Construct S matrix
% for i = 1:nOut
%     S((i-1)*N+1:i*N-1, (i-1)*N+2:i*N) = eye(N-1);
%     S(i*N, i*N) = 1;
% end

% Construct H matrix
% for i = 1:nOut
%     H((i-1)*N + 1:i*N, i) = ones(N, 1);
% end