%% MIMO DMC simulation 
% -- Author: Dingran Yuan
% -- Date: 29/03/2022

clearvars
clc;

%% Model simulation (MIMO state space model)
A = [-2,1;1,-1];
B = [5,0;0,2];
C = [1,0;0,1];
D = [1,1;2,1];

sys = ss(A,B,C,D);
step(sys, 0:1:50);

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
Q_arg = [1 5]; % weight for system error
R_arg = [3 4]; % weight for control input
Out_ref = [20 5]; % reference of the output
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

t_sim = 100.0; % Simulate for 20s
mpc_time = 0.0;
u = [0;0]; % Control signal
X_curr = [0;0]; % initial state

Y_array = X_curr; % Store output for plotting
figure

while(mpc_time < t_sim)
    for i = 1:nOut
        YP((i-1)*P+1:i*P) = YN((i-1)*N+1:(i-1)*N+P);
    end
    du = D * (wr - YP);

    YN1 = YN + AN * du;
    Y1 = [YN1(1);YN1(N+1)]; % The current estimation after control input
    
    u = u + du;
    [X_curr, Y_curr] = RKsolver(sys, u, X_curr, Ts);
    e1 = X_curr - Y1; % calculate the error with the real output
    Y_cor = YN1 + H*e1;
    YN = S * Y_cor;
    Y_array = [Y_array, [YN(1);YN(N+1)]];
    mpc_time = mpc_time + Ts;

    %% visualise by step
    subplot(1,2,1)
    plot(linspace(0, mpc_time, size(Y_array,2)), Y_array(1,:), 'b-');
    grid on
    grid minor
    yline(Out_ref(1), 'r--', 'Reference')
    xlabel('{t} (sec)');
    ylabel('Output value');
    drawnow
    subplot(1,2,2)
    plot(linspace(0, mpc_time, size(Y_array,2)), Y_array(2,:), 'b-');
    grid on
    grid minor
    xlabel('{t} (sec)');
    ylabel('Output value');
    yline(Out_ref(2),'r--','Reference')
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
    [t, X_tmp] = ode45(F, [0, Ts], x0);
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