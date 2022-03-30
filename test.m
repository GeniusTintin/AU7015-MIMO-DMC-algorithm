%% step response of SISO system
sys = tf(1, [1,1]);
subplot(1,2,1)
step(sys);

sys2 = ss(-1,1,1,0);
subplot(1,2,2)
step(sys2)

%% step response of MIMO system
m = 5; k = 2; c = 0.1;
A=[0 1;-k/m, -c/m];
B=[1 1/m; 1/m 1];
C=[1 0; 0 1];
D=0;
sys_OL=ss(A,B,C,D);
step(sys_OL)

eig(A)

%% Stable MIMO 

J = [8 -3 -3; -3 8 -3; -3 -3 8];
F = 0.2*eye(3);
A = -J\F;
B = inv(J);
C = eye(3);
D = 0;
sys = ss(A,B,C,D);
step(sys)
eig(A)
%%

A = [-7,0;0,-10];
B = [5,0;0,2];
C = [1,-4;-4,0.5];
D = [0,-2;2,0];
% ts = 0.2;
sys = ss(A,B,C,D);
step(sys)

isstable(sys)

%% ode testing
theta = 1;
F  = @(t, x) -x + theta;
[t, x] = ode45(F, [0 10], 0);
plot(t, x)
x(end)