% Declare
%% Quadrotor state variables(13*1)
syms x y z vx vy vz ...
     wx wy wz e1 e2 e3 e4
X = [x y z vx vy vz ...
     wx wy wz e1 e2 e3 e4];


%% other parameters
syms m g 
e3_vec = [0 0 1].';
syms Jx Jy Jz
J = [Jx 0 0;
     0 Jy 0;
     0  0 Jz];

syms r c 
allocation_matrix = [   1 1 1 1;
             -r r r -r;
             r -r r -r;
             c c -c -c];
% control_input
syms f1 f2 f3 f4 
F_vec = [f1 f2 f3 f4].';
E_diag = diag([e1 e2 e3 e4]);
control_input = allocation_matrix*E_diag*F_vec;
f = control_input(1);
M = control_input(2:4);

%rotation matrix
syms r1 r2 r3 r4 r5 r6 r7 r8 r9 
R = [r1 r2 r3;
     r4 r5 r6;
     r7 r8 r9];

w = [wx wy wz].';

%% dynamics
a = f*R*e3_vec/m - g*e3_vec;
dw = inv(J)*(M - cross(w,J*w));

%% measurement model (6*1) 
Y = [x y z wx wy wz].';

%% X_dot
% quadrotor state_dot(13*1)
X_dot = [ vx vy vz a(1) a(2) a(3) ...
          dw(1) dw(2) dw(3) ... 
          0 0 0 0].'; 

%% observability matrix
% observability matrix first round 
obs_row1 = jacobian(Y,X);
% observability matrix second round
Lie1 = obs_row1 * X_dot; 
obs_row2 = jacobian(Lie1,X);
% observability matrix third round
Lie2 = obs_row2 * X_dot;
obs_row3 = jacobian(Lie2,X);

observability_matrix = [obs_row1;obs_row2;obs_row3];
rank(observability_matrix)