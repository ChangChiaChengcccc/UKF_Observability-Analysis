% Declare
%% Quadrotor1 state variables(13*1)
syms x y z vx vy vz ...
     wx wy wz e1 e2 e3 e4
X = [x y z vx vy vz ...
     wx wy wz e1 e2 e3 e4];

%% Quadrotor2 state variables(13*1)
syms x2 y2 z2 vx2 vy2 vz2 ...
     wx2 wy2 wz2 e12 e22 e32 e42
X2 = [x2 y2 z2 vx2 vy2 vz2 ...
     wx2 wy2 wz2 e12 e22 e32 e42];

%% Payload state variables(15*1)
syms xp yp zp vxp vyp vzp ...
     wxp wyp wzp fx fy fz fx2 fy2 fz2
Xp = [xp yp zp vxp vyp vzp ...
     wxp wyp wzp fx fy fz fx2 fy2 fz2];

%% X overall state variables(41*1)
X = [X X2 Xp];

%% other parameters
syms m m2 mp g 
e3_vec = [0 0 1].';
syms Jx Jy Jz Jx2 Jy2 Jz2 Jxp Jyp Jzp
J = [Jx 0 0;
     0 Jy 0;
     0  0 Jz];
J2 = [Jx2 0 0;
     0 Jy2 0;
     0  0 Jz2];
Jp = [Jxp 0 0;
     0 Jyp 0;
     0  0 Jzp];

syms r c 
allo_mat = [   1 1 1 1;
             -r r r -r;
             r -r r -r;
             c c -c -c];
% control_input
syms f1 f2 f3 f4 
F_vec = [f1 f2 f3 f4].';
E_diag = diag([e1 e2 e3 e4]);
control_input = allocation_matrix*E_diag*F_vec;
f_c = control_input(1);
M = control_input(2:4);

% control_input2
syms f12 f22 f32 f42 
F_vec2 = [f12 f22 f32 f42].';
E_diag2 = diag([e12 e22 e32 e42]);
control_input2 = allocation_matrix*E_diag2*F_vec2;
f2_c = control_input2(1);
M2 = control_input2(2:4);

%rotation matrix
syms r1 r2 r3 r4 r5 r6 r7 r8 r9 
R = [r1 r2 r3;
     r4 r5 r6;
     r7 r8 r9];
%% pre dynamics
w = [wx wy wz].';
F = [fx fy fz].';
F2 = [fx2 fy2 fz2].';
syms tx ty tz tx2 ty2 tz2 
Tau = [tx ty tz].';
Tau2 = [tx2 ty2 tz2].';
F_mat = [0 0 0;
         0 0 0.5;
         0 -0.5 0];
F2_mat = [0 0 0;
         0 0 -0.5;
         0 0.5 0];
Mp = F_mat*F + F2_mat*F2;

%% dynamics
a = f_c*R*e3_vec/m - g*e3_vec + F; 
dw = inv(J)*(M - cross(w,J*w) + Tau);

a2 = f2_c*R*e3_vec/m2 - g*e3_vec + F2; 
dw2 = inv(J2)*(M2 - cross(w,J2*w) + Tau2); 

ap = -(F+F2)/mp - g*e3_vec;
dwp = inv(Jp)*(Mp - cross(w,Jp*w) -(Tau+Tau2) ); 


%% measurement model (24*1) 
Y = [x y z wx wy wz tx ty tz x2 y2 z2 wx2 wy2 wz2 tx2 ty2 tz2 xp yp zp wxp wyp wzp].';

%% X_dot
% overall state_dot(*41)
X_dot = [ vx vy vz a(1) a(2) a(3) ...
          dw(1) dw(2) dw(3) ... 
          0 0 0 0 ...
          vx2 vy2 vz2 a2(1) a2(2) a2(3) ...
          dw2(1) dw2(2) dw2(3) ... 
          0 0 0 0 ...
          vxp vyp vzp ap(1) ap(2) ap(3) ...
          dwp(1) dwp(2) dwp(3) ... 
          0 0 0 ... %F
          0 0 0].'; %F2

%% observability matrix
% observability matrix first round 
obs_row1 = jacobian(Y,X);
% observability matrix second round
Lie1 = obs_row1 * X_dot; 
obs_row2 = jacobian(Lie1,X);
% observability matrix third round
Lie2 = obs_row2 * X_dot;
obs_row3 = jacobian(Lie2,X);
% observability matrix fourth round
Lie3 = obs_row3 * X_dot;
obs_row4 = jacobian(Lie3,X);


observability_matrix = [obs_row1;obs_row2;obs_row3;obs_row4];
rank(observability_matrix)