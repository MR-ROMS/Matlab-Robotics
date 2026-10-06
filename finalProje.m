%ROMARIC SILANANG DOROBA : 2110225556
%WİS ALİBRAHİM 20010225561
clc; clear; close all;
cube_center = [0.55, 0, 0.55];
side        = 0.3;
h_side      = side/2;

run('YASKAWA_GP10_DataFile.m');
robot = importrobot('YASKAWA_GP10.slx');
robot.DataFormat = 'column';

ik      = inverseKinematics("RigidBodyTree", robot);
weights = [0.25 0.25 0.25 1 1 1];
inGuess = homeConfiguration(robot);

%Points
xc = cube_center(1); yc = cube_center(2); zc = cube_center(3);
P1 = [xc+h_side, yc-h_side, zc+h_side];  
P2 = [xc+h_side, yc+h_side, zc+h_side];   
P3 = [xc-h_side, yc+h_side, zc+h_side];  
P4 = [xc-h_side, yc-h_side, zc+h_side];   

%Orientations
R_down = roty(180); 

R12 = R_down *roty(-45);
R23 = R_down * rotx(45);
R34 = R_down * roty(45);
R41 = R_down * rotx(45);

% Waypoints
T_P1 = [R12, P1'; 0 0 0 1];   
T_P2 = [R_down, P2'; 0 0 0 1];   
T_P3 = [R34, P3'; 0 0 0 1];   
T_P4 = [R41, P4'; 0 0 0 1];  

%Time Setup
dt = 0.05;
t1 = 0:dt:3;     
t2 = 3:dt:5;     
t3 = 5:dt:6;    
t4 = 6:dt:8;     
t5 = 8:dt:9;  
t6 = 9:dt:11;    
t7 = 11:dt:12; 
t8 = 12:dt:14;  
t9 = 14:dt:16; 

%Waypoints
[q_P1, ~] = ik('Body6', T_P1, weights, inGuess); q_P1 = unwrapJoints(q_P1, inGuess);
[q_P2, ~] = ik('Body6', T_P2, weights, q_P1);    q_P2 = unwrapJoints(q_P2, q_P1);
[q_P3, ~] = ik('Body6', T_P3, weights, q_P2);    q_P3 = unwrapJoints(q_P3, q_P2);
[q_P4, ~] = ik('Body6', T_P4, weights, q_P3);    q_P4 = unwrapJoints(q_P4, q_P3);

%Joint Space trajectories
[qseg1, qdseg1, qddseg1] = cubicpolytraj([inGuess q_P1], [t1(1) t1(end)], t1);
[qseg3, qdseg3, qddseg3] = cubicpolytraj([q_P2    q_P4], [t3(1) t3(end)], t3);
[qseg5, qdseg5, qddseg5] = cubicpolytraj([q_P1    q_P3], [t5(1) t5(end)], t5);
[qseg7, qdseg7, qddseg7] = cubicpolytraj([q_P4    q_P2], [t7(1) t7(end)], t7);
[qseg9, qdseg9, qddseg9] = cubicpolytraj([q_P3 inGuess], [t9(1) t9(end)], t9);

%Cartesian trajectories
[T12, v12, ~] = transformtraj(T_P1, T_P2, [t2(1) t2(end)], t2);
[T41, v41, ~] = transformtraj(T_P4, T_P1, [t4(1) t4(end)], t4);
[T34, v34, ~] = transformtraj(T_P3, T_P4, [t6(1) t6(end)], t6);
[T23, v23, ~] = transformtraj(T_P2, T_P3, [t8(1) t8(end)], t8);

%Inverse Kinematics and jacobians
[~,~,nsteps] = size(T12);
qseg2 = zeros(6,nsteps); qdseg2 = zeros(6,nsteps); qddseg2 = zeros(6,nsteps);
initialguess = q_P1;
for i = 1:nsteps
    [qseg2(:,i),~] = ik('Body6',T12(:,:,i),weights,initialguess);
    initialguess = qseg2(:,i);
    J = geometricJacobian(robot,qseg2(:,i),'Body6');
    qdseg2(:,i) = J\v12(:,i);
end

% P4 -> P1
[~,~,nsteps] = size(T41);
qseg4 = zeros(6,nsteps); qdseg4 = zeros(6,nsteps); qddseg4 = zeros(6,nsteps);
initialguess = q_P4;
for i = 1:nsteps
    [qseg4(:,i),~] = ik('Body6',T41(:,:,i),weights,initialguess);
    initialguess = qseg4(:,i);
    J = geometricJacobian(robot,qseg4(:,i),'Body6');
    qdseg4(:,i) = J\v41(:,i);
end

%P3 -> P4
[~,~,nsteps] = size(T34);
qseg6 = zeros(6,nsteps); qdseg6 = zeros(6,nsteps); qddseg6 = zeros(6,nsteps);
initialguess = q_P3;
for i = 1:nsteps
    [qseg6(:,i),~] = ik('Body6',T34(:,:,i),weights,initialguess);
    initialguess = qseg6(:,i);
    J = geometricJacobian(robot,qseg6(:,i),'Body6');
    qdseg6(:,i) = J\v34(:,i);
end

% P2 -> P3
[~,~,nsteps] = size(T23);
qseg8 = zeros(6,nsteps); qdseg8 = zeros(6,nsteps); qddseg8 = zeros(6,nsteps);
initialguess = q_P2;
for i = 1:nsteps
    [qseg8(:,i),~] = ik('Body6',T23(:,:,i),weights,initialguess);
    initialguess = qseg8(:,i);
    J = geometricJacobian(robot,qseg8(:,i),'Body6');
    qdseg8(:,i) = J\v23(:,i);
end

%Delete duplicates
qseg1(:,end)=[]; qdseg1(:,end)=[]; qddseg1(:,end)=[];
qseg2(:,end)=[]; qdseg2(:,end)=[]; qddseg2(:,end)=[];
qseg3(:,end)=[]; qdseg3(:,end)=[]; qddseg3(:,end)=[];
qseg4(:,end)=[]; qdseg4(:,end)=[]; qddseg4(:,end)=[];
qseg5(:,end)=[]; qdseg5(:,end)=[]; qddseg5(:,end)=[];
qseg6(:,end)=[]; qdseg6(:,end)=[]; qddseg6(:,end)=[];
qseg7(:,end)=[]; qdseg7(:,end)=[]; qddseg7(:,end)=[];
qseg8(:,end)=[]; qdseg8(:,end)=[]; qddseg8(:,end)=[];

q   = [qseg1 qseg2 qseg3 qseg4 qseg5 qseg6 qseg7 qseg8 qseg9];
qd  = [qdseg1 qdseg2 qdseg3 qdseg4 qdseg5 qdseg6 qdseg7 qdseg8 qdseg9];
qdd = [qddseg1 qddseg2 qddseg3 qddseg4 qddseg5 qddseg6 qddseg7 qddseg8 qddseg9];

N     = size(q, 2);
t_vec = (0:N-1)' * dt;
qi    = [t_vec, q', qd', qdd'];

Kp = 50000;
Kd = 8000;

assignin('base', 'qi',    qi);
assignin('base', 'Kp',    Kp);
assignin('base', 'Kd',    Kd);
assignin('base', 'robot', robot);
assignin('base', 'cube_center', cube_center);
assignin('base', 'side', side);

disp("Simulation starting soon...");
sim("YASKAWA_GP10cont.slx", N*dt);