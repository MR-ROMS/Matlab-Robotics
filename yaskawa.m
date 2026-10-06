clc; clear; close all;
cube_center = [0.65 0 0.7]; 
side = 0.3; h_side = side/2; dir = 1;
% ------------------ ROBOT SETUP ------------------
robot = importrobot("YASKAWA_GP10.slx");
robot.DataFormat = 'column';
ik = inverseKinematics("RigidBodyTree",robot);
inGuess = homeConfiguration(robot);

weights = [0.25 0.25 0.25 1 1 1];

%inGuess = homeConfiguration(robot);  
% ------------------ WAYPOINTS ------------------
xc = cube_center(1); yc = cube_center(2); zc = cube_center(3);
cube = [
    xc - h_side, yc - h_side*dir, zc + h_side;   % P1
    xc + h_side, yc - h_side*dir, zc + h_side;   % P2
    xc + h_side, yc - h_side*dir, zc - h_side;   % P3
    xc + h_side, yc + h_side*dir, zc - h_side;   % P4
    xc + h_side, yc + h_side*dir, zc + h_side;   % P5
];

 
% ------------------ ORIENTATIONS ------------------
R_i= rotx(180);  
%R = R_base * rotx(90*dir);   % torch pointing down
% R = R_base * rotx(180) * rotx(-90);

R_list = {
    roty(180),           % P1: top face   → tool points -Z (down)
    roty(90),           % P2: top face   → same
    roty(90),            % P3: side face  → tool points +X
    roty(90)*rotz(90),   % P4: side face  → tool points +Y  
    rotx(180),           % P5: top face   → back to top
};

% Build transforms for each cube corner
end_effector = 0.1; 
% T_offset = trvec2tform([0 0 tool_offset]);  % offset along Z
n = size(cube,1);
T_all = zeros(4,4,n);
 for i = 1:n
     R_i = R_list{i};
    normal_offset = [0 0 -end_effector];  % ← negative: back off from surface
    T_offset = trvec2tform(normal_offset);
    T_all(:,:,i) = trvec2tform(cube(i,:)) * rotm2tform(R_i) * T_offset;
end

% After fixing R, solve P1 first and check visually
% [sol, ~] = ik('Body6', T_all(:,:,1), weights, q_home);
% show(robot, sol);   % inspect before running the full trajectory


% ------------------ TRAJECTORY ------------------
dt = 0.05;
Q = [];
iguess = inGuess;  % your elbow-up seed

% 1) HOME → via → P1
[sol, info] = ik('Body6', T_all(:,:,1), weights, iguess);
q1 = sol; 
% q_via = inGuess;   % safe mid-posture

t1 = 0:dt:3;
q_seg1 = cubicpolytraj([inGuess, iguess, q1], [0 1.5 3], t1);
Q = [Q, q_seg1];
iguess = q1;

% Helper — use this pattern for ALL Cartesian segments
function q_seg = cartesian_seg(ik, body, T_traj, weights, iguess_in)
    n = size(T_traj, 3);
    q_seg = zeros(6, n);
    iguess = iguess_in;
    q_prev = iguess_in;
    for i = 1:n
        [sol,~] = ik(body, T_traj(:,:,i), weights, iguess);
        delta = sol - q_prev;
        sol = sol - round(delta/(2*pi)) * 2*pi;  % unwrap
        iguess = sol; q_prev = sol;
        q_seg(:,i) = sol;
    end
end

% 2) P1→P2
t2 = (3+dt):dt:4;
[T2,~,~] = transformtraj(T_all(:,:,1), T_all(:,:,2), [t2(1) t2(end)], t2);
q_seg2 = cartesian_seg(ik, 'Body6', T2, weights, iguess);
Q = [Q, q_seg2]; iguess = q_seg2(:,end);

% 3) P2→P3
t3 = (4+dt):dt:5;
[T3,~,~] = transformtraj(T_all(:,:,2), T_all(:,:,3), [t3(1) t3(end)], t3);
q_seg3 = cartesian_seg(ik, 'Body6', T3, weights, iguess);
Q = [Q, q_seg3]; iguess = q_seg3(:,end);

% 4) P3→P4
t4 = (5+dt):dt:6;
[T4,~,~] = transformtraj(T_all(:,:,3), T_all(:,:,4), [t4(1) t4(end)], t4);
q_seg4 = cartesian_seg(ik, 'Body6', T4, weights, iguess);
Q = [Q, q_seg4]; iguess = q_seg4(:,end);

% 5) P4→P5
t5 = (6+dt):dt:7;
[T5,~,~] = transformtraj(T_all(:,:,4), T_all(:,:,5), [t5(1) t5(end)], t5);
q_seg5 = cartesian_seg(ik, 'Body6', T5, weights, iguess);
Q = [Q, q_seg5]; iguess = q_seg5(:,end);

% 6) P5→HOME
t6 = (7+dt):dt:9;
q_seg6 = cubicpolytraj([Q(:,end), iguess, inGuess], [t6(1) mean([t6(1) t6(end)]) t6(end)], t6);
Q = [Q, q_seg6];

% ------------------ UNWRAP ENTIRE TRAJECTORY ------------------
for j = 1:6
    Q(j,:) = unwrap(Q(j,:));
end

t_total = 0:dt:(size(Q,2)-1)*dt;
qi = [t_total', Q'];

disp("Simulation starting soon...");
sim("YASKAWA_GP10.slx", t_total(end));
