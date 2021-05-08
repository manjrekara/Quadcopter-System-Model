% ECE 3111 Final Project Matlab File
% by Atharva Manjrekar

% Collaboraters: 
%Ajani Blackwood
%Tariq Ali
%David Szymanowski

clc
clear;
close all;
%% Problem 1
%Symbolic variables

% XYZ are frame coordinates
syms x y z;

% Linear Velocities
syms Vx Vy Vz;

% Accelearation V'(x), V'(y), V'(z)
syms Vx_dot Vy_dot Vz_dot;

% Angular Velocities
syms Wx Wy Wz;

% Angular Acceleratio
syms Wx_dot Wy_dot Wz_dot;

% Roll Pitch Yaw
syms Roll Pitch Yaw;

% System Rotor Voltages
syms Vr1 Vr2 Vr3 Vr4;

% Constants for mass and other parameters
Mass = 0.5;             %Unit:kg          
L = 0.25;               %Unit:m           %symbol L
K = 3e-6;               %symbol K
Kd = 0.25;              %Unit Kg/s        %symbol Kd
Drag_Coeff = 1e-7;      %Unit Nms^2       %symbol b

g = 9.81;               %Unit m/s^2       %symbol g
Ixx = 5e-3;             %Unit Kgm^2       %Ixx
Iyy = Ixx;              %Iyy
Izz = 1e-2;             %Unit             %Izz
Cm = 1e4;               %Unit V^-2 S^-2   %Cm

% U values
% We know U1+U2+U3+U4 = (Vr1^2)+(Vr2^2)+(Vr3^2)+(Vr4^2)
syms U1 U2 U3 U4 
% U1 = Vr1^2;
% U2 = Vr2^2;
% U3 = Vr3^2;
% U4 = Vr4^2;

% Linearized EOM
Xdot = Vx;                                                  %Equation 5

Ydot = Vy;                                                  %Equation 6

Zdot = Vz;                                                  %Equation 7

Vx_dot = (-Kd/Mass)*Vx + (g*Pitch);                         %Equation 8

Vy_dot = (-Kd/Mass)*Vy - (g*Roll);                          %Equation 9

% Vz_dot = (-Kd/Mass)*Vz - g + (K*Cm/Mass)*(U1+U2+U3+U4);    
% 
% Vz_dot = (-Kd/Mass)*Vz - 9.81 + ((3E-6)*(1E4)/0.5)*(163.5);     
% 
% Vz_dot = (-Kd/Mass)*Vz - 9.81 + 9.81;                      

Vz_dot = (-Kd/Mass)*Vz;                                     %Equation 10

Roll_dot = Wx;                                              %Equation 11

Pitch_dot = Wy;                                             %Equation 12

Yaw_dot = Wz;                                               %Equation 13
    
Wx_dot = ((L*K*Cm) / Ixx) * (U1-U3);                        %Equation 14 

Wy_dot = ((L*K*Cm) / Iyy) * (U2-U4);                        %Equation 15

Wz_dot = ((L*K*Cm) / Izz) * (U1-U2-U3-U4);                  %Equation 16

%% Problem 2 Part 1
%Matrix pieces

Mat_X = transpose([x y z Vx Vy Vz Roll Pitch Yaw Wx Wy Wz])

Mat_A = [0   0    0    1    0    0    0    0    0    0    0    0    % X dot
         0   0    0    0    1    0    0    0    0    0    0    0    % Y dot
         0   0    0    0    0    1    0    0    0    0    0    0    % Z dot
         0   0    0 -0.5    0    0    0    g    0    0    0    0    % Vx dot
         0   0    0    0 -0.5    0   -g    0    0    0    0    0    % Vy dot
         0   0    0    0    0 -0.5    0    0    0    0    0    0    % Vz dot  
         0   0    0    0    0    0    0    0    0    1    0    0    % Roll dot
         0   0    0    0    0    0    0    0    0    0    1    0    % Pitch dot
         0   0    0    0    0    0    0    0    0    0    0    1    % Yaw dot
         0   0    0    0    0    0    0    0    0    0    0    0    % Wx dot
         0   0    0    0    0    0    0    0    0    0    0    0    % Wy dot
         0   0    0    0    0    0    0    0    0    0    0    0]   % Wz dot

fprintf('The size of the matrix A is: %s', num2str(size(Mat_A)))

Mat_U = transpose([U1 U2 U3 U4]);
                       
Mat_B = [ 0     0     0     0           % X dot
          0     0     0     0           % Y dot
          0     0     0     0           % Z dot
          0     0     0     0           % Vx dot
          0     0     0     0           % Vy dot
        .06    .06   .06   .06          % Vz dot 
          0     0     0     0           % Roll dot
          0     0     0     0           % Pitch dot
          0     0     0     0           % Yaw dot
         1.5    0   -1.5    0           % Wx dot
          0    1.5    0   -1.5          % Wy dot
         .75  -.75   .75  -.75 ]        % Wz dot
                       
fprintf('The size of the matrix B is: %s', num2str(size(Mat_B)))

Mat_Xdot = (Mat_A*Mat_X) + (Mat_B*Mat_U)
fprintf('The size of the matrix Xdot is: %s', num2str(size(Mat_Xdot)))
             
%Output Matrix y = transpose([x y z φ θ ψ])

%Mat_Y = transpose([x y z Roll Pitch Yaw])

Mat_C = [1 0 0 0 0 0 0 0 0 0 0 0            % X dot
         0 1 0 0 0 0 0 0 0 0 0 0            % Y dot
         0 0 1 0 0 0 0 0 0 0 0 0            % Z dot
         0 0 0 0 0 0 1 0 0 0 0 0            % Roll dot
         0 0 0 0 0 0 0 1 0 0 0 0            % Pitch dot
         0 0 0 0 0 0 0 0 1 0 0 0]           % Yaw dot

fprintf('The size of the matrix C is: %s', num2str(size(Mat_C)))
% Mat_Y = Mat_C * Mat_X
%% Problem 2 Part 2           

% Matrix GIJ
% Gij (s) = Ci(sI −A)−1Bj + D
syms s

SI_minusA = s* eye(12) - Mat_A;

Mat_Gij = (Mat_C *(SI_minusA)^-1 *(Mat_B))
fprintf('The size of the matrix Gij is: %s', num2str(size(Mat_Gij)))
%% Problem 3

% Four Rotor Voltages:
U = (9.8 / (K*Cm/Mass));            % U = 163.33 V
U1 = (9.8 / (K*Cm/Mass)) / 4;       % U1 = 40.833 V
U2 = U1 
U3 = U1
U4 = U1

V1 = sqrt(U1)                       % V1 = 6.39 V
V2 = sqrt(U2)
V3 = sqrt(U3)
V4 = sqrt(U4)

syms Zs s Vzs Us

% Equation #2
Vzs = s*Zs;

% Equation #3
Vz_dot = -0.5 * Vzs + (K*Cm / Mass) * U1;

% Equation #2 into #3
Us = (1/(4*(K*Cm / Mass))) * (Vz_dot + 0.5*Vzs)
Zs = Us * (4*(K*Cm / Mass) / (s^2 + 0.5*s))

Transfer_func = (expand(Zs / Us))

% Transffer function to plot

TF = tf([0.24], [1 0.5 0])

%% Problem 4 PID controller

figure(1)
margin(TF)

figure(2)
rlocus(TF)

figure(3)
step(TF)









