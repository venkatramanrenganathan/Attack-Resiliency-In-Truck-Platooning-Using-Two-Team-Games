%
% ACC 2019 Conference
% Venkatraman Renganathan
% UT Dallas - Last Edited on 24th July, 2018
% 

%% Code to clear all before starting fresh 
clear all; close all; clc;

%% Set Truck Details

N             = 10;      % #of total following trucks
truck_vector  = 1:N;     % Vector for indexing each truck
mass          = 750;     % Mass of each truck - 750 kg
C_d           = 0.3;     % Drag coefficient 
A_f           = 1.3;     % Frontal area of the vehicle 
rho           = 1.2;     % Specific mass of air 
f_r           = 0.01;    % Rolling resistance coefficient
g             = 9.81;    % Gravity acceleration 
theta         = 0;       % Road Grade (elevation/descent)
v_ss          = 20;      % Steady state velocity 20 m/s
a_ss          = 0;       % Steady state acceleration 0 m/s^2
k_p           = 50;      % Proportional gain 
k_i           = 0;       % Integral gain 
k_d           = 700;     % Derivative gain


%% Calculating desired constant details

gamma         = C_d * A_f * rho;
F_roll        = mass * g * f_r;          % Rolling resistance force 
F_g           = mass * g * sin(theta);   % Gravity Force