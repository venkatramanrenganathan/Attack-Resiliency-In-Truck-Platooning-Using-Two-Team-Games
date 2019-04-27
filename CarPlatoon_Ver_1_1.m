% Ver 1_1:
%           - Base on simplified platoon control
%

addpath(genpath('Helpers'));


%% Simulation parameters

n       = 6;                  % Number of cars
q0      = rand(2*n,1);        % Random initial positions
phi0    = 2*pi*rand(n,1);     % Random initial steering angles
theta0  = 2*pi*rand(n,1);     % Random initial heading angles


%% Control parameters

vSat        = [3, 5];            % Allowed speed range
omegaSat    = [-pi/4, pi/4];     % Allowed heading angle rate of change
ds          = 0.2;               % Desired distance to the leading car
p           = [1; 0];            % Desired travel direction 
L           = 0.1;               % Wheelbase length
kv          = 1;                 % Linear velocity control gain
kw          = 3;                 % Angular velocity control gain


%% Simulate the model

% close all

T = [0, 100];     % Simulation time
simItr = 100;   % Simulation iterations 

Tvec = linspace(T(1), T(2), simItr); % Simulation time samples

% Parameters passed down to the ODE solver
par.n        = n;
par.ds       = ds;
par.p        = p;
par.L        = L;
par.kv       = kv; 
par.kw       = kw;

state0 = [q0; theta0; phi0];           % Initial state

% Simulate the ODE system
opt = odeset('AbsTol', 1.0e-05, 'RelTol', 1.0e-05);
[t,stateMat] = ode45(@CarSys_Ver1_1, Tvec, state0, opt, par);


%%

% % A rough plot of trajectories
% figure;
% hold on
% for i = 1 : n
%     scatter(stateMat(:,2*i-1), stateMat(:,2*i),10);
% end
% hold off
% axis equal
% 
% % A rough plot of final positions
% figure;
% hold on
% for i = 1 : n
%     scatter(stateMat(end,2*i-1), stateMat(end,2*i),50);
% end
% hold off
% axis equal


%% Make movie and plot the results

close all

% Simulation parameters
plotParam.N         = n;
plotParam.stateMat  = stateMat;
plotParam.L         = L;

% Make movie 
fileName = 'MovieCar';
MovieCarVer2_1_1(fileName, plotParam)