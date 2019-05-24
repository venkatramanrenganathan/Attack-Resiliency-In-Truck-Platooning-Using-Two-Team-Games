%
% GameSec 2018 Conference
% Venkatraman Renganathan
% UT Dallas - Last Edited on 11th June, 2018
% 

%% Get Attacker Position Details
clear all; close all; clc;
N                                  = 5;                                           % Number of total following trucks
truck_vector                       = 1:N;
total_attackers                    = 3;
[attacker_position, gamma_record]  = get_attacker_positions_Ver2(N);
bad_trucks                         = sort(attacker_position(1:total_attackers),1);
critical_gamma                     = max(gamma_record(1:total_attackers));         % parameter quantifying attacker's capability
min_gamma                          = critical_gamma + 1;
max_gamma                          = 100000;
gamma                              = max_gamma;    % We can set it to max_gamma or min_gamma
bad_trucks                         = sort(bad_trucks);
truck_vector(bad_trucks)           = 0;
good_trucks                        = truck_vector(truck_vector ~= 0);

%% Dynamical model of Platooning Trucks (except leader)
% Continuous time model
%Ac   = [zeros(N), eye(N); zeros(N), zeros(N)];
Ac   = diag(ones(2*N-1,1),-1) - eye(2*N);
Bc   = [zeros(N); diag(ones(N-1,1), -1) - eye(N)];
B_k  = Bc;
F_k  = Bc;
   
for i = 1:length(good_trucks)    
    F_k(:,good_trucks(i)) = zeros(2*N,1);
end
for i = 1:length(bad_trucks)
    B_k(:,bad_trucks(i)) = zeros(2*N,1); 
end

% Sampling time
T_sample = 0.5;

% Discrete time model
A   = eye(2*N) + T_sample * Ac;
B   = T_sample * Bc;
B_k = T_sample * B_k;
F_k = T_sample * F_k;

% Cost function parameters
Kc    = 3;
K     = Kc / T_sample;
rho   = 0.01;              % penalty scaling constant
Qf    = 100 * eye(2*N);    % Terminal State Penalty Matrix 
Q     = eye(2*N);          % State Penalty Matrix
R_u   = rho * eye(N);      % Defender Control Penalty Matrix
R_v   = rho * eye(N);      % Attacker Control Penalty Matrix
D     = 1;                 % Desired Safe inter-vehicular distance
E     = D.*ones(1,2*N);    % Vector specifying the distance offset
T     = K + 5;             % Total Simulation Steps

%%
% Zero Sum Case - Value Function is of the form x^T P_t x + 2 q_t' x + r_t
%
 
QfN = zeros(2*N,2*N);
QN  = zeros(2*N,2*N);

for player = 1:N
    player_states = [player, player+N];
    QfN(player_states,player_states) = Qf(player_states,player_states);
    QN(player_states,player_states) = Q(player_states,player_states);
end

PN = nan(2*N,2*N,K+1);
qN = nan(2*N,K+1);
rN = nan(K+1);

PN(:,:,K+1) = QfN;             % P_T = Q_T
qN(:,K+1)   = zeros(2*N,1);    % q_T = 0
rN(K+1)     = 0;               % r_T = 0

theta     = nan(N,2*N,K);
phi       = nan(N,2*N,K);
lambda_1  = nan(N,K);
lambda_2  = nan(N,K);
alpha     = nan(N,N,K);
alpha_inv = nan(N,N,K);
beta      = nan(N,N,K);
beta_inv  = nan(N,N,K);

for k=K:-1:1   
    
    alpha(:,:,k)     = R_u + B_k' * PN(:,:,k+1) * B_k;        
    alpha_inv(:,:,k) = inv(alpha(:,:,k));

    beta(:,:,k)      = -gamma^2 * R_v + F_k' * PN(:,:,k+1) * F_k;
    beta_inv(:,:,k)  = inv(beta(:,:,k));

    mu_k             = alpha_inv(:,:,k) * B_k' * PN(:,:,k+1) * F_k * beta_inv(:,:,k) * F_k';
    zeta_k           = beta_inv(:,:,k) * F_k' * PN(:,:,k+1) * B_k * alpha_inv(:,:,k) * B_k';
    kappa            = eye(N) - mu_k * PN(:,:,k+1) * B_k;
    eta              = eye(N) - zeta_k * PN(:,:,k+1) * F_k;

    theta(:,:,k)     = kappa \ ((mu_k * PN(:,:,k+1) - alpha_inv(:,:,k) * B_k' * PN(:,:,k+1)) * A);
    lambda_1(:,k)    = kappa \ ((mu_k - alpha_inv(:,:,k) * B_k') * qN(:,k+1));

    phi(:,:,k)       = eta \ ((zeta_k * PN(:,:,k+1) - beta_inv(:,:,k) * F_k' * PN(:,:,k+1)) * A);
    lambda_2(:,k)    = eta \ ((zeta_k - beta_inv(:,:,k) * F_k') * qN(:,k+1));

    % Calculation of recursion coefficients

    PN(:,:,k) = Q + A' * PN(:,:,k+1) * A + theta(:,:,k)'*alpha(:,:,k)*theta(:,:,k) + ...
                phi(:,:,k)'*beta(:,:,k)*phi(:,:,k) + 2*theta(:,:,k)'*B_k'*PN(:,:,k+1)*A + ...
                2*phi(:,:,k)'*F_k'*PN(:,:,k+1)*A + theta(:,:,k)'*B_k'*PN(:,:,k+1)*F_k*phi(:,:,k) + ...
                phi(:,:,k)'*F_k'*PN(:,:,k+1)*B_k*theta(:,:,k);

    qN(:,k)   = (lambda_1(:,k)'*alpha(:,:,k)*theta(:,:,k) + lambda_2(:,k)' * beta(:,:,k)*phi(:,:,k) + ...
                lambda_1(:,k)'*B_k'*PN(:,:,k+1) * (A + F_k * phi(:,:,k)) + ...
                lambda_2(:,k)'*F_k'*PN(:,:,k+1) * (A + B_k * theta(:,:,k)) + ...
                qN(:,k+1)'*(A + B_k * theta(:,:,k) + F_k * phi(:,:,k)) - E)';

    rN(k)     = rN(k+1) + lambda_1(:,k)'*alpha(:,:,k)*lambda_1(:,k) + lambda_2(:,k)' * beta(:,:,k)*lambda_2(:,k) + ...
                lambda_1(:,k)'*B_k'*PN(:,:,k+1) * F_k * lambda_2(:,k) + ...
                lambda_2(:,k)'*F_k'*PN(:,:,k+1) * B_k * lambda_1(:,k) + ...
                2*K*D^2 + 2*qN(:,k+1)'*(B_k * lambda_1(:,k) + F_k * lambda_2(:,k));
    
end

%% optimal feedback law
xN      = nan(2*N,T);
xN(:,1) = [zeros(N,1); -0.1 *ones(N,1)];
uN      = zeros(N,T);
vN      = zeros(N,T);

% Defender input, u_k     = theta * x + lambda_1
% Attacker input, v_k     = phi   * x + lambda_2
% State Update,   x_{k+1} = A     * x + B_k * u_k + F_k * v_k

for t = 1:T-1
    if t <= K
        uN(:,t) = theta(:,:,t) * xN(:,t) + lambda_1(:,t); 
        vN(:,t) = phi(:,:,t) * xN(:,t) + lambda_2(:,t);
    end 
    xN(:,t+1) = A * xN(:,t) + B_k * uN(:,t) + F_k * vN(:,t);   
end

% Plot only corresponding controls
uN(bad_trucks,:)  = []; 
vN(good_trucks,:) = [];

% figure
% plot(1:T,xN(N,:))
% grid on
% xlabel('time')
% title('Distance to The Preceding Truck')
% a = findobj(gcf, 'type', 'axes');
% h = findobj(gcf, 'type', 'line');
% set(h, 'linewidth', 4);
% set(a, 'linewidth', 4);
% set(a, 'FontSize', 24);


figure
plot(1:T,uN)
grid on
xlabel('time')
title('Good Trucks Control (acceleration)')
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);

figure
plot(1:T,vN)
grid on
xlabel('time')
title('Bad Trucks Control (acceleration)')
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);

% Truck trajectory
wN = [diag(ones(N,1), -1) - eye(N+1)] \ [-ones(1,2*T); [zeros(N,T), xN(N+1:end,:)]];
yN = -ones(2*T,1)*(0:N) + cumsum(wN');

figure
f1 = plot(1:2*T,yN(:,1), 'k');
hold on
f2 = plot(1:2*T,yN(:,good_trucks+1), 'b');
hold on
f3 = plot(1:2*T,yN(:,bad_trucks+1), 'r');
hold off
grid on
f = [f1;f2(1);f3(1)];
legend(f,'Leader', 'Good Trucks', 'Bad Trucks');
title('Truck trajectories (Nash equilibrium)');
xlim([0 2*T])
xlabel('time')
ylabel('truck position')
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);