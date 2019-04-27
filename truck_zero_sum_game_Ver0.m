clear all; close all; clc;
%% Dynamical model
% Trucks (except leader)
N = 5;
good_N = [1 3 5];
bad_N = [2 4];
% Continuous time model
Ac = [zeros(N), eye(N); zeros(N), zeros(N)];
Bc = [zeros(N); diag(ones(N-1,1), -1) - eye(N)];
B_k = Bc;
F_k = Bc;
B_k(:,2) = zeros(2*N,1);
B_k(:,4) = zeros(2*N,1);
F_k(:,1) = zeros(2*N,1);
F_k(:,3) = zeros(2*N,1);
F_k(:,5) = zeros(2*N,1);

% Sampling time
Ts = 0.2;

% Discrete time model
A = eye(2*N) + Ts*Ac;
B = Ts*Bc;
B_k = Ts*B_k;
F_k = Ts*F_k;

% Cost function
% J = x(K)' Qf x(K) + sum k=0ˆ(K?1) (x(k)' Q x(k) + u(k)' R u(k))
Kc = 3;
K = Kc / Ts;
Qf = 100 * eye(2*N);
Q = eye(2*N);
R = 0.01 * eye(1);
T = K + 5;

%%
% Zero Sum Case
%

QfN = zeros(2*N,2*N,N);
QN = zeros(2*N,2*N,N);

for player = 1:N
    player_states = [player, player+N];
    QfN(player_states,player_states,player) = Qf(player_states,player_states);
    QN(player_states,player_states,player) = Q(player_states,player_states);
end

PN = nan(2*N,2*N,K+1,N);

for player = 1:N
    PN(:,:,K+1,player) = QfN(:,:,player);
end

GammaN = nan(N,2*N,K);
Phi = nan(2,2*N);
Psi = nan(2,2*N);
F = zeros(N);
G = zeros(N,2*N);
S = zeros(2*N);

for k=K:-1:1
    
    for player=1:N
        
        Phi = inv( R + B_k(:,player)' * PN(:,:,k+1,player) * B_k(:,player)) * B_k(:,player)' * PN(:,:,k+1,player);
        Psi = inv( -R + F_k(:,player)' * PN(:,:,k+1,player) * F_k(:,player)) * F_k(:,player)' * PN(:,:,k+1,player);
        
        for pl2 = 1:length(good_N)
            F(player,good_N(pl2)) = Phi * B_k(:,good_N(pl2));
            G(player,:) = -Phi * A;
        end
        
        for pl2 = 1:length(bad_N)
            F(player,bad_N(pl2)) = Psi * F_k(:,bad_N(pl2));
            G(player,:) = -Psi * A;
        end
        
        F(player,player) = 1;
        
    end
    
    GammaN(:,:,k) = F\G;
    
    for player = 1:N
        S = A;
        for pl2 = 1:N
            S = S + B(:,pl2) * GammaN(pl2,:,k);
        end
        PN(:,:,k,player) = QN(:,:,player) + GammaN(player,:,k)' * R * GammaN(player,:,k) + ...
        S' * PN(:,:,k+1,player) * S;
    end
end

%% optimal feedback law
xN = nan(2*N,T);
xN(:,1) = [zeros(N,1); -0.5 *ones(N,1)];
uN = zeros(N,T);
good_uN = uN;
bad_uN = uN;
for t = 1:T-1
    if t <= K
        uN(:,t) = GammaN(:,:,t) * xN(:,t); % u = kx + Fv, Fv not coded yet                                            
    end 
    xN(:,t+1) = A * xN(:,t) + B * uN(:,t);
end


figure
plot(1:T,xN(1:N,:))
grid on
xlabel('time')
legend('Truck 1', 'Truck 2', 'Truck 3', 'Truck 4','Truck 5');
title('Position')

figure
plot(1:T,uN)
grid on
xlabel('time')
legend('Truck 1', 'Truck 2', 'Truck 3', 'Truck 4','Truck 5');
title('Control (acceleration)')

% Truck trajectory
%wN = [diag(ones(N,1), -1) - eye(N+1)] \ [- [ones(1,T), 0.5*ones(1,T)]; [zeros(N,T), xN(N+1:end,:)]];
wN = inv([diag(ones(N,1), -1) - eye(N+1)]) * [- [ones(1,T), 0.5*ones(1,T)]; [zeros(N,T), xN(N+1:end,:)]];
yN = -ones(2*T,1)*(0:N) + cumsum(wN');

figure
plot(1:2*T,yN);
grid on
title('Truck trajectories (Nash equilibrium)')
xlim([T-5 2*T])
xlabel('time')