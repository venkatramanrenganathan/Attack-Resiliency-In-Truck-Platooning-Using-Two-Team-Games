clear all; close all; clc;
%% Dynamical model
% Trucks (except leader)
N = 5;
% Continuous time model
Ac = [zeros(N), eye(N); zeros(N), zeros(N)];
Bc = [zeros(N); diag(ones(N-1,1), -1) - eye(N)];

% Sampling time
Ts = 0.2;

% Discrete time model
A = eye(2*N) + Ts*Ac;
B = Ts*Bc;

%% One?player LQR
% Cost function
% J = x(K)' Qf x(K) + sum k=0ˆ(K?1) (x(k)' Q x(k) + u(k)' R u(k))
Kc = 3;
K = Kc / Ts;
Qf = 100 * eye(2*N);
Q = eye(2*N);
R = 0.01 * eye(1);
%% Optimal gain
P = nan(2*N,2*N,K+1);
Gamma = nan(N,2*N,K);

P(:,:,K+1) = Qf;
for k=K:-1:1
    P(:,:,k) = Q + A' * P(:,:,k+1) * A - ...
               A' * P(:,:,k+1) * B * inv(kron(R,eye(N)) + ...
               B' * P(:,:,k+1) * B) * B' * P(:,:,k+1) * A;
    Gamma(:,:,k) = - inv(kron(R,eye(N)) + ...
                   B' * P(:,:,k+1) * B) * B' * P(:,:,k+1) * A;
end

figure
plot(squeeze(reshape(P,[(2*N)^2,1,K+1]))')
grid on
title('one player LQR - P matrix')

figure
plot(squeeze(reshape(Gamma,[2*(N)^2,1,K]))')
grid on
title('one player LQR - \Gamma matrix')

%% Optimal gain (steady state)
Pbar = dare(A,B,Q,kron(R,eye(N)),zeros(2*N,N), eye(2*N));
Gammabar = - inv(R + B' * Pbar * B) * B' * Pbar * A;

figure
hold on
scatter(ones((2*N)^2,1),Pbar(:))
hold off

figure
hold on
scatter(ones(2*N^2,1),Gammabar(:))
hold off

% optimal feedback law
T = K+5;
axt = 1:T;
x = nan(2*N,T);
x(:,1) = [zeros(N,1); -0.5 *ones(N,1)];
xbar = nan(2*N,T);
xbar(:,1) = [zeros(N,1); -0.5 *ones(N,1)];
u = zeros(N,T);
ubar = zeros(N,T);

for t = 1:T-1
    if t <= K
        u(:,t) = Gamma(:,:,t) * x(:,t);
        ubar(:,t) = Gammabar * xbar(:,t);
    end
    x(:,t+1) = A * x(:,t) + B * u(:,t);
    xbar(:,t+1) = A * xbar(:,t) + B * ubar(:,t);
end

figure(3)
plot(1:T,x(1:N,:),'*',1:T,u, '--' )
grid on
xlabel('time')
title('Position and control (acceleration)')

% Truck trajectory
w = (diag(ones(N,1), -1) - eye(N+1)) \ ...
    [-ones(1,T), -0.5*ones(1,T); zeros(N,T), x(N+1:end,:)];
y = -ones(2*T,1)*(0:N) + cumsum(w');
figure
plot(1:2*T,y);
grid on
title('Truck trajectories')
xlim([T-5 2*T])
xlabel('time')

%% N?player LQR
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
F = zeros(N);
G = zeros(N,2*N);
S = zeros(2*N);

for k=K:-1:1
    for player=1:N
        Phi = inv( R + B(:,player)' * PN(:,:,k+1,player) * B(:,player)) * B(:,player)' * PN(:,:,k+1,player);
        for pl2 = 1:N
            F(player,pl2) = Phi * B(:,pl2);
            G(player,:) = -Phi * A;
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

figure
for player=1:min(N,6)
    subplot(2,3,player)
    plot(squeeze(reshape(PN(:,:,:,player),[(2*N)^2,1,K+1]))');
end
subplot(2,3,2)
title('N player LQR - P matrix of the first 6 players')

figure(12)
for player=1:min(N,6)
    subplot(2,3,player)
    plot(squeeze(reshape(GammaN(player,:,:),[2*N,1,K]))');
end
subplot(2,3,2)
title('N player LQR - \Gamma matrix of the first 6 players')

%% optimal feedback law
xN = nan(2*N,T);
xN(:,1) = [zeros(N,1); -0.5 *ones(N,1)];
uN = zeros(N,T);
for t = 1:T-1
    if t <= K
        uN(:,t) = GammaN(:,:,t) * xN(:,t);
    end
    xN(:,t+1) = A * xN(:,t) + B * uN(:,t);
end

figure
plot(1:T,xN(1:N,:),'*',1:T,uN, '--' )
grid on
xlabel('time')
title('Position and control (acceleration)')

% Truck trajectory
%wN = [diag(ones(N,1), -1) - eye(N+1)] \ [- [ones(1,T), 0.5*ones(1,T)]; [zeros(N,T), xN(N+1:end,:)]];
wN = inv([diag(ones(N,1), -1) - eye(N+1)]) * [- [ones(1,T), 0.5*ones(1,T)]; [zeros(N,T), xN(N+1:end,:)]];
yN = -ones(2*T,1)*(0:N) + cumsum(wN');

figure(14)
plot(1:2*T,yN);
grid on
title('Truck trajectories (Nash equilibrium)')
xlim([T-5 2*T])
xlabel('time')