% Ver 1_1: 
%   - Based on car model
%
%
function dstate = CarSys_Ver1_1(t,state,par)

% Parameters passed down to the ODE solver
n       = par.n;            % Number of cars
ds      = par.ds;           % Desired distance to the car in front
p       = par.p;            % Desired travel direction 
L       = par.L;            % Car wheelbase
kv      = par.kv;           % Linear velocity control gain
kw      = par.kw;           % Angular velocity control gain


%% Simulation parameters

q = state(1:2*n);           % Position vector
theta = state(2*n+1: 3*n);  % Steering vector
phi = state(3*n+1: end);    % Heading vector

G  = zeros(2*n, n);         % Steering matrix
Gp = zeros(2*n, n);         % Perpendicular steering matrix
H  = zeros(2*n, n);         % Heading matrix
Hp = zeros(2*n, n);         % Perpendicular heading matrix
R  = [0 -1; 1 0];           % 90 degree roation matrix

delta = theta + phi;        % Sum of heading and steering angles

% Steering vector
g = [cos(delta).'; sin(delta).'];
g = g(:);

% Heading vector
h = [cos(theta).'; sin(theta).'];
h = h(:);

% Steering and heading matrices, and their perpendiculars
for i = 1 : n    
    G(2*i-1:2*i,i)  =  g(2*i-1:2*i);
    Gp(2*i-1:2*i,i) =  R * g(2*i-1:2*i);
    H(2*i-1:2*i,i)  =  h(2*i-1:2*i);
    Hp(2*i-1:2*i,i) =  R * h(2*i-1:2*i);
end

PhiS = sin( diag(phi) );


%% Platoon control

% Preallocate control vector
u = zeros(2*n,1);

for i = 1 : n
    if i == 1 % Set control for agent 1 along the desired direction of motion
        u(2*i-1:2*i) = p;
    else % For other agents go to the desired point behind the leading neighbor
        u(2*i-1:2*i) = (q(2*i-3:2*i-2) - ds .* p) - q(2*i-1:2*i);
    end
end


%% Car control

v       = kv .* G.' * u;        % Wheel linear speed
omega   = kw .* Gp.' * u;       % Wheel steering speed

dq      = G * v;
dtheta  = (1/L) * PhiS * v;
dphi    = omega;
dstate  = [dq; dtheta; dphi];   % Derivative of state


end





































































