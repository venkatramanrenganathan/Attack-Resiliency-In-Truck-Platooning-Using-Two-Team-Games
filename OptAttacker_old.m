%%Defenders assumed at all nodes
%%Attackers assigned by greedy algorithm to location with highest crit
%%gamma corresponding to best position for attack

clear all; close all; clc;

N = 10;      % number of states
T = 100;     % time period. Set to ~10000 for infinite horizon approximation
Q = eye(N);  % state cost
R_u = eye(N);  % attacker input cost
R_v = eye(N);  % attacker input cost


%System Dynamics Matrix - line graph
A_base = eye(N) + diag(ones(N-1,1),-1) + diag(ones(N-1,1),1);
A      = A_base /(0.1 + max(abs(eig(A_base))));

%Input Matrix - assumes input to all nodes
B = ones(N); 

%Calculation initialization
gamma_low         = zeros(N);    
gamma_initial     = 100;                   % gamma initial value
gamma_high        = gamma_initial*ones(N);
gamma             = zeros(N);
attacker_position = zeros(N,1);            %final order of attackers
gamma_record      = NaN(N,1);
F_base            = zeros(N);

for player_1 = 1:N
    
    for player_2 = 1:N
    
        F = F_base;
        
        if sum(F(player_2,:))~=0
            continue
        end
        
        F(player_2,player_1) = 1;
        check_flag = 1;
        
        while (check_flag > 0)
            
            gamma(player_2,player_1) = (gamma_low(player_2,player_1) + gamma_high(player_2,player_1))/2;
            Z = Q;
            attacker_check_flag = 0;
            defender_check_flag = 0;
            
            for time = 1:T
                
                Z = Q + A'*Z*inv(eye(N) + ...
                    (B*B' - (gamma(player_2,player_1)^-2)*F*F')*Z)*A;
                
                if (min(eig(gamma(player_2,player_1)^2 * R_v - F'*Z*F)) <= 0)
                    attacker_check_flag = 1;
                    break
                end
                if min(eig(R_u + B'*Z*B)) <= 0
                    defender_check_flag = 1;
                    break
                end
                
            end
            
            if ((isnan(trace(Z)) == 1) || (attacker_check_flag == 1) || (defender_check_flag == 1))
                gamma_low(player_2,player_1) = gamma(player_2,player_1);
            else
                if (abs(gamma_high(player_2,player_1) - gamma(player_2,player_1)) < 0.0001)
                    check_flag = 0;
                end
                if (abs(gamma_low(player_2,player_1)-gamma_initial) < 0.0001)
                    check_flag = -1;
                end
                gamma_high(player_2,player_1) = gamma(player_2,player_1);
            end
            
        end
        
        if (attacker_position(player_1) == 0)
            attacker_position(player_1) = player_2;
            gamma_record(player_1) = gamma(player_2,player_1);            
        elseif ((gamma_record(player_1) < gamma(player_2,player_1)) && (check_flag == 0))
            attacker_position(player_1) = player_2;
            gamma_record(player_1) = gamma(player_2,player_1);            
        end
        
    end  
    F_base(attacker_position(player_1),player_1)=1;
    gamma
end
attacker_position
figure
imagesc(gamma(:,:))
colorbar()
xlabel('Attacker number');
ylabel('Attacker position');
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);