function [attacker_position] = get_attacker_positions_Ver0()

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
attacker_position = zeros(N,1);            % final order of attackers
gamma_record      = NaN(N,1);
F_base            = zeros(N);

for player_1 = 1:N
    
    if (player_1 > 1)
        F_base(attacker_position(player_1-1),player_1-1) = 1;
    end
    
    F = F_base;
    
    for player_2 = 1:N
        
        if (sum(F(player_2,:)) ~= 0)
            continue
        end
        
        F(player_2,player_1) = 1;
        check_flag = 1;
        
        while (check_flag > 0)
            
            gamma(player_2,player_1) = (gamma_low(player_2,player_1) + gamma_high(player_2,player_1))/2;            
            P = Q;
            attacker_check_flag = 0;
            defender_check_flag = 0;
            
            for time = 1:T
                
                alpha     = R_u + B' * P * B;        
                alpha_inv = inv(alpha);
                beta      = -gamma(player_2,player_1)^2 * R_v + F' * P * F;
                beta_inv  = inv(beta);
                mu_k      = alpha \ (B' * P * F * (beta \ F'));
                zeta_k    = beta \ ( F' * P * B * (alpha \ B'));
                kappa_inv = inv(eye(N) - mu_k * P * B);
                eta_inv   = inv(eye(N) - zeta_k * P * F);
                theta     = kappa_inv * (mu_k * P - alpha_inv * B' * P) * A;
                phi       = eta_inv * (zeta_k * P - beta_inv * F' * P) * A;                                
                P         = Q + A'*P*A + theta'*alpha*theta + phi'*beta*phi + ...
                            2*theta'*B'*P*A + 2*phi'*F'*P*A + ...
                            theta'*B'*P*F*phi + phi'*F'*P*B*theta;
%                 P         = Q + A'*P*inv(eye(N) + (B*B' - (gamma(player_2,player_1)^-2)*F*F')*P)*A;
                
                if (min(eig(gamma(player_2,player_1)^2 * R_v - F'*P*F)) <= 0)
                    attacker_check_flag = 1;
                    break
                end
                if min(eig(R_u + B'*P*B)) <= 0
                    defender_check_flag = 1;
                    break
                end
                
            end
            
            if ((isnan(trace(P)) == 1) || (attacker_check_flag == 1) || (defender_check_flag == 1))
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
end



end

