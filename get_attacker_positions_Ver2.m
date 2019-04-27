function [F_record gamma] = get_attacker_positions_Ver2(N)

%Attacker assignment for dynamics and fixed defenders

%% Initialization

T = 1000;   %time period of evaluation - reducing T => faster evaluation, less accurate for infinite horizon approximation

%System Dynamics matrix - Assign system dynamics to 'A'
%Defender location inputs - Assign defender node inputs to 'B'
A_base = eye(N)+diag(ones(N-1,1),-1)+diag(ones(N-1,1),1);
A      = A_base/(1.1*max(abs(eig(A_base))));
B      = zeros(N,N);

%Cost parameters
Q             = eye(N);
Q_inv         = inv(Q);
R             = eye(N);
R_inv         = inv(R);
gamma_initial = 1000;
gamma_low     = zeros(N,N);
gamma_high    = gamma_initial*ones(N,N);
gamma         = NaN(N,N);
epsilon       = 0.001;                   % evaluation accuracy - adjust to higher values for less accurate convergence
F_record      = zeros(N,N);              % Final Attacker node assignment matrix

%% Evaluation

for player = 1:N
    
     for position = 1:N
         
          F = F_record;
          
          if sum(F(position,:)) == 0
               F(position,player) = 1;
          else
               continue;
          end
          
          check_flag = 1;
          
          while check_flag > 0
              
               Z_inv                  = Q_inv;
               gamma_check_flag       = 1;
               gamma(player,position) = (gamma_low(player,position) + ...
                                        gamma_high(player,position))/2;
                                    
               for t=T:-1:0
                   
                    Z_inv = Q_inv - Q_inv*A'*inv(Z_inv + B*B' - (gamma(player,position)^(-2))*(F*F') + ...
                            A*Q_inv*A')*A*Q_inv;
                       
                    if min(eig(eye(N) + F'*inv(Z_inv*(gamma(player,position)^2) - F*F')*F)) <= 0
                         gamma_check_flag = 0;
                         break;
                    elseif min(eig(R_inv - R_inv*B'*inv(Z_inv + B*R_inv*B')*B*R_inv)) <= 0
                         gamma_check_flag = 0;
                         break;
                    elseif min(abs(eig(Z_inv))) < epsilon
                         gamma_check_flag = 0;
                         break;
                    end
                    
               end
               if gamma_check_flag == 0
                    gamma_low(player,position) = gamma(player,position);
               else
                    gamma_high(player,position) = gamma(player,position);
               end
               if abs(gamma_high(player,position) - gamma_low(player,position)) <= epsilon
                    check_flag = 0;
               end
               gamma;
               
          end
          
     end
     
     gmax = max(gamma(player,:));
     for position = 1:N
          if gmax == gamma(player,position)
               F_record(position,player) = 1;
               break;
          end
     end
     
end


end % end of function

