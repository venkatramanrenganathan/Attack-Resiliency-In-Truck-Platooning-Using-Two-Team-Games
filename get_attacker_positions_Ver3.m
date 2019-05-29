function [F_record,g] = get_attacker_positions_Ver3(N,N_att)

%output function is the strongest attacker position matrix F for discrete time and gamma for evaluation

%% Initialization

T = 1000;   %time period of evaluation - reducing T => faster evaluation, less accurate for infinite horizon approximation

%System Dynamics matrix - Assign system dynamics to 'A'
%Defender location inputs - Assign defender node inputs to 'B'
A_base = eye(N)+diag(ones(N-1,1),-1)+diag(ones(N-1,1),1);
A = 0.99*A_base/max(abs(eig(A_base)));
B_record = zeros(N,N);

%Cost parameters
Q             = eye(N);
R             = eye(N);
g_high = 1000;
g_low     = 0;
% evaluation accuracy - adjust to higher values for less accurate convergence
g_conv_limit = 10^-5;
P_conv_limit = 10^-6;
T_limit = 10000;

%% Evaluation

%Value of gamma
F_record      = zeros(N,N_att);
gamma_record = NaN(N,N_att);
for i = 1:N_att
     %      parfor j = 1:N %If using parallel processing
     for j = 1:N         %If not using parallel processing
          F = F_record;
          B = B_record;
          if sum(F(j,:)) ~= 0
               continue
          else
               F(j,i) = 1;
               gl = g_low;
               gu = g_high;
               g_converge = 0;
               while g_converge == 0
                    g = (gl+gu)/2;
                    P = Q;
                    Phis = zeros(size(P));
                    P_converge = 0;
                    t = 0;
                    while P_converge == 0
                         t = t+1;
                         mu = 
                         alpha = inv(B'*P*B + R);
                         beta = inv(F'*P*F - gamma^2*eye(N_att)
                         
                         
                         L_mat = g^2*eye(N_att) - F'*P*F;
                         M_mat = R + B'*P*B;
                         if min(eig(L_mat)) < 0 || min(eig(M_mat)) < 0
                              P_converge = 2;
                              break;
                         end
                         P = Q + A'*P*inv(eye(N) + (B*B' - (g^(-2))*(F*F'))*P)*A;
                         if t > T_limit
                              P_converge = 2;
                              break;
                         end
                         if norm(P-Phis,'fro') <= P_conv_limit
                              P_converge = 1;
                         else
                              Phis = P;
                         end
                    end
                    if P_converge == 2
                         gl = g;
                    elseif P_converge == 1
                         gu = g;
                         gamma_record(j,i) = g;
                    end
                    if abs(gu-gl) < g_conv_limit
                         g_converge = 1;
                    end
               end
          end
     end
     [~,pos] = max(gamma_record(:,i));
     F_record(pos,i) = 1;
     
end

g = max(max(gamma_record));

%attacker position wrt cost for fixed gamma
cost_record = zeros(N,N_att);
F_record = zeros(N,N_att);

for i = 1:N_att
     %      parfor j = 1:N %If using parallel processing
     for j = 1:N         %If not using parallel processing
          F = F_record;
          B = B_record;
          if sum(F(j,:)) ~= 0
               continue
          else
               F(j,i) = 1;
               P = Q;
               Phis = zeros(size(P));
               P_converge = 0;
               t = 0;
               while P_converge == 0
                    t = t+1;
                    L_mat = g^2*eye(N_att) - F'*P*F;
                    M_mat = R + B'*P*B;
                    if min(eig(L_mat)) < 0 || min(eig(M_mat)) < 0
                         P_converge = 2;
                         break;
                    end
                    P = Q + A'*P*inv(eye(N) + (B*B' - (g^(-2))*(F*F'))*P)*A;
                    if t > T_limit
                         P_converge = 2;
                         break;
                    end
                    if norm(P-Phis,'fro') <= P_conv_limit
                         P_converge = 1;
                    else
                         Phis = P;
                    end
               end
               if P_converge == 1
                    cost_record(j,i) = trace(Phis);
               end
          end
     end
     [~,pos] = max(cost_record(:,i));
     F_record(pos,i) = 1;
end
end