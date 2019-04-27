clear;
clc;
%Attacker assignment for dynamics and fixed defenders

%% Initialization

n=10;     %number of states
T=1000;   %time period of evaluation - reducing T => faster evaluation, less accurate for infinite horizon approximation

%System Dynamics matrix - Assign system dynamics to 'A'
A_base=eye(n)+diag(ones(n-1,1),-1)+diag(ones(n-1,1),1);
A=A_base/(1.1*max(abs(eig(A_base))));
% A_base(1,n)=1;A_base(n,1)=1;
% eig(A_base)'

%Defender location inputs - Assign defender node inputs to 'B'
B = zeros(n,n);

%Cost parameters
Q=eye(n);
Qinv = inv(Q);
R=eye(n);
Rinv = inv(R);

gl = zeros(n,n);
gu = 1000*ones(n,n);
g = NaN(n,n);

acc = 0.001;   %evaluation accuracy - adjust to higher values for less accurate convergence

Frec = zeros(n,n);  %Final Attacker node assignment matrix

%% Evaluation

for i=1:n
     F_i = Frec(:,:);
     
     for j=1:n
          F=F_i;
          if sum(F(j,:))==0
               F(j,i)=1;
          else
               continue;
          end
          check = 1;
          while check>0
               Zinv = Qinv;
               check2=1;
               g(i,j)=(gl(i,j)+gu(i,j))/2;
               for t=T:-1:0
                    Zinv = Qinv - Qinv*A'*inv(Zinv + B*B' - (g(i,j)^(-2))*(F*F') + A*Qinv*A')*A*Qinv;
                    if min(eig(eye(n) + F'*inv(Zinv*(g(i,j)^2) - F*F')*F))<=0
                         check2=0;
                         break;
                    elseif min(eig(Rinv - Rinv*B'*inv(Zinv + B*Rinv*B')*B*Rinv)) <=0
                         check2=0;
                         break;
                    elseif min(abs(eig(Zinv)))<acc
                         check2=0;
                         break;
                    end
               end
               if check2==0
                    gl(i,j) = g(i,j);
               else
                    gu(i,j)=g(i,j);
               end
               if abs(gu(i,j)-gl(i,j))<=acc
                    check=0;
               end
               g(:,:)
          end
     end
     gmax = max(g(i,:));
     for j = 1:n
          if gmax == g(i,j)
               Frec(j,i)=1;
               break;
          end
     end
end

Frec