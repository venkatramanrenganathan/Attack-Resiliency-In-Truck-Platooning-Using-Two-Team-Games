%%Defenders assumed at all nodes
%%Attackers assigned by greedy algorithm to location with highest crit
%%gamma corresponding to best position for attack

n=10;       %number of states

T=100;        %time period of evaluation. Set to ~10000 for infinite horizon approximation
Q=eye(n);       %state cost
R=eye(n);       %input cost

%System Dynamics Matrix - line graph
A_base=eye(n)+diag(ones(n-1,1),-1)+diag(ones(n-1,1),1);
A=A_base/(0.1+max(abs(eig(A_base))));

%Input Matrix - assumes input to all nodes
B = ones(n);

%Calculation initialization
gl=zeros(n);
ginit = 100;
gu=ginit*ones(n);
g=zeros(n);
pos = zeros(n,1);   %final order of attack
grec = NaN(n,1);
F_base = zeros(n);

for i = 1:n
    for j = 1:n
        F = F_base;
        if sum(F(j,:))~=0
            continue
        end
        F(j,i)=1;
        check=1;
        while check>0
            g(j,i) = (gl(j,i)+gu(j,i))/2;
            Z=Q;
            check2=0;
            check3=0;
            for t=1:T
                Z=Q+A'*Z*inv(eye(n)+ (B*B' - (g(j,i)^-2)*F*F')*Z)*A;
                if min(eig(eye(n)*g(j,i)^2 - F'*Z*F)) <=0
                    check2=1;
                    break
                end
                if min(eig(R + B'*Z*B)) <=0
                    check3=1;
                    break
                end
            end
            if isnan(trace(Z))==1 || check2==1 || check3==1
                gl(j,i)=g(j,i);
            else
                if abs(gu(j,i)-g(j,i))<0.0001
                    check=0;
                end
                if abs(gl(j,i)-ginit)<0.0001
                    check=-1;
                end
                gu(j,i)=g(j,i);
            end
        end
        if pos(i)==0
            pos(i)=j;
            grec(i)=g(j,i);
            dummy=1;
        elseif grec(i)<g(j,i) && check==0
            pos(i)=j;
            grec(i)=g(j,i);
            dummy=1;            
        end
    end
    F_base(pos(i),i)=1;
    g
end

figure
imagesc(g(:,:))
colorbar()
xlabel('Attacker number')
ylabel('Attacker position')
