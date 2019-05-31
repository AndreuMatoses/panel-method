function[solution]=panel_solverL(panels,V_inf,alpha)
%% Conditions of the solution
gamma=-50; %initial guess

solution.alpha=alpha;
solution.V_inf=V_inf;
solution.NofPanels=length(panels.mid_points);

%% coeficients for the integrals w.r.t s and n
N=length(panels.mid_points);
A=zeros(N);
B=zeros(N);
C=zeros(N);
D=zeros(N);
E=zeros(N);
F=zeros(N);
G=zeros(N);

for i=1:N
    for j=1:N
        xci=panels.mid_points(i,1);
        yci=panels.mid_points(i,2);
        xbj=panels.vertex(j,1);
        ybj=panels.vertex(j,2);
        phij=panels.phi(j);
        phii=panels.phi(i);
        Sj=panels.long(j);
        
        A(i,j)=-(xci-xbj)*cosd(phij)-(yci-ybj)*sind(phij);
        B(i,j)=(xci-xbj)^2+(yci-ybj)^2;
        C(i,j)=sind(phii-phij);
        D(i,j)=cosd(phii-phij);
        F(i,j)=log(1+(Sj^2+2*A(i,j)*Sj)/(B(i,j)));
        E(i,j)=(xci-xbj)*sind(phij)-(yci-ybj)*cosd(phij);
        G(i,j)=atan((E(i,j)*Sj)/(A(i,j)*Sj+B(i,j)));
        
    end
end

%% linear sistem coeff in two steps where R is calculated with assumed gamma
Dn= -((C.*F)./2)+D.*G-eye(N)*pi;
Dn=Dn./(2*pi);

Qn=((D.*F)./2)+C.*G;
Qn=Qn./(2*pi);

%% Solving the LS for a converging gamma untill kutta condition is met
while (1)
    R=V_inf*sind(alpha-panels.phi(:,1))+Qn*gamma*ones(N,1);
    
    q=linsolve(Dn,R);
    
    solution.lambda=q(1:N);
    solution.gamma=gamma;
    
    solution.lambda_balance=sum(solution.lambda.*panels.long);
    
    %% Vt calculation for Cp
    Vt=zeros(N,1);
    for i=1:N
        Vs=zeros(N,1);
        Vg=zeros(N,1);
        for j=1:N
            Vs(j)=solution.lambda(j)/(2*pi)*(D(i,j)*F(i,j)/2+C(i,j)*G(i,j));
            if i==j
                Vg(j)=pi; %%%
            else
                Vg(j)=C(i,j)*F(i,j)/2-D(i,j)*G(i,j);
            end
        end
        Vt(i)=V_inf*cosd(panels.phi(i)-alpha)-sum(Vs)+solution.gamma/(2*pi)*sum(Vg);
    end
    error=abs(Vt(1)+Vt(N));
    relaxation=0.1;
    if error<(0.001)
        break;
    else
        gamma=gamma+error*relaxation;
    end
end
solution.gamma=gamma;
solution.V=Vt;

%% Cp and Cl calculation
solution.Cp=1-(solution.V./V_inf).^2;
solution.Cpupper=-solution.Cp(panels.upperindexes);
solution.Cplower=-solution.Cp(panels.lowerindexes);

solution.cl_kutta=2*solution.gamma/V_inf/abs((panels.vertex(panels.LE_vertex,1)-panels.vertex(panels.TE_vertex,1)))*sum(panels.long);

% solution.cx=trapz(panels.mid_points(:,1),solution.Cp);
% solution.cy=trapz(panels.mid_points(:,2),solution.Cp);
% solution.cl_I=-cx*sind(alpha)-cy*cosd(alpha);
% solution.cd_I=cx*cosd(alpha)+cy*sind(alpha);
end