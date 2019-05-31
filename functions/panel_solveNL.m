function[solution]=panel_solveNL(panels,V_inf,alpha)

%% Conditions of the solution
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
%% linear sistem and solution
Dn= -((C.*F)./2)+D.*G-eye(N)*pi;
Dn=Dn./(2*pi);

R=V_inf*sind(alpha-panels.phi(:,1));

q=linsolve(Dn,R);
solution.lambda=q(1:N);

solution.lambda_balance=sum(solution.lambda.*panels.long);

%% Cp annd Cl calculation
Vt=zeros(N,1);
for i=1:N
    
    Vs=zeros(N,1);
    for j=1:N
        Vs(j)=solution.lambda(j)/(2*pi)*(D(i,j)*F(i,j)/2+C(i,j)*G(i,j));
    end
    Vt(i)=V_inf*cosd(panels.phi(i)-alpha)-sum(Vs);
end
solution.V=Vt;

solution.Cp=1-(solution.V./V_inf).^2;
solution.Cpupper=-solution.Cp(panels.upperindexes);
solution.Cplower=-solution.Cp(panels.lowerindexes);

% solution.cx=trapz(panels.mid_points(:,1),solution.Cp);
% solution.cy=trapz(panels.mid_points(:,2),solution.Cp);
% solution.cl_I=-cx*sind(alpha)-cy*cosd(alpha);
% solution.cd_I=cx*cosd(alpha)+cy*sind(alpha);
end