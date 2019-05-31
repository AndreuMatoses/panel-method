function[field]=fields_NL(panels,solution,tamany,N)
%tamany=[times x thickness, times y thickness], N=n of points for
%discretization

interval_x=(max(panels.vertex(:,1))-min(panels.vertex(:,1)));
interval_y=(max(panels.vertex(:,2))-min(panels.vertex(:,2)));

max_x=max(panels.vertex(:,1))+(tamany(1)/2-0.5)*interval_x;
min_x=min(panels.vertex(:,1))-(tamany(1)/2-0.5)*interval_x;
max_y=max(panels.vertex(:,2))+(tamany(2)/2-0.5)*interval_y;
min_y=min(panels.vertex(:,2))-(tamany(2)/2-0.5)*interval_y;

field.x=repmat(linspace(min_x,max_x,N),N,1);
field.y=repmat(linspace(min_y,max_y,N)',1,N);
field.u=zeros(N);
field.v=zeros(N);
% field.x=linspace(min_x,max_x,N);
% field.y=linspace(min_y,max_y,N)';

for i=1:N
    for j=1:N
        xci=field.x(i,j);
        yci=field.y(i,j);
        xbj=panels.vertex(:,1); %vector
        ybj=panels.vertex(:,2); %vector
        phij=panels.phi;
        Sj=panels.long;
        A=-(xci-xbj).*cosd(phij)-(yci-ybj).*sind(phij);
        B=(xci-xbj).^2+(yci-ybj).^2;
        C=sind(-phij);
        D=cosd(-phij);
        F=log(1+(Sj.^2+2.*A.*Sj)./(B));
        E=(xci-xbj).*sind(phij)-(yci-ybj).*cosd(phij);
        G=atan((E.*Sj)./(A.*Sj+B));
        
        % check if points are insed the profile and make thme v=0
        if inpolygon(field.x(i,j),field.y(i,j),panels.vertex(:,1),panels.vertex(:,2))
            field.u(i,j)=0;
            field.v(i,j)=0;
        else
            field.u(i,j)=(solution.V_inf*cosd(-solution.alpha)-sum((solution.lambda./(2*pi)).*((D.*F)./2+C.*G)));
            field.v(i,j)=-solution.V_inf*sind(-solution.alpha)-sum((solution.lambda./(2*pi)).*(-(C.*F)./2+D.*G));
        end
    end
end
field.V=sqrt(field.u.^2+field.v.^2);
field.Cp=1-(field.V./solution.V_inf).^2;
field.p=field.Cp*0.5*1.225*solution.V_inf^2+101325; %atmos conditions
end

