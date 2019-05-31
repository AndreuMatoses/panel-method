function[]=plot_solution(solution,panels)

% subplot(1,2,1)
[theta,rho] = cart2pol(panels.mid_points(:,1),panels.mid_points(:,2));
theta(theta<=0)=2*pi+theta(theta<=0);
theta=rad2deg(theta);

theta_ideal=[0:45:360];
Cp_ideal=[1 -1 -3 -1 1 -1 -3 -1 1];

leg={};
% subplot(2,1,1)
for i=1:length(solution)
plot(theta,solution(i).Cp,'-o',theta_ideal,Cp_ideal,'^r','LineWidth',1.5)
hold on
leg{i}=['Cp for alpha= ' num2str(solution(i).alpha)];
end
leg{1}='Cp numerical result';
leg{2}='Analytical value of Cp';
title('Cp as a function of polar angle')
xlabel('theta (deg)')
ylabel('Cp')
xlim([0 360])
legend(leg)
hold off

% subplot(2,1,2)
% leg={};
% for i=1:length(solution)
% plot(panels.mid_points(panels.upperindexes,1),-solution(i).Cp(panels.upperindexes),'-^',panels.mid_points(panels.lowerindexes,1),-solution(i).Cp(panels.lowerindexes),'-^')
% hold on
% leg{2*i-1}=['Alpha= ' num2str(solution(i).alpha) 'upper'];
% leg{2*i}=['Alpha= ' num2str(solution(i).alpha) ' lower'];
% end
% title('Cp as a function of x/c')
% xlabel('x/c')
% legend(leg)
% hold off

end