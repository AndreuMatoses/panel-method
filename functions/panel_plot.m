function[figura]=panel_plot(panels)

figura=plot([panels.vertex(:,1);panels.vertex(1,1)],[panels.vertex(:,2);panels.vertex(1,2)],panels.vertex(:,1),panels.vertex(:,2),'o',panels.mid_points(:,1),panels.mid_points(:,2),'*');
text(panels.mid_points(:,1),panels.mid_points(:,2),num2str([1:length(panels.mid_points(:,2))]'))
ax=gca;
ax.YGrid = 'on';
pbaspect([1 1 1]);
axis normal
% xlim([-25,25])
% ylim([-25,25])
title('Panels')
legend('Linear Panels','Vertex','Mid-Points','Location','best')
end