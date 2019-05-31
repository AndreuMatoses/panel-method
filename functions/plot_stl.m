function[]=plot_stl(model,section)

patch(model,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15,'FaceAlpha',0.5);
     
% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');
axis('image');
view([45 35]);
hold on;
%plot the section plane
w = null(section.normal); % Find two orthonormal vectors which are orthogonal to v
max_points=max(section.points);
min_points=min(section.points);
max_points=max_points(section.normal==0);
min_points=min_points(section.normal==0);
[P,Q] = meshgrid(-[max_points(1),min_points(1)],[max_points(2),min_points(2)]); % Provide a gridwork (you choose the size)
X = section.origin(1)+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
Y = section.origin(2)+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
Z = section.origin(3)+w(3,1)*P+w(3,2)*Q;
surf(X,Y,Z,'FaceAlpha',0.5,'FaceColor','r')
%plot intersection points
plot3(section.points(:,1),section.points(:,2),section.points(:,3),'*')
hold off;
end