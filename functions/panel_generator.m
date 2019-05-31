function[panels]=panel_generator(points,n_panels,slope_sensitive)
%points=two x,y columns vector, sorted
%n_panels=number of panels, int
%slope_sensitive=from 0 to 1, how sensible is panel length to slope of curve
[~,Imin]=min(points(:,1));
lower_points=points(1:Imin-1,:); %points before LE
upper_points=points(Imin:end,:); %points after LE + LE 
%plot(upper_points(:,1),upper_points(:,2),lower_points(:,1),lower_points(:,2))

%delete redundant points
[~,no_rep_up]=unique(upper_points(:,1),'stable');
upper_points=upper_points(no_rep_up,:);

[~,no_rep_down]=unique(lower_points(:,1),'stable');
lower_points=lower_points(no_rep_down,:);

%% a lot of random points to calculate approximated distances and the derivative of a 'line' going throug the points upper/lower
discretization=500; % times n of panels
random_x=linspace(min(points(:,1)),max(points(:,1)),n_panels*discretization+1)';
lower_line=[random_x(end:-1:2) interp1(lower_points(:,1),lower_points(:,2),random_x(end:-1:2),'spline','extrap')];
upper_line=[random_x(1:end-1) interp1(upper_points(:,1),upper_points(:,2),random_x(1:end-1),'spline','extrap')]; %basically upper_points but with a lot of points


whole_line=[lower_line;upper_line;lower_line(1,:)]; %repeat initail point at the end to close the cirlce
%plot(whole_line(:,1),whole_line(:,2),upper_points(:,1),upper_points(:,2),'o',lower_points(:,1),lower_points(:,2),'o')

distances_vector=whole_line(2:end,:)-whole_line(1:end-1,:);
distances_abs=sqrt(distances_vector(:,1).^2+distances_vector(:,2).^2);
total_distance=sum(distances_abs); %add the last point to the first point distance
%for the circle, total_distance should be perimeter
fprintf('Error of total ditance suposing circle of r1: %g %%, perimeter: %g \n',((total_distance-2*pi*1)/(2*pi*1)*100),total_distance)

derivative=distances_vector(:,2)./distances_vector(:,1); %For the sphere should be like the plot of the tangent from 0 to 2pi
%plot(derivative)

%% Using the abs of derivative as a density function for the lenght fo the panels:
derivative=abs(derivative);
% Group derivative in n_panels groups
cumulative_derivative=ones(n_panels,1);
for i=1:n_panels
    cumulative_derivative(i,1)=sum(derivative(discretization*2*(i-1)+1:discretization*2*i,1));
end
%plot(cumulative_derivative)
%applaying the sensitivity to density function parameter, to flatten density, approaching it to the mean 
cumulative_derivative=cumulative_derivative+(mean(cumulative_derivative)-cumulative_derivative)*(1-slope_sensitive);
%normalize the function so sum=1, which means that the sum(derivative) of 1
%means a panel lenght of total_distance
cumulative_derivative=1./cumulative_derivative;
cumulative_derivative=cumulative_derivative/sum(cumulative_derivative);
% figure(2)
% plot(cumulative_derivative)
% ylim([0 0.25])
panels_long_ideal=total_distance*cumulative_derivative;
%by sight, as the curve is always longer than the lines
%% find points that accomplish best the desired distance

pan_index=1;
correction=1; %by sight, as the curve is always longer than the lines, ex:0.99
panels.vertex(pan_index,:)=whole_line(1,:);
for k=1:length(whole_line)
    if pdist([whole_line(k,:);panels.vertex(pan_index,:)])>=panels_long_ideal(pan_index)*correction
       pan_index=pan_index+1;
       panels.vertex(pan_index,:)=whole_line(k,:);
    end
end
%check if the last vertex is the same as the first
if panels.vertex(end,:)==panels.vertex(1,:)
   panels.vertex=panels.vertex(1:end-1,:);%delete last vertex
   pan_index=pan_index-1;
end
fprintf('It has been possible to make %g panels\n',pan_index)

%% midpoints
panels.mid_points=[(panels.vertex(2:end,:)+panels.vertex(1:end-1,:))/2;(panels.vertex(end,:)+panels.vertex(1,:))/2];

%% real panel parameters
temp_vertex=[panels.vertex;panels.vertex(1,:)];
distances_vector=temp_vertex(2:end,:)-temp_vertex(1:end-1,:);
panels.vectors=distances_vector;
panels.long=sqrt(distances_vector(:,1).^2+distances_vector(:,2).^2);%;pdist([panels.vertex(1,:);panels.vertex(end,:)]);
panels.phi=atan2d(panels.vectors(:,2),panels.vectors(:,1));
panels.phi(panels.phi<0)=360+panels.phi(panels.phi<0);
% panels.beta=panels.phi+90;
% panels.beta(panels.beta>360)=panels.beta(panels.beta>360)-360;
% panels.normal=[cosd(panels.beta) sind(panels.beta)];
panels.TE_vertex=1;
[~,panels.LE_vertex]=min(panels.vertex(:,1));
[~,panels.Top_vertex]=max(panels.vertex(:,2));
[~,panels.Bot_vertex]=min(panels.vertex(:,2));
panels.lowerindexes=panels.TE_vertex:(panels.LE_vertex-1);
panels.upperindexes=panels.LE_vertex:length(panels.mid_points);
panels.rightindexes=[1:panels.Bot_vertex-1 panels.Top_vertex:length(panels.mid_points)]';
panels.leftindexes=[panels.Bot_vertex:panels.Top_vertex-1]';
panels.upperx=(panels.mid_points(panels.upperindexes,1)-panels.vertex(panels.LE_vertex))/(panels.vertex(panels.TE_vertex)-panels.vertex(panels.LE_vertex));
panels.lowerx=(panels.mid_points(panels.lowerindexes,1)-panels.vertex(panels.LE_vertex))/(panels.vertex(panels.TE_vertex)-panels.vertex(panels.LE_vertex));
end