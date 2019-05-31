function[section]=angular_sorting(section)
x = section.points(:,1);
y = section.points(:,2);
z = section.points(:,3);
xCenter = mean(x);
yCenter = mean(y);
zCenter = mean(z);
if (find(section.normal)==2)
    [theta,~] = cart2pol((x-xCenter),(z-zCenter));
    [~,Imax]=max(x);
    [theta0,~]=cart2pol((x(Imax)-xCenter),(z(Imax)-zCenter));
elseif (find(section.normal)==3)
    [theta,~] = cart2pol((x-xCenter),(y-yCenter));
    [~,Imax]=max(x);
    [theta0,~]=cart2pol((x(Imax)-xCenter),(y(Imax)-yCenter));
else
    [theta,~] = cart2pol((y-yCenter),(z-zCenter));
    [~,Imax]=max(y);
    [theta0,~]=cart2pol((y(Imax)-yCenter),(z(Imax)-zCenter));
end
%initial sortin angle:
theta=theta-theta0;
theta(theta<=0)=theta(theta<=0)+2*pi;
theta(theta<=0)=theta(theta<=0)+2*pi;

[sorttheta, sortedIndexes] = sort(theta,'descend');
% lower=sorttheta>mod(-pi,360);
% upper=sorttheta<=mod(-pi,360);
section.points = [x(sortedIndexes),y(sortedIndexes),z(sortedIndexes)];
% section.points=[section.points(upper,:);section.points(lower,:)];
end