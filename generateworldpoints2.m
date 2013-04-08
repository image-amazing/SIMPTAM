function [ points ] = generateworldpoints2()
%GENERATEWORLDPOINTS Generates some world points

pointCount = 0;

Corners = [-20 -5 -20; -20 5 -20; 20 -5 -20]';
wall1 = pointcloudplane(Corners, 10);

for i = 1:size(wall1,2)
    pointCount = pointCount + 1;
    points(pointCount).id = pointCount;
    points(pointCount).location = wall1(:,i);
end

Corners = [20 -5 -20; 20 5 -20; 20 -5 20]';
wall1 = pointcloudplane(Corners, 10);

for i = 1:size(wall1,2)
    pointCount = pointCount + 1;
    points(pointCount).id = pointCount;
    points(pointCount).location = wall1(:,i);
end

Corners = [-20 -5 20; -20 5 20; 20 -5 20]';
wall1 = pointcloudplane(Corners, 10);

for i = 1:size(wall1,2)
    pointCount = pointCount + 1;
    points(pointCount).id = pointCount;
    points(pointCount).location = wall1(:,i);
end

Corners = [-20 -5 -20; -20 5 -20; -20 -5 20]';
wall1 = pointcloudplane(Corners, 10);

for i = 1:size(wall1,2)
    pointCount = pointCount + 1;
    points(pointCount).id = pointCount;
    points(pointCount).location = wall1(:,i);
end








end

