function [ points ] = generateworldpoints3(r,h)
%GENERATEWORLDPOINTS Generates some world points

points = [];

for i = 1:400
    x = rand*2*r - r;
    y = rand*2*h - h;
    z = (2*randi([0 1]) - 1)*sqrt(r^2 - x^2);
    
    points(i).id = i;
    points(i).location = [x y z 1]';
end

    







end

