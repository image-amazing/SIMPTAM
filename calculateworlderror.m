function [ error ] = calculateworlderror(World, EstWorld )
%CALCULATEWORLDERROR Summary of this function goes here
%   Detailed explanation goes here

error = 0;

for i = 1:length(EstWorld.points)
    pointid = EstWorld.points(i).id;
    if (World.points(pointid).id ~= pointid)
        display('There seems to be a problem');
    else
        error = error + norm(World.points(pointid).location-EstWorld.points(i).location);
%         display('=======');
%         display(World.points(pointid).location);
%         display(EstWorld.points(i).location);
    end
end



end

