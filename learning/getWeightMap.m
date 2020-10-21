 function [weightMap] = getWeightMap(edge, mask)
     %maxDist = sqrt(size(edge,1)^2 + size(edge,1)^2);
     maxDist = sqrt(size(edge,1)^1.2 + size(edge,1)^1.2);
     Bdi = bwdist(edge);
     Bdi = (maxDist - Bdi);
     Bdi = Bdi./maxDist;
     weightMap = mask.*Bdi;

end
