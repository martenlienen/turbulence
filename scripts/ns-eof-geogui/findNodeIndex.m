function [ index ] = findNodeIndex( indexSet, nodeList )
%FINDNODEINDEX Summary of this function goes here
%   Detailed explanation goes here

for i = 1:1:size(nodeList,1)
    if nodeList(i,3)==indexSet(1) && nodeList(i,4)==indexSet(2)
        index = i;
        break;
    end
end

end

