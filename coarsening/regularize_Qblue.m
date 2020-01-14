function [coordinates, newElements, varargout] = regularize_Qblue(coordinates, elements, irregular, varargin)
global nB
nE = size(elements,1);
%*** Obtain geometric information on edges
[edge2nodes,irregular2edges,element2edges,boundary2edges{1:nargin-3}]...
    = provideGeometricData(irregular,elements,varargin{1:end});
%*** Mark edges for refinement
hash = [1,1,1,1;1,1,0,0;0,0,1,1];
[map,value] = hash2map((0:15)',hash);
edge2newNode = zeros(1,size(edge2nodes,1));
edge2newNode(irregular2edges(:,1)) = 1;
swap = 1;
while ~isempty(swap)
    markedEdges = edge2newNode(element2edges);
    dec = sum(markedEdges.*(ones(nE,1)*2.^(0:3)),2);
    val = value(dec+1);
    [idx,jdx] = find(~markedEdges & map(dec+1,:));
    swap = idx+(jdx-1)*nE;
    edge2newNode(element2edges(swap)) = 1;
end
%*** Generate new nodes on edges
edge2newNode(irregular2edges(:,1)) = 0;
edge2newNode(edge2newNode~=0) = size(coordinates,1)...
                                    + (1:nnz(edge2newNode));
idx = find(edge2newNode);
coordinates(edge2newNode(idx),:) = (coordinates(edge2nodes(idx,1),:)...
                                + coordinates(edge2nodes(idx,2),:))/2;
edge2newNode(irregular2edges(:,1)) = irregular(:,3);                            
%*** Refine boundary conditions
varargout = cell(nargout-2,1);
for j = 1:nargout-2
    boundary = varargin{j};
    if ~isempty(boundary)
        newNodes = edge2newNode(boundary2edges{j})';
        markedEdges = find(newNodes);
        if ~isempty(markedEdges)
            boundary = [boundary(~newNodes,:); ...
                boundary(markedEdges,1),newNodes(markedEdges); ...
                newNodes(markedEdges),boundary(markedEdges,2)];
        end
    end
    varargout{j} = boundary;
end
%*** Provide new nodes for refinement of elements
newNodes = edge2newNode(element2edges);
%*** Determine type of refinement for each red element
none   = find(val == 0);
red    = find(val == 1);
bluer = find(val == 2);
bluel = find(val == 3);
%*** Generate new interior nodes if red elements are refined
idx = [red,bluer,bluel];
midNodes = zeros(nE,1);
midNodes(idx) = size(coordinates,1)+(1:length(idx));
coordinates = [coordinates; ...
    ( coordinates(elements(idx,1),:) ...
    + coordinates(elements(idx,2),:) ...
    + coordinates(elements(idx,3),:) ...
    + coordinates(elements(idx,4),:) )/4];
%*** Generate element numbering for refined mesh
rdx = zeros(nE,1);
rdx(none)    = 1;
rdx(red)     = 4;
rdx = [1;1+cumsum(rdx)];
bdx = zeros(nE,1);
bdx([bluer,bluel]) = 3;
bdx = rdx(end)+[0;0+cumsum(bdx)];
%*** Generate new red elements
tmp = [elements,midNodes,newNodes];
%*** Generate new red elements first
newElements = 1+zeros(bdx(end)-1,4);
newElements(rdx(none),:) = elements(none,:);
newElements([rdx(red),1+rdx(red),2+rdx(red),3+rdx(red)],:) ...
    = [tmp(red,[1,6,5,9]);tmp(red,[2,7,5,6]);...
    tmp(red,[3,8,5,7]);tmp(red,[4,9,5,8]);];
%*** New blue elements
newElements([bdx(bluer),1+bdx(bluer),2+bdx(bluer)],:) ...
    = [tmp(bluer,[3,4,5,7]);tmp(bluer,[1,6,5,4]);tmp(bluer,[6,2,7,5])];
newElements([bdx(bluel),1+bdx(bluel),2+bdx(bluel)],:) ...
    = [tmp(bluel,[1,2,5,9]);tmp(bluel,[3,8,5,2]);tmp(bluel,[8,4,9,5])];
nB = size(newElements,1)-rdx(end)+1;
end