function [coordinates, newElements] = regularize_Tgreen(coordinates, elements, irregular)
global nG
nE = size(elements,1);
% obtain geometric information
[~,element2edges,~] = provideGeometricData([elements;irregular],zeros(0,4));
irregular2edges = element2edges(nE+1:end,:);
element2edges = element2edges(1:nE,:);
edge2elements = createEdge2Elements(element2edges);
gdx = edge2elements(irregular2edges(:,1));
%*** sort all to the same position
for k =1:2
jdx = find(elements(gdx)~=irregular(:,1));
elements(gdx(jdx),:) = elements(gdx(jdx),[2,3,1]);
end
green = nE+1:2:(nE+2*size(gdx,1));
if isempty(green)
    newElements  = zeros(nE,3);
else
    newElements = zeros(green(end)+1,3);
end
%*** define new elements
newElements(1:nE,:) = elements;
newElements([green,green+1],:) = [elements(gdx,[3,1]),irregular(:,3); irregular(:,3),elements(gdx,[2,3])];
newElements(gdx,:) = [];
nG = 2*size(gdx,1);
end

