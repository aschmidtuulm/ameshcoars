function [coordinates,elements,constrains,varargout] = ...
    QcoarsenR(N0,coordinates,elements,constrains,varargin)

%QcoarsenR: local coarsening of quadrilateral mesh  emerged from  
%          QrefineR where marked elements are coarsened by inverting the red
%          refinement
%
%Usage:
%
% [coordinates,elements4,irregular,dirichlet,neumann] ...
%    = QcoarsenR(N0,coordinates,elements4,irregular,dirichlet,neumann,marked)
% or
%
% [coordinates,elements4,irregular] ...
%    = QcoarsenR(N0,coordinates,elements4,irregular,marked)
%
%Comments:
%
%    QcoarsenR expects as input a mesh described by the 
%    fields coordinates, elements4, irregular, dirichlet (optional) and neumann 
%    (optional). N0 specifies the number of coordinates in the initial mesh
%    as the routine does not eliminate nodes from the initial mesh.
%    The vector marked contains the indices of elements which shall be
%    coarsened. A family of elements is coarsened to their common father
%    element when at least one of the elements is marked for coarsening.
%    1-Irregularity of the mesh is ensured by the 1-Irregular Rule.
% 
%    The function returns the coarsened mesh in terms of the same data as
%    for the input.
%
%Remark:
%
%    This program is a supplement to the paper 
%    >> CRITERIA FOR NON-RECURSIVE LOCAL COARSENING IN 2D <<
%    by S. Funken, and A. Schmidt. The reader should 
%    consult that paper for more information.   
%
%Authors:
% 
%    S. Funken, A. Schmidt  02-08-19

nC = size(coordinates,1);
nE = size(elements,1);
marked = varargin{end};
%*** get admissible center nodes
temp = sortrows(elements,3);
adm = temp(temp(1:end-3,3)==temp(4:end,3),3);
adm = adm(adm>N0); % do not delete nodes of initial triangulation
%*** determine elements belonging to admissable center nodes
if isempty(adm)
    varargout = varargin;
    return
end
i2fi = 1:nC; i2fi(adm) = -i2fi(adm);
idx = -min(i2fi(elements),[],2); kdx = find(idx > 0);
temp = sortrows([idx(kdx),kdx]);
node2elements = reshape(temp(:,2),4,[])';
%*** consider only marked elements
marked_elements = zeros(1,nE); marked_elements(marked) = 1;
idx = sum(marked_elements(node2elements),2)>1; % >1 coarsen, if at least one element is marked; ==4 if all elements are marked
node2elements = node2elements(idx,:);
%*** adjust order
i1 = [3,4,4]; i2 = [2,2,3];
for k = 1:3
    idx = find(elements(node2elements(:,i1(k)),4) == ...
        elements(node2elements(:,i2(k)-1),2));
    tmp = node2elements(idx,i2(k));
    node2elements(idx,i2(k)) = node2elements(idx,i1(k));
    node2elements(idx,i1(k)) = tmp;
end
%*** obtain geometric information on edges
[edge2nodes,constrains2edges,element2edges,boundary2edges{1:nargin-5}] ...
    = provideGeometricData(constrains,elements,varargin{1:end-1});
adm2edges = [element2edges(node2elements(:,1),[4,1]), ...
    element2edges(node2elements(:,2),[4,1]), ...
    element2edges(node2elements(:,3),[4,1]), ...
    element2edges(node2elements(:,4),[4,1])];
markedEdges = zeros(size(edge2nodes,1),1);
markedEdges(constrains2edges(:,1)) = 1;
%*** follow 3-neighbor rule
boundinfo = edge2nodes(accumarray([element2edges(:); constrains2edges(:)],1,[length(edge2nodes) 1])==1,:);
idx = find(sum(reshape(markedEdges(adm2edges),[],8),2)); % follow 1-Irregularity
foo = 1;
while ~isempty(idx) || foo
    adm(idx)=[];
    node2elements(idx,:)=[];
    %*** mark nodes to be eliminated
    markedNodes = zeros(nC,1);
    markedNodes(adm) = 1;
    markedNodes = markedNodes+accumarray(elements(node2elements,2),ones(size(elements(node2elements,2))),[nC 1]);
    markedNodes(boundinfo(:,2))= 2*markedNodes(boundinfo(:,2)); % count boundary edges twice 
    markedNodes(constrains(:,3)) = 2*markedNodes(constrains(:,3));% count constraints twice
    idx = find(sum(reshape(markedNodes(elements(node2elements,4)),[],4),2)<6);
    foo = 0;
end
%*** create new elements
idx = setdiff(1:nE,unique(node2elements));
newelements = [elements(idx,:);...
    reshape(elements(node2elements(:),1),[],4)];
%*** update constrains, allow redundancies
constrains = [constrains; ...
    [elements(node2elements(:,1),1),elements(node2elements(:,2),1),elements(node2elements(:,1),2)];...
    [elements(node2elements(:,2),1),elements(node2elements(:,3),1),elements(node2elements(:,2),2)];...
    [elements(node2elements(:,3),1),elements(node2elements(:,4),1),elements(node2elements(:,3),2)];...
    [elements(node2elements(:,4),1),elements(node2elements(:,1),1),elements(node2elements(:,4),2)] ] ;
constrains = sortrows(constrains,3);
idx = find(constrains(1:end-1,3) == constrains(2:end,3));
constrains([idx(:);idx(:)+1],:) = [];
node2bound = zeros(nC,1);
node2bound(unique(boundinfo(:)),:) = 1;
idx = (markedNodes & node2bound);
node2constr = zeros(nC,1);
node2constr(constrains(:,3),1) = 1:size(constrains,1);
constrains(node2constr(idx),:) = [];
%*** update coordinates and element numbering
i2fi = zeros(nC,1);
i2fi(newelements) = 1;
coordinates = coordinates(find(i2fi),:);
i2fi = cumsum(i2fi);
elements = reshape(i2fi(newelements),[],4);
constrains = reshape(i2fi(constrains),[],3);
%*** correct boundaries
for j = 1:nargout-3
    boundary = varargin{j};
    if ~isempty(boundary)
        node2boundary = zeros(nC,2);
        node2boundary(boundary(:,1),1) = 1:size(boundary,1);
        node2boundary(boundary(:,2),2) = 1:size(boundary,1);
        idx = ( markedNodes & node2boundary(:,2) );
        boundary(node2boundary(idx,2),2) = boundary(node2boundary(idx,1),2);
        boundary(node2boundary(idx,1),2) = 0;
        varargout{j} = i2fi(boundary(find(boundary(:,2)),:));
    else
        varargout{j} = [];
    end
end
%*** sort such that oldest node is first
[~,jdx] = min(elements,[],2);
map = [1,2,3,4; 2,3,4,1; 3,4,1,2; 4,1,2,3];
for k=1:4
    idx = jdx == k;
    elements(idx,:) = elements(idx,map(k,:));
end
end

