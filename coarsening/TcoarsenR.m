function [coordinates,elements,irregular,varargout] = ...
    TcoarsenR(N0,coordinates,elements,irregular,varargin)

%TcoarsenR: local coarsening of triangular mesh emerged from  
%          TrefineR where marked elements are coarsened by inverting the red
%          refinement
%
%Usage:
%
% [coordinates,elements3,irregular,dirichlet,neumann] ...
%    = TcoarsenR(N0,coordinates,elements3,irregular,dirichlet,neumann,marked)
% or
%
% [coordinates,elements3,irregular] ...
%    = TcoarsenR(N0,coordinates,elements3,irregular,marked)
%
%Comments:
%
%    TcoarsenR expects as input a mesh described by the 
%    fields coordinates, elements3, irregular, dirichlet (optional) and neumann 
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
%*** sort such that oldest node is first
[~,jdx] = min(elements,[],2);
map = [1,2,3;2,3,1;3,1,2];
for k=1:3
    idx = jdx ==k;
    elements(idx,:) = elements(idx,map(k,:));
end
%*** obtain geometric information
[edge2nodes,element3edges,~] = provideGeometricData([elements;...
    irregular],zeros(0,4),varargin{1:end-1});
irregular2edges = element3edges(nE+1:end,:);
element3edges = element3edges(1:nE,:);
edge2element = createEdge2Elements(element3edges);
%*** determine elements to coarsen
minNode = elements(:,1);
value2edges = minNode.*ones(nE,3);
minEdgeVal = accumarray(element3edges(:),value2edges(:),[size(edge2nodes,1) 1],@min);
midElements = find(~((minEdgeVal(element3edges(:,1))==minEdgeVal(element3edges(:,2))).*...
    (minEdgeVal(element3edges(:,2))==minEdgeVal(element3edges(:,3)))));
idx = sum(elements(midElements,:)>N0,2)<3; % do not delete nodes from initial triangulation
midElements(idx) = [];
midEdges = element3edges(midElements,:);
elem2coarse = reshape(edge2element(midEdges(:),:),[],6);
idx = any(elem2coarse==0,2);
elem2coarse(idx,:)= []; midElements(idx) = [];
elem2coarse_tmp = elem2coarse';
idx = (reshape(midElements',1,[]) .* ones(6,size(midElements,1)));
elem2coarse = [reshape(elem2coarse_tmp(elem2coarse_tmp~=idx),3,[])', midElements];
%*** adjust order of elements
idx = find(elements(elem2coarse(:,1),2) == ...
    elements(elem2coarse(:,3),3));
tmp = elem2coarse(idx,3);
elem2coarse(idx,3) = elem2coarse(idx,2);
elem2coarse(idx,2) = tmp;
%*** consider unfavorable numbering in initial mesh
jdx = elem2coarse(:,1:3);
dummy = minEdgeVal(element3edges(jdx(:),:));
idx = find(~((dummy(:,1)==dummy(:,2)).* (dummy(:,2)==dummy(:,3)).* (dummy(:,1)==dummy(:,3))));
row_idx = mod(idx,size(elem2coarse,1));
row_idx(row_idx==0) = size(elem2coarse,1);
elem2coarse(row_idx,:)= [];
%*** consider only marked elements
marked_elements = zeros(1,nE); marked_elements(marked) = 1;
idx = sum(marked_elements(elem2coarse),2)>1; % >1 coarsen, if at least one element is marked; ==4 if all elements are marked
elem2coarse = elem2coarse(idx,:);
%*** follow 1-irregular rule
adm2edges = [element3edges(elem2coarse(:,1),[1,3]), ...
    element3edges(elem2coarse(:,2),[1,3]), element3edges(elem2coarse(:,3),[1,3])];
markedEdges = zeros(size(edge2nodes,1),1);
markedEdges(irregular2edges(:,1)) = 1;
elem2coarse(find(sum(reshape(markedEdges(adm2edges),[],6),2)),:)=[]; 
%*** mark nodes to be eliminated
markedNodes = zeros(nC,1);
markedNodes(elements(elem2coarse(:,1:3),2)) = 1;
%*** create new elements
idx = setdiff(1:nE,unique(elem2coarse));
newelements = [elements(idx,:); reshape(elements(elem2coarse(:,1:3),1),[],3)];
%*** update constrains, allow redundancies
irregular = [irregular; ...
    [elements(elem2coarse(:,1),1),elements(elem2coarse(:,2),1),elements(elem2coarse(:,1),2)];...
    [elements(elem2coarse(:,2),1),elements(elem2coarse(:,3),1),elements(elem2coarse(:,2),2)];...
    [elements(elem2coarse(:,3),1),elements(elem2coarse(:,1),1),elements(elem2coarse(:,3),2)]];
irregular = sortrows(irregular,3);
idx = find(irregular(1:end-1,3) == irregular(2:end,3));
irregular([idx(:);idx(:)+1],:) = [];
boundinfo = edge2nodes(accumarray([element3edges(:);irregular2edges(:)],1,[length(edge2nodes) 1])==1,:);
node2bound = zeros(nC,1);
node2bound(unique(boundinfo(:)),:) = 1;
node2constr = zeros(nC,1);
node2constr(irregular(:,3),1) = 1:size(irregular,1);
idx = (markedNodes & node2bound);
irregular(node2constr(idx),:) = [];
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
        varargout{j} = boundary(find(boundary(:,2)),:);
    else
        varargout{j} = [];
    end
end
%*** update coordinates and element numbering
i2fi = zeros(nC,1);
i2fi(newelements) = 1;
coordinates = coordinates(find(i2fi),:);
i2fi = cumsum(i2fi);
elements = reshape(i2fi(newelements),[],3);
irregular = reshape(i2fi(irregular),[],3);
for j = 1:nargout-3
    varargout{j} = i2fi(varargout{j});
end
end

