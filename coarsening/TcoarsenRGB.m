function [coordinates,elements,varargout] = TcoarsenRGB(N0,coordinates,elements,varargin)
%TcoarsenRGB: local coarsening of triangular mesh emerged from  
%          TrefineRGB where marked elements are coarsened by inverting the
%          red-green-blue refinement
%
%Usage:
%
% [coordinates,elements3,dirichlet,neumann] ...
%    = TcoarsenRGB(N0,coordinates,elements3,dirichlet,neumann,marked)
% or
%
% [coordinates,elements3] ...
%    = TcoarsenRGB(N0,coordinates,elements3,marked)
%
%Comments:
%
%    TcoarsenRGB expects as input a mesh described by the 
%    fields coordinates, elements3, dirichlet (optional) and neumann 
%    (optional). N0 specifies the number of coordinates in the initial mesh
%    as the routine does not eliminate nodes from the initial mesh.
%    The vector marked contains the indices of elements which shall be
%    coarsened. A family of elements is coarsened to their common father
%    element when at least one of the elements is marked for coarsening.
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

nC = size(coordinates,1); nE = size(elements,1);
marked = varargin{end};
%*** obtain geometric information
[~,element2edges] = provideGeometricData(elements,zeros(0,4));
edge2elements = createEdge2Elements_adv(element2edges);
%*** determine admissible nodes for coarsening
midelement = find(sum(abs(element2edges(4:end,:)-[element2edges(3:end-1,1),...
   element2edges(1:end-3,2),element2edges(2:end-2,3)]),2)==0)+3;
newest_nodes = zeros(nC,1); marked_nodes = zeros(nC,1);
newest_nodes(elements(:,3)) = 1; newest_nodes(1:N0)=0;
marked_nodes(elements(marked,:)) = 1;
%*** define color of nodes
mid_nodes = elements(midelement,:);
red_node = accumarray(reshape(mid_nodes,[],1),1);
valence = accumarray(elements(:),ones(3*nE,1),[nC,1])-[red_node;zeros(nC-size(red_node,1),1)];
newest_marked = elements(marked,3);
green_node = setdiff(newest_marked(newest_marked>N0),find(red_node));
green_val = valence(green_node);
%*** determine arising hanging nodes
blocked_nodes = ~(((valence == 2) | (valence == 4)) & newest_nodes & marked_nodes);
old = -1;
while ~isempty(find(old ~=blocked_nodes))
    old = blocked_nodes;
    [row,~] = find(reshape(blocked_nodes(mid_nodes),[],3));
    blocked_nodes(elements(midelement(row),3)) = 1; % mark reference edge
end
hanging_nodes = blocked_nodes & newest_nodes;
%*** determine coarsening pattern regarding hanging nodes
red_elements = reshape(hanging_nodes(mid_nodes),[],3);
val = sum(red_elements .* [1,2,4],2);
none = find(val==0);
gdx = find(val == 4);
b_rdx = find(val==5);
b_ldx = find(val==6);
%*** define new elements of former red pattern regarding hanging nodes
mark_red = reshape(mid_nodes,[],1);
mark_mid = edge2elements(element2edges(midelement,:));
elements(midelement(none)-1,:)=[elements(mark_mid(none,2),1),elements(mark_mid(none,3),2),elements(mark_mid(none,1),3)];
elements([midelement(gdx)-3,midelement(gdx)-2],:) = [elements(mark_mid(gdx,1),3),elements(mark_mid(gdx,2),[1,2]);...
    elements(mark_mid(gdx,3),2),elements(mark_mid(gdx,1),3),elements(mark_mid(gdx,3),1)];
elements([midelement(b_ldx)-3,midelement(b_ldx)-2,midelement(b_ldx)-1],:) = [elements(mark_mid(b_ldx,2),2),elements(mark_mid(b_ldx,1),[3,1]);...
    elements(mark_mid(b_ldx,2),:);elements(mark_mid(b_ldx,3),2),elements(mark_mid(b_ldx,1),3),elements(mark_mid(b_ldx,3),1)];
elements([midelement(b_rdx)-3,midelement(b_rdx)-2,midelement(b_rdx)-1],:) = [elements(mark_mid(b_rdx,1),3),elements(mark_mid(b_rdx,2),[1,2]);...
    elements(mark_mid(b_rdx,3),:);elements(mark_mid(b_rdx,1),3),elements(mark_mid(b_rdx,3),[1,3])];
%*** define new elements of former green pattern neighbouring a red pattern
green_remove = mark_red(red_node(mark_red)==1 & valence(mark_red)==4 & ~hanging_nodes(mark_red)); 
non_redelement = setdiff(1:size(elements,1),[midelement-3;midelement-2;midelement-1;midelement]);
dummy = zeros(nC,1); dummy(green_remove) = 1;
green_elem = find(dummy(elements(non_redelement,3)));
elements(non_redelement(green_elem(1:2:end-1)),:) = [elements(non_redelement(green_elem(1:2:end-1)),2),elements(non_redelement(green_elem(2:2:end)),1:2)];
%*** define new elements of former green pattern 
mark_green = green_node(green_val == 4 | green_val == 2);
mark_green_bound = green_node(green_val == 2);
if ~isempty(mark_green)
    dummy = zeros(nC,1); dummy(mark_green) = 1;
    green_bound = zeros(nC,1); green_bound(mark_green_bound) = 1;
    greengreen = find(dummy(elements(:,3)));
    check_neighbour = greengreen(1:end-1) == greengreen(2:end)-1;
    for i = 1:size(check_neighbour,1)-1
        if check_neighbour(i) == 1
            check_neighbour(i+1) = 0;
        end
    end
    greengreen(logical([~check_neighbour;1])) = [];
else
    greengreen = []; green_bound = [];
end
elements(greengreen,:) = [elements(greengreen,2),elements(greengreen+1,1:2)];
%*** delete old triangulation
elements([midelement(none)-3;midelement(none)-2;midelement(none);midelement(gdx)-1;midelement(gdx);midelement(b_ldx);midelement(b_rdx);...
    non_redelement(green_elem(2:2:end))';greengreen+1],:) = [];
%*** correct boundaries
bound_nodes = intersect(find(green_bound),elements(:));
not_to_eliminate = zeros(nC,1); not_to_eliminate(bound_nodes) = 1;
    for j = 1:nargout-2
        boundary = varargin{j};
        if ~isempty(boundary)
            node2boundary = zeros(nC,2);
            node2boundary(boundary(:,1),1) = 1:size(boundary,1);
            node2boundary(boundary(:,2),2) = 1:size(boundary,1);
            idx = (~blocked_nodes) & logical(node2boundary(:,2))& logical(~not_to_eliminate);
            boundary(node2boundary(idx,2),2) = boundary(node2boundary(idx,1),2);
            boundary(node2boundary(idx,1),2) = 0;
            varargout{j} = boundary(find(boundary(:,2)),:);
        else
            varargout{j} = [];
        end
    end
%*** update coordinates and elements
i2fi = zeros(nC,1);
i2fi(elements) = 1;
coordinates = coordinates(find(i2fi),:);
i2fi = cumsum(i2fi);
elements = reshape(i2fi(elements),[],3);
for i = 1:nargout-2
    varargout{i} = i2fi(varargout{i});
end
end