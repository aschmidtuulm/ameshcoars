function [coordinates,elements,varargout] = TcoarsenNVB(N0,coordinates,elements,varargin) 

%coarsenNVB: coarsening of finite element mesh by reverse newest vertex 
%            bisection
%
%Usage:
%
%[coordinates,elements,dirichlet,neumann] ...
%    = TcoarsenNVB(n0,coordinates,elements,dirichlet,neumann,marked)
%
%Comments:
%
%    TcoarsenNVB expects as input a finite element mesh described by the fields 
%    coordinates, elements, dirichlet and neumann. Vertices with index smaller
%    than N0 belong to the initial mesh and may not be removed. The vector
%    MARKED contains the indices of elements which shall be coarsened (if 
%    possible).
%
%
%Remark:
%
%    This program is a supplement to the paper 
%    >> Efficient Implementation of Adaptive P1-FEM in Matlab <<
%    by S. Funken, D. Praetorius, and P. Wissgott. The reader should 
%    consult that paper for more information.   
%
%Authors:
% 
%    S. Funken, D. Praetorius, P. Wissgott  10-07-08

nC = size(coordinates,1); 
nE = size(elements,1);
%*** Obtain geometric information on neighbouring elements
I = elements(:);
J = reshape(elements(:,[2,3,1]),3*nE,1);
nodes2edge = sparse(I,J,1:3*nE);
mask = nodes2edge>0;
[foo{1:2},idxIJ] = find( nodes2edge );
[foo{1:2},neighbourIJ] = find( mask + mask.*sparse(J,I,[1:nE,1:nE,1:nE]') );
element2neighbours(idxIJ) = neighbourIJ - 1;
element2neighbours = reshape(element2neighbours,nE,3);
%*** Determine which nodes (created by refineNVB) are deleted by coarsening
marked = zeros(nE,1);
marked(varargin{end}) = 1;
newestNode = unique(elements((marked & elements(:,3)>N0),3));
valence = accumarray(elements(:),1,[nC 1]);
markedNodes = zeros(nC,1); 
markedNodes(newestNode((valence(newestNode) == 2 | valence(newestNode) == 4))) = 1;
%*** Collect pairs of brother elements that will be united
idx = find(markedNodes(elements(:,3)) & (element2neighbours(:,3) > (1:nE)'))';
markedElements = zeros(nE,1); 
markedElements(idx) = 1;
for element = idx
    if markedElements(element)
        markedElements(element2neighbours(element,3)) = 0;
    end
end
idx = find(markedElements);
%*** Coarsen two brother elements
brother = element2neighbours(idx,3);
elements(idx,[1 3 2]) = [elements(idx,[2 1]) elements(brother,1)];
%*** Delete redundant nodes
activeNodes = find(~markedNodes);
coordinates = coordinates(activeNodes,:);
%*** Provide permutation of nodes to correct further data
coordinates2newCoordinates = zeros(1,nC);
coordinates2newCoordinates(activeNodes) = 1:length(activeNodes);
%*** Delete redundant elements + correct elements
elements(brother,:) = [];                 
elements = coordinates2newCoordinates(elements);
%*** Delete redundant boundaries + correct boundaries
for j = 1:nargout-2;
    boundary = varargin{j};
    if ~isempty(boundary)
      node2boundary = zeros(nC,2);
      node2boundary(boundary(:,1),1) = 1:size(boundary,1);
      node2boundary(boundary(:,2),2) = 1:size(boundary,1);
      idx = ( markedNodes & node2boundary(:,2) );
      boundary(node2boundary(idx,2),2) = boundary(node2boundary(idx,1),2);
      boundary(node2boundary(idx,1),2) = 0;
      varargout{j} = coordinates2newCoordinates(boundary(find(boundary(:,2)),:));
    else
      varargout{j} = [];
    end
end