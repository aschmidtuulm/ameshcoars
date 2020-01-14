function [coordinates,elements,varargout] = ...
    QcoarsenRB(N0,coordinates,elements,varargin)
%QcoarsenRB: local coarsening of quadrilateral mesh emerged from  
%          QrefineRB where marked elements are coarsened by inverting the
%          red-blue refinement
%
%Usage:
%
% [coordinates,elements4,dirichlet,neumann] ...
%    = QcoarsenRB(N0,coordinates,elements4,dirichlet,neumann,marked)
% or
%
% [coordinates,elements4] ...
%    = QcoarsenRB(N0,coordinates,elements4,marked)
%
%Comments:
%
%    QcoarsenRB expects as input a mesh described by the 
%    fields coordinates, elements4, irregular, dirichlet (optional) and neumann 
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

%*** remove green elements
[coordinates,elements,irregular,boundary{1:nargin-4},marked] = coarse_blueelements(coordinates,elements,varargin{1:end});
%*** coarse irregular mesh
[coordinates,elements,irregular,boundary{1:nargin-4}] = ...
    QcoarsenR(N0,coordinates,elements,irregular,boundary{1:end},marked);
%*** closure --> remove irregular edges with green patterns
[coordinates, elements,varargout{1:nargin-4}] = regularize_Qblue(coordinates, elements, irregular,boundary{1:end});
end