function [coordinates,elements,varargout] = ...
    TcoarsenRG(N0,coordinates,elements,varargin)
%TcoarsenRG: local coarsening of triangular mesh emerged from  
%          TrefineRG where marked elements are coarsened by inverting the
%          red-green refinement
%
%Usage:
%
% [coordinates,elements3,dirichlet,neumann] ...
%    = TcoarsenRG(N0,coordinates,elements3,dirichlet,neumann,marked)
% or
%
% [coordinates,elements3] ...
%    = TcoarsenRG(N0,coordinates,elements3,marked)
%
%Comments:
%
%    TcoarsenRG expects as input a mesh described by the 
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

%*** remove green elements
[coordinates,elements,irregular,marked] = coarse_greenelements(coordinates,elements,varargin{end});
%*** coarse irregular mesh
[coordinates,elements,irregular,varargout{1:nargin-4}] = ...
    TcoarsenR_2neighbor(N0,coordinates,elements,irregular,varargin{1:end-1},marked);
%*** closure --> remove irregular edges with green patterns
[coordinates, elements] = regularize_Tgreen(coordinates, elements, irregular);
end