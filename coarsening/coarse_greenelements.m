function [coordinates,elements,irregular,varargout] = coarse_greenelements(coordinates,elements,varargin)
global nG
marked = varargin{end};
nE = size(elements,1);
varargout{nargout-3} = marked(marked<=nE-nG); % delete marking of mesh closure
%*** determine green elements
gdx = (nE-nG+1):2:nE;
%*** define irregular edges
irregular = [elements(gdx,2), elements(gdx+1,2), elements(gdx+1,1)];
%*** define new elements
elements = [elements; elements(gdx,2), elements(gdx+1,2), elements(gdx+1,3)];
elements([gdx,gdx+1],:) =[];
nG = 0;
end

