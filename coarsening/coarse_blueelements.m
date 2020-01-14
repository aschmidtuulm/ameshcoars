function [coordinates,elements,irregular,varargout] = coarse_blueelements(coordinates,elements,varargin)
global nB
marked = varargin{end};
nE = size(elements,1);
varargout{nargout-3} = marked(marked<=nE-nB); % delete marking of mesh closure
bdx = (nE-nB+1):3:nE;
%*** Obtain geometric information on edges
[edge2nodes,~,element2edges,boundary2edges{1:nargin-4}] ...
    = provideGeometricData(zeros(0,3),elements,varargin{1:end-1});
boundinfo = edge2nodes(accumarray(element2edges(:),1,[length(edge2nodes) 1])==1,:);
irregular = [elements(bdx+1,1), elements(bdx+2,[2,1]);...
    elements(bdx+2,2), elements(bdx,[1,4])];
%*** delete duplicates
irregular = sortrows(irregular,3);
idx = find(irregular(1:end-1,3)== irregular(2:end,3));
irregular([idx,idx+1],:) = [];
%*** delte boundary edges & correct boundaries
[~,idx] = setdiff(irregular(:,3), boundinfo(:,2));
%*** correct boundaries
jdx = setdiff(1:size(irregular,1),idx);
for j = 1:nargout-4
    boundary = varargin{j};
    if ~isempty(boundary)
        for i = 1:length(jdx) 
            [kdx,~] = find(boundary==irregular(jdx(i),3));
            boundary = [boundary; [boundary(kdx(2),1), boundary(kdx(1),2)]];
            boundary(kdx,:)=[];
        end
    end
    varargout{j} = boundary;
end
irregular = irregular(idx,:);
%*** define new elements
elements = [elements; elements(bdx,[1,2]), elements(bdx+1,1), elements(bdx+2,2)];
elements([bdx,bdx+1,bdx+2],:) =[];
nB = 0;
%*** sort
[~,jdx] = min(elements,[],2);
idx = jdx == 3;
elements(idx,:) = elements(idx,[3,4,1,2]);
end