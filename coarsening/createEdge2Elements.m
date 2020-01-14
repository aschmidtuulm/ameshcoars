function edge2elements = createEdge2Elements(element2edges)

nT = size(element2edges,1);  % number of elements (triangles)
nE = max(element2edges(:));  % number of edges
edge2elements = zeros(nE,2);
ind = sortrows(reshape([element2edges(:), ...
               reshape( (1:nT)'*ones(1,3),[],1)],[],2));
flag = [false;ind(1:end-1,1) == ind(2:end,1)];
edge2elements(ind(~flag,1),1) = ind(~flag,2);
edge2elements(ind(flag,1),2) = ind(flag,2);
