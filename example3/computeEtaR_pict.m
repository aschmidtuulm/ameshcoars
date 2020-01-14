function etaR = computeEtaR_pict(pict,elements,coordinates)
    etaR= zeros(1,size(elements,1));
   
    
    max_x = max(coordinates(:,1));
    min_x = min(coordinates(:,1));
    max_y = max(coordinates(:,2));
    min_y = min(coordinates(:,2));
    
    [m,n] = size(pict);
    
    for i=1:size(elements,1)
        min_tri_x = min(coordinates(elements(i,:),1));
        max_tri_x = max(coordinates(elements(i,:),1));
        min_tri_y = min(coordinates(elements(i,:),2));
        max_tri_y = max(coordinates(elements(i,:),2));
        min_tri_x_idx = floor((min_tri_x - min_x)/((max_x-min_x)/n))+1;
        max_tri_x_idx = ceil((max_tri_x - min_x)/((max_x-min_x)/n));
        min_tri_y_idx = floor((min_tri_y - min_y)/((max_y-min_y)/m))+1;
        max_tri_y_idx = ceil((max_tri_y - min_y)/((max_y-min_y)/m));
        
        xpts = linspace(min_x + min_tri_x_idx*((max_x-min_x)/n),min_x + max_tri_x_idx*((max_x-min_x)/n),max_tri_x_idx - min_tri_x_idx + 1);
        ypts = linspace(min_y + min_tri_y_idx*((max_y-min_y)/m),min_y + max_tri_y_idx*((max_y-min_y)/m),max_tri_y_idx - min_tri_y_idx + 1);
        [ypts,xpts] = meshgrid(ypts,xpts);
        [xn,xm] = size(xpts);
        xpts = reshape(xpts,1,xn*xm);
        ypts = reshape(ypts,1,xn*xm);
        %p = [xpts;ypts];
        %inside = PointInTriangle2(p,elements(i,:),coordinates);
        inside = inpolygon(xpts,ypts,coordinates(elements(i,:),1),coordinates(elements(i,:),2));
        if sum(inside) > 0
            %min_tri_y_idx:max_tri_y_idx
            %min_tri_x_idx:max_tri_x_idx
            submat = pict(min_tri_y_idx:max_tri_y_idx,min_tri_x_idx:max_tri_x_idx);
            [subm,subn] = size(submat);
            submat = reshape(submat',subm*subn,1)';
            etaR(i) = max(submat(inside>0))-min(submat(inside>0));
        else
            etaR(i) = 0;
        end
    etaR = etaR';
end