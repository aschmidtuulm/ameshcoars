%test_moving_circle: refine and coarsen given mesh along a moving circle
%                   with different mesh refinement strategies
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
%    S. Funken, A. Schmidt  14-01-20

clear all, close all, format compact
figure(1), clf

addpath('../refinement')
addpath('../coarsening')

%% Define initial mesh and circle
c = [0,0;1,0;1,1;0,1;2,0;2,1];
e3 = [3,1,2;1,3,4;2,6,3;6,2,5];
e4 = [1,2,3,4;2,5,6,3];
b = [1,2;2,5;5,6;6,3;3,4;4,1];
N0 = size(c,1);
N_min = 1e2;
N_max = 1e4;

C = [0.5,0.7];%[0.8,0.7];%[0.5,0.7];
R = 0.4;%0.5;%0.15;
h = 2.5e-3;

s = linspace(0,1.5,30);

%% *** Do TrefineNVB
coordinates = c; elements = e3; boundary = b;
if ~exist('TNVB', 'dir')
    mkdir('TNVB')
end

for k=1:length(s)
    while 1
        [marked,~] = markCircle(coordinates,elements,~[],[C(1)+s(k),C(2)-0.6*s(k)],R,h);
        marked = find(marked);
        [coordinates,elements,boundary] ...
            = TrefineNVB(coordinates,elements,boundary,marked);
        if isempty(marked)|| (size(coordinates,1)>N_max)
            break
        end
    end
    
    
    clf
    patch('Faces',elements,'Vertices',coordinates,'Facecolor','none')
    axis equal
    axis off
    print(figure(1),['TNVB\MeshTrefineNVB' num2str(k)],'-depsc2')
    
    %% *** Do TcoarsenNVB
    while 1
        marked = 1:size(elements,1);
        [coordinates,elements,boundary] = TcoarsenNVB(N0,coordinates,elements,boundary,marked);
        if size(coordinates,1)<N_min
            break
        end
    end
    
end

% *** Do TrefineRGB
coordinates = c; elements = e3; boundary = b;
if ~exist('TRGB', 'dir')
    mkdir('TRGB')
end

for k=1:length(s)
    while 1
        [marked,~] = markCircle(coordinates,elements,~[],[C(1)+s(k),C(2)-0.6*s(k)],R,h);
        marked = find(marked);
        [coordinates,elements,boundary] ...
            = TrefineRGB(coordinates,elements,boundary,marked);
        if isempty(marked)|| (size(coordinates,1)>N_max)
            break
        end
    end
    
    
    clf
    patch('Faces',elements,'Vertices',coordinates,'Facecolor','none')
    axis equal
    axis off
    %print(figure(1),['TRGB\MeshTrefineRGB' num2str(k)],'-depsc2')
    
    % *** Do TcoarsenRGB
    while 1
        marked = 1:size(elements,1);
        [coordinates,elements,boundary] = TcoarsenRGB(N0,coordinates,elements,boundary,marked);
        if size(coordinates,1)<N_min
            break
        end
    end
    
end


%% *** Do TrefineRG
coordinates = c; elements = e3; boundary = b;
if ~exist('TRG', 'dir')
    mkdir('TRG')
end

for k=1:length(s)
    while 1
        [marked,~] = markCircle(coordinates,elements,~[],[C(1)+s(k),C(2)-0.6*s(k)],R,h);
        marked = find(marked);
        [coordinates,elements,boundary] ...
            = TrefineRG(coordinates,elements,boundary,marked);
        if isempty(marked)|| (size(coordinates,1)>N_max)
            break
        end
    end
    
    
    clf
    patch('Faces',elements,'Vertices',coordinates,'Facecolor','none')
    axis equal
    axis off
    print(figure(1),['TRG\MeshTrefineRG' num2str(k)],'-depsc2')
    
    %% *** Do TcoarsenRG
    while 1
        marked = 1:size(elements,1);
        [coordinates,elements,boundary] = ...
            TcoarsenRG(N0,coordinates,elements,boundary,marked);
        if size(coordinates,1)<N_min
            break
        end
    end
    
end

%% *** Do TrefineR
coordinates = c; elements = e3; irregular = zeros(0,3); boundary = b;
if ~exist('TR', 'dir')
    mkdir('TR')
end

for k=1:length(s)
    while 1
        [marked,~] = markCircle(coordinates,elements,[],[C(1)+s(k),C(2)-0.6*s(k)],R,h);
        marked = find(marked);
        [coordinates,elements,irregular,boundary] ...
            = TrefineR(coordinates,elements,irregular,boundary,marked);
        if isempty(marked)|| (size(coordinates,1)>N_max)
            break
        end
    end
    
    
    clf
    patch('Faces',elements,'Vertices',coordinates,'Facecolor','none')
    axis equal
    axis off
    print(figure(1),['TR\MeshTrefineR' num2str(k)],'-depsc2')
    
    %% *** Do TcoarsenR
    while 1
        marked = 1:size(elements,1);
        [coordinates,elements,irregular,boundary] =...
            TcoarsenR(N0,coordinates,elements,irregular,boundary,marked);
        if size(coordinates,1)<N_min
            break
        end
    end
    
end

%% *** Do QrefineR
coordinates = c; elements = e4; irregular = zeros(0,3); boundary = b;
if ~exist('QR', 'dir')
    mkdir('QR')
end

for k=1:length(s)
    while 1
        [~,marked] = markCircle(coordinates,[],elements,[C(1)+s(k),C(2)-0.6*s(k)],R,h);
        marked = find(marked);
        [coordinates,elements,irregular,boundary] ...
            = QrefineR(coordinates,elements,irregular,boundary,marked);
        if isempty(marked)|| (size(coordinates,1)>N_max)
            break
        end
    end
    
    
    clf
    patch('Faces',elements,'Vertices',coordinates,'Facecolor','none')
    axis equal
    axis off
    print(figure(1),['QR\MeshQrefineR' num2str(k)],'-depsc2')
    
    %% *** Do QcoarsenR
    while 1
        marked = 1:size(elements,1);
        [coordinates,elements,irregular,boundary] = ...
            QcoarsenR(N0,coordinates,elements,irregular,boundary,marked);
        if size(coordinates,1)<N_min
            break
        end
    end
    
end

%% *** Do QrefineRB
coordinates = c; elements = e4; boundary = b;
if ~exist('QRB', 'dir')
    mkdir('QRB')
end

for k=1:length(s)
    while 1
        [~,marked] = markCircle(coordinates,[],elements,[C(1)+s(k),C(2)-0.6*s(k)],R,h);
        marked = find(marked);
        [coordinates,elements,boundary] ...
            = QrefineRB(coordinates,elements,boundary,marked);
        if isempty(marked)|| (size(coordinates,1)>N_max)
            break
        end
    end
    
    
    clf
    patch('Faces',elements,'Vertices',coordinates,'Facecolor','none')
    axis equal
    axis off
    print(figure(1),['QRB\MeshQrefineRB' num2str(k)],'-depsc2')
    
    %% *** Do QcoarsenRB
    while 1
        marked = 1:size(elements,1);
        [coordinates,elements,boundary] =...
            QcoarsenRB(N0,coordinates,elements,boundary,marked);
        if size(coordinates,1)<N_min
            break
        end
    end
    
end

%% *** Do QrefineRG
coordinates = c; elements3 = zeros(0,3); elements4 = e4; boundary = b;
if ~exist('QRG', 'dir')
    mkdir('QRG')
end

for k=1:length(s)
    while 1
        [marked3,marked4] = markCircle(coordinates,elements3,elements4,[C(1)+s(k),C(2)-0.6*s(k)],R,h);
        marked3 = find(marked3);
        marked4 = find(marked4);
        %*** Refine mesh
        [coordinates,elements4,marked,irregular]...
            =recoarseedges_tri(coordinates,elements3,elements4,marked3,marked4);
        
        [coordinates,elements4,irregular] ...
            = QrefineR(coordinates,elements4,irregular,marked);
        
        [coordinates,elements4,elements3]...
            =regularizeedges_tri(coordinates,elements4,irregular);
        if (isempty(marked4) && isempty(marked3)) || size(coordinates,1)>N_max
            break
        end
    end
    
    
    clf
    patch('Faces',elements4,'Vertices',coordinates,'Facecolor','none')
    hold on
    patch('Faces',elements3,'Vertices',coordinates,'Facecolor','none')
    axis equal
    axis off
    print(figure(1),['QRG\MeshQrefineRG' num2str(k)],'-depsc2')
    
    %% *** Do QcoarsenRG
    while 1
        marked = 1:size(elements4,1);
        [coordinates,elements3,elements4] = QcoarsenRG(N0,coordinates,elements3,elements4,marked);
        if size(coordinates,1)<N_min
            break
        end
    end
    
end