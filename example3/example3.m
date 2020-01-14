%example3:  triangulate a GIF using a gradient based indicator via
%           refinement and coarsening
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

clear all, clc,close all


addpath ../refinement
addpath ../coarsening

v=VideoReader('giphy.mp4');

numberofframes = v.Duration*v.FrameRate;
video = cell(int8(numberofframes));
m = v.Height;
n = v.Width;
i =1;
while hasFrame(v)
    video{i} = rgb2gray(readFrame(v));
    %imshow(video{i});
    i= i+1;
end

pic = video{1};
pic = pic(m:-1:1,:);

%*** initial coordinates, elements, boundary
c = [1,1;n,1;n,m;1,m];
e3 = [2,4,1;4,2,3];
e4 = [1,2,3,4];
b = [1,2;2,3;3,4;4,1];
t = 10;
tol = 0.03;
sigma = 0.25;
N0 = size(c,1);
N_min = 1e3;
N_max = 1e4;

f = figure;
set(f, 'Units', 'normalized', 'Position', [0.2, 0.1, 0.8, 0.75]);

%% *** Do TrefineRGB
coordinates = c; elements = e3; boundary = b;
if ~exist('TRGB', 'dir')
    mkdir('TRGB')
end
it = 0;
while 1
    it = it+1;
    fprintf('******** TrefineRGB ******** STEP  %d **********\n',it);
    %*** compute refinement indicator
    etaR = computeEtaR_pict(pic,elements,coordinates);
    if etaR <= tol^2
        break
    end
    %*** mark elements
    marked = markElementsDoerfler_pict(etaR,0.5);
    %*** refine mesh
    [coordinates,elements,boundary]= TrefineRGB(coordinates,elements,boundary,marked);
    % Stopping criterion
    if isempty(marked)|| (size(coordinates,1)>N_max)
        break
    end
end
clf
patch('Faces',elements,'Vertices',coordinates,'FaceColor','none')
axis off
print(figure(1),['TRGB\MeshTrefineRGB' num2str(1)],'-depsc2')

for i=2:int8(numberofframes)
    pic = video{i};
    pic = pic(m:-1:1,:);
    it = 2*it;
    % coarsen
    while it
        it = it - 1;
        fprintf('******** TcoarsenRGB ******** STEP  %d **********\n',it);
        %*** Mark elements for coarsening
        marked = 1:size(elements,1);
        %*** Try to coarse mesh
        [coordinates,elements,boundary] = ...
            TcoarsenRGB(N0,coordinates,elements,boundary,marked);
        if size(coordinates,1)<N_min || it == 0
            break
        end
    end
    % refine
    it = 0;
    while 1
        it = it+1;
        fprintf('******** TrefineRGB ******** STEP  %d **********\n',it);
        %*** compute refinement indicator
        etaR = computeEtaR_pict(pic,elements,coordinates);
        if etaR <= tol^2
            break
        end
        %*** mark elements
        marked = markElementsDoerfler_pict(etaR,0.5);
        %*** refine mesh
        [coordinates,elements,boundary]= TrefineRGB(coordinates,elements,boundary,marked);
        if isempty(marked)|| (size(coordinates,1)>N_max)
            break
        end
    end
    clf
    patch('Faces',elements,'Vertices',coordinates,'FaceColor','none')
    axis off
    print(figure(1),['TRGB\MeshTrefineRGB' num2str(i)],'-depsc2')
end

%% *** Do TrefineNVB
coordinates = c; elements = e3; boundary = b;
if ~exist('TNVB', 'dir')
    mkdir('TNVB')
end
it = 0;
while 1
    it = it + 1;
    fprintf('******** TrefineNVB ******** STEP  %d **********\n',it);
    %*** compute refinement indicator
    etaR = computeEtaR_pict(pic,elements,coordinates);
    if etaR <= tol^2
        break
    end
    %*** mark elements
    marked = markElementsDoerfler_pict(etaR,0.5);
    %*** refine mesh
    [coordinates,elements,boundary]= TrefineNVB(coordinates,elements,boundary,marked);
    % Stopping criterion
    if isempty(marked)|| (size(coordinates,1)>N_max)
        break
    end
end
clf
patch('Faces',elements,'Vertices',coordinates,'FaceColor','none')
axis off
print(figure(1),['TNVB\MeshTrefineNVB' num2str(1)],'-depsc2')

for i=2:int8(numberofframes)
    pic = video{i};
    pic = pic(m:-1:1,:);
    it = 2*it;
    % coarsen
    while it
        it = it -1;
        fprintf('******** TcoarsenNVB ******** STEP  %d **********\n',it);
        %*** Mark elements for coarsening
        marked = 1:size(elements,1);
        %*** Try to coarse mesh
        [coordinates,elements,boundary] = ...
            TcoarsenNVB(N0,coordinates,elements,boundary,marked);
        if size(coordinates,1)<N_min || it == 0
            break
        end

    end
    % refine
    it = 0;
    while 1
        it = it + 1;
        fprintf('******** TrefineNVB ******** STEP  %d **********\n',it);
        %*** compute refinement indicator
        etaR = computeEtaR_pict(pic,elements,coordinates);
        if etaR <= tol^2
            break
        end
        %*** mark elements
        marked = markElementsDoerfler_pict(etaR,0.5);
        %*** refine mesh
        [coordinates,elements,boundary]= TrefineNVB(coordinates,elements,boundary,marked);
        if isempty(marked)|| (size(coordinates,1)>N_max)
            break
        end
    end
    clf
    patch('Faces',elements,'Vertices',coordinates,'FaceColor','none')
    axis off
    print(figure(1),['TNVB\MeshTrefineNVB' num2str(i)],'-depsc2')
end


%% *** Do TrefineRG
coordinates = c; elements = e3; boundary = b; pic = video{1}; pic = pic(m:-1:1,:);
if ~exist('TRG', 'dir')
    mkdir('TRG')
end
it = 0;
while 1
    it = it +1;
    fprintf('******** TrefineRG ******** STEP  %d **********\n',it);
    %*** compute refinement indicator
    etaR = computeEtaR_pict(pic,elements,coordinates);
    if etaR <= tol^2
        break
    end
    %*** mark elements
    marked = markElementsDoerfler_pict(etaR,0.5);
    %*** refine mesh
    [coordinates,elements,boundary]= TrefineRG(coordinates,elements,boundary,marked);
    % Stopping criterion
    if isempty(marked)|| (size(coordinates,1)>N_max)
        break
    end
end
clf
patch('Faces',elements,'Vertices',coordinates,'FaceColor','none')
axis off
print(figure(1),['TRG\MeshTrefineRG' num2str(1)],'-depsc2')

for i=2:int8(numberofframes)
    pic = video{i};
    pic = pic(m:-1:1,:);
    it = 2 * it;
    while it
        it = it - 1;
        fprintf('******** TcoarsenRG ******** STEP  %d **********\n',it);
        %*** Mark elements for coarsening
        marked = 1:size(elements,1);
        %*** Try to coarse mesh
        [coordinates,elements,boundary] = ...
            TcoarsenRG(N0,coordinates,elements,boundary,marked);
        if size(coordinates,1)<N_min || it == 0
            break
        end

    end
    it = 0;
    while 1
        it = it + 1;
        fprintf('******** TrefineRG ******** STEP  %d **********\n',it);
        %*** compute refinement indicator
        etaR = computeEtaR_pict(pic,elements,coordinates);
        if etaR <= tol^2
            break
        end
        %*** mark elements
        marked = markElementsDoerfler_pict(etaR,0.5);
        %*** refine mesh
        [coordinates,elements,boundary]= TrefineRG(coordinates,elements,boundary,marked);
        if isempty(marked)|| (size(coordinates,1)>N_max)
            break
        end
    end
    clf
    patch('Faces',elements,'Vertices',coordinates,'FaceColor','none')
    axis off
    print(figure(1),['TRG\MeshTrefineRG' num2str(i)],'-depsc2')
end

%% *** Do QrefineRG
coordinates = c; elements3 = zeros(0,3); elements4 = e4; boundary = b;  pic = video{1}; pic = pic(m:-1:1,:);
if ~exist('QRG', 'dir')
    mkdir('QRG')
end
it = 0;
while 1
    it = it + 1;
    fprintf('******** QrefineRG ******** STEP  %d **********\n',it);
    %*** compute refinement indicator
    etaR3 = computeEtaR_pict(pic,elements3,coordinates);
    etaR4 = computeEtaR_pict(pic,elements4,coordinates);
    %*** mark elements
    mark3 = markElementsDoerfler_pict(etaR3,0.5);
    mark4 = markElementsDoerfler_pict(etaR4,0.5);
    %*** refine mesh
    [coordinates,elements4,marked,irregular] = ...
        recoarseedges_tri(coordinates,elements3,elements4,mark3,mark4);
    [coordinates,elements4,irregular,boundary] = ...
        QrefineR(coordinates,elements4,irregular,boundary,marked);
    [coordinates,elements4,elements3] = ...
        regularizeedges_tri(coordinates,elements4,irregular);
    % Stopping criterion
    if isempty(marked)|| (size(coordinates,1)>N_max)
        break
    end
end
clf
patch('Faces',elements4,'Vertices',coordinates,'FaceColor','none')
hold on
patch('Faces',elements3,'Vertices',coordinates,'FaceColor','none')
axis off
print(figure(1),['QRG\MeshQrefineRG' num2str(1)],'-depsc2')

for i=2:int8(numberofframes)
    pic = video{i};
    pic = pic(m:-1:1,:);
    it = 2*it;
    % coarsen
    while 1
        it = it - 1;
        fprintf('******** QcoarsenRG ******** STEP  %d **********\n',it);
        %*** Mark elements for coarsening
        marked = 1:size(elements,1);
        %*** Try to coarsen mesh
        [coordinates,elements3,elements4,boundary] = ...
            QcoarsenRG(N0,coordinates,elements3,elements4,boundary,marked);
        if size(coordinates,1)<N_min || it == 0
            break
        end

    end
    it = 0;
    % refine
    while 1
        it = it + 1;
        fprintf('******** QrefineRG ******** STEP  %d **********\n',it);
        %*** compute refinement indicator
        etaR3 = computeEtaR_pict(pic,elements3,coordinates);
        etaR4 = computeEtaR_pict(pic,elements4,coordinates);
        %*** mark elements
        mark3 = markElementsDoerfler_pict(etaR3,0.5);
        mark4 = markElementsDoerfler_pict(etaR4,0.5);
        %*** refine mesh
        [coordinates,elements4,marked,irregular] = ...
            recoarseedges_tri(coordinates,elements3,elements4,mark3,mark4);
        [coordinates,elements4,irregular,boundary] = ...
            QrefineR(coordinates,elements4,irregular,boundary,marked);
        [coordinates,elements4,elements3] = ...
            regularizeedges_tri(coordinates,elements4,irregular);
        % Stopping criterion
        if isempty(marked)|| (size(coordinates,1)>N_max)
            break
        end
    end
    clf
    patch('Faces',elements4,'Vertices',coordinates,'FaceColor','none')
    hold on
    patch('Faces',elements3,'Vertices',coordinates,'FaceColor','none')
    axis off
    print(figure(1),['QRG\MeshQrefineRG' num2str(i)],'-depsc2')
end

%% *** Do QrefineRB
coordinates = c; elements = e4; boundary = b;  pic = video{1}; pic = pic(m:-1:1,:);
if ~exist('QRB', 'dir')
    mkdir('QRB')
end
it = 0;
while 1
    it = it + 1;
    fprintf('******** QrefineRB ******** STEP  %d **********\n',it);
    %*** compute refinement indicator
    etaR = computeEtaR_pict(pic,elements,coordinates);
        if etaR <= tol^2
        break
    end
    %*** mark elements
    marked = markElementsDoerfler_pict(etaR,0.5);
    %*** refine mesh
    [coordinates,elements,boundary]= QrefineRB(coordinates,elements,boundary,marked);
    % Stopping criterion
         if isempty(marked)|| (size(coordinates,1)>N_max)
            break
        end

end
clf
patch('Faces',elements,'Vertices',coordinates,'FaceColor','none')
axis off
print(figure(1),['QRB\MeshQrefineRB' num2str(1)],'-depsc2')

for i=2:int8(numberofframes)
    pic = video{i};
    pic = pic(m:-1:1,:);
    it = 2*it;
    while 1 
        it = it -1;
        fprintf('******** QcoarsenRB ******** STEP  %d **********\n',it);
        %*** Mark elements for coarsening
        marked = 1:size(elements,1);
        %*** Try to coarse mesh
        [coordinates,elements,boundary] = ...
            QcoarsenRB(N0,coordinates,elements,boundary,marked);
        if size(coordinates,1)<N_min || it == 0
            break
        end

    end
it = 0;
    while 1
        it = it + 1;
        fprintf('******** QrefineRB ******** STEP  %d **********\n',it);
        %*** compute refinement indicator
        etaR = computeEtaR_pict(pic,elements,coordinates);
        %*** mark elements
        marked = markElementsDoerfler_pict(etaR,0.5);
        %*** refine mesh
        [coordinates,elements,boundary]= QrefineRB(coordinates,elements,boundary,marked);
            % Stopping criterion
         if isempty(marked)|| (size(coordinates,1)>N_max)
            break
         end     
    end
    clf
    patch('Faces',elements,'Vertices',coordinates,'FaceColor','none')
    axis off
    print(figure(1),['QRB\MeshQrefineRB' num2str(i)],'-depsc2')
end

%% *** Do TrefineR
coordinates = c; elements = e3; irregular = zeros(0,3); boundary = b;  pic = video{1}; pic = pic(m:-1:1,:);
if ~exist('TR', 'dir')
    mkdir('TR')
end
it = 0;
while 1
    it = it + 1;
    fprintf('******** TrefineR ******** STEP  %d **********\n',it);
    %*** compute refinement indicator
    etaR = computeEtaR_pict(pic,elements,coordinates);
    if etaR <= tol^2
        break
    end
    %*** mark elements
    marked = markElementsDoerfler_pict(etaR,0.5);
    %*** refine mesh
    [coordinates,elements,irregular,boundary]= TrefineR(coordinates,elements,irregular,boundary,marked);
    % Stopping criterion
    if isempty(marked)|| (size(coordinates,1)>N_max)
        break
    end
end
clf
patch('Faces',elements,'Vertices',coordinates,'FaceColor','none')
axis off
print(figure(1),['TR\MeshTrefineR' num2str(1)],'-depsc2')

for i=2:int8(numberofframes)
    pic = video{i};
    pic = pic(m:-1:1,:);
    it = 2* it;
    while 1
        it = it - 1;
        fprintf('******** TcoarsenR ******** STEP  %d **********\n',it);
        %*** Mark elements for coarsening
        marked = 1:size(elements,1);
        %*** Try to coarse mesh
        [coordinates,elements,irregular,boundary] = ...
            TcoarsenR(N0,coordinates,elements,irregular,boundary,marked);
        if size(coordinates,1)<N_min || it == 0
            break
        end

    end
    it = 0;
    while 1
        it = it + 1;
        fprintf('******** TrefineR ******** STEP  %d **********\n',it);
        %*** compute refinement indicator
        etaR = computeEtaR_pict(pic,elements,coordinates);
        if etaR <= tol^2
            break
        end
        %*** mark elements
        marked = markElementsDoerfler_pict(etaR,0.5);
        %*** refine mesh
        [coordinates,elements,irregular,boundary]= TrefineR(coordinates,elements,irregular,boundary,marked);
        % Stopping criterion
        if isempty(marked)|| (size(coordinates,1)>N_max)
            break
        end
    end
    clf
    patch('Faces',elements,'Vertices',coordinates,'FaceColor','none')
    axis off
    print(figure(1),['TR\MeshTrefineR' num2str(i)],'-depsc2')
end

%% *** Do QrefineR
coordinates = c; elements = e4; irregular = zeros(0,3); boundary = b;  pic = video{1}; pic = pic(m:-1:1,:);
if ~exist('QR', 'dir')
    mkdir('QR')
end
it = 0;
while 1
    it = it + 1;
    fprintf('******** QrefineR ******** STEP  %d **********\n',it);
    %*** compute refinement indicator
    etaR = computeEtaR_pict(pic,elements,coordinates);
    if etaR <= tol^2
        break
    end
    %*** mark elements
    marked = markElementsDoerfler_pict(etaR,0.5);
    %*** refine mesh
    [coordinates,elements,irregular,boundary]= QrefineR(coordinates,elements,irregular,boundary,marked);
    % Stopping criterion
    if isempty(marked)|| (size(coordinates,1)>N_max)
        break
    end
end
clf
patch('Faces',elements,'Vertices',coordinates,'FaceColor','none')
axis off
print(figure(1),['QR\MeshQrefineR' num2str(1)],'-depsc2')

for i=2:int8(numberofframes)
    pic = video{i};
    pic = pic(m:-1:1,:);
    it = 2*it;
    % coarsen
    while 1
        it = it -1;
        fprintf('******** QcoarsenR ******** STEP  %d **********\n',it);
        %*** Mark elements for coarsening
        marked = 1:size(elements,1);
        %*** Try to coarse mesh
        [coordinates,elements,irregular,boundary] = ...
            QcoarsenR(N0,coordinates,elements,irregular,boundary,marked);
        if size(coordinates,1)<N_min || it == 0
            break
        end

    end
    % refine
    it = 0;
    while 1
        it = it + 1;
        fprintf('******** QrefineR ******** STEP  %d **********\n',it);
        %*** compute refinement indicator
        etaR = computeEtaR_pict(pic,elements,coordinates);
        if etaR <= tol^2
            break
        end
        %*** mark elements
        marked = markElementsDoerfler_pict(etaR,0.5);
        %*** refine mesh
        [coordinates,elements,irregular,boundary]= QrefineR(coordinates,elements,irregular,boundary,marked);
        % Stopping criterion
        if isempty(marked)|| (size(coordinates,1)>N_max)
            break
        end
    end
    clf
    patch('Faces',elements,'Vertices',coordinates,'FaceColor','none')
    axis off
    print(figure(1),['QR\MeshQrefineR' num2str(i)],'-depsc2')
end
