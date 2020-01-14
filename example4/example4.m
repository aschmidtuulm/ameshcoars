%example4: coarsen a fine mesh locally (through a discrete point set)
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

close all, clear all, clc

addpath('../refinement')
addpath('../coarsening')

%% Define initial mesh
c = [0,0;1,0;1,1;0,1;2,0;2,1];
e3 = [3,1,2;1,3,4;2,6,3;6,2,5];
e4 = [1,2,3,4;2,5,6,3];
b = [1,2;2,5;5,6;6,3;3,4;4,1];
N0 = size(c,1);

% Set parameters
N_max = 1e3;

% define discrete points (here a disc)
phi = -pi:pi/50:pi;
r = 0.2:1/50:0.4;
[r,phi] = meshgrid(r,phi);
s = r.*cos(phi);
t = r.*sin(phi);
points = [s(:)+1,t(:)+0.5];


%% *** TrefineRGB
coordinates = c; elements = e3; boundary = b; c_old = 0;
% Refine uniformly
while 1
    marked = 1:size(elements,1);
    [coordinates,elements,boundary] ...
        = TrefineRGB(coordinates,elements,boundary,marked);
    if isempty(marked)|| (size(coordinates,1)>N_max)
        break
    end
end

% Coarsen at discrete points
while 1
    mark3 = point2element(coordinates,elements,points);
    [coordinates,elements,boundary] = TcoarsenRGB(N0,coordinates,elements,boundary,mark3);
    if size(coordinates,1) == c_old
        break
    end
    c_old = size(coordinates,1);
end

% plot results
clf
patch('Faces',elements,'Vertices',coordinates,'Facecolor','none')
axis equal
axis off
pause(1)

%% *** TrefineNVB
coordinates = c; elements = e3; boundary = b; c_old = 0;
% Refine uniformly
while 1
    marked = 1:size(elements,1);
    [coordinates,elements,boundary] ...
        = TrefineNVB(coordinates,elements,boundary,marked);
    if isempty(marked)|| (size(coordinates,1)>N_max)
        break
    end
end

% Coarsen at discrete points
while 1
    mark3 = point2element(coordinates,elements,points);
    [coordinates,elements,boundary] = TcoarsenNVB(N0,coordinates,elements,boundary,mark3);
    if size(coordinates,1) == c_old
        break
    end
    c_old = size(coordinates,1);
end

% plot results
clf
patch('Faces',elements,'Vertices',coordinates,'Facecolor','none')
axis equal
axis off
pause(1)

%% *** TrefineRG
coordinates = c; elements = e3; boundary = b; c_old = 0;
% Refine uniformly
while 1
    marked = 1:size(elements,1);
    [coordinates,elements,boundary] ...
        = TrefineRG(coordinates,elements,boundary,marked);
    if isempty(marked)|| (size(coordinates,1)>N_max)
        break
    end
end

% Coarsen at discrete points
while 1
    mark3 = point2element(coordinates,elements,points);
    [coordinates,elements,boundary] = TcoarsenRG(N0,coordinates,elements,boundary,mark3);
    if size(coordinates,1) == c_old
        break
    end
    c_old = size(coordinates,1);
end

% plot results
clf
patch('Faces',elements,'Vertices',coordinates,'Facecolor','none')
axis equal
axis off
pause(1)

%% *** TrefineR
coordinates = c; elements = e3; irregular = zeros(0,3); boundary = b; c_old = 0;
% Refine uniformly
while 1
    marked = 1:size(elements,1);
    [coordinates,elements,irregular,boundary] ...
        = TrefineR(coordinates,elements,irregular,boundary,marked);
    if isempty(marked)|| (size(coordinates,1)>N_max)
        break
    end
end

% Coarsen at discrete points
while 1
    mark3 = point2element(coordinates,elements,points);
    [coordinates,elements,irregular,boundary] = ...
        TcoarsenR(N0,coordinates,elements,irregular,boundary,mark3);
    if size(coordinates,1) == c_old
        break
    end
    c_old = size(coordinates,1);
end

% plot results
clf
patch('Faces',elements,'Vertices',coordinates,'Facecolor','none')
axis equal
axis off
pause(1)

%% *** QrefineRB
coordinates = c; elements = e4; boundary = b; c_old = 0;
% Refine uniformly
while 1
    marked = 1:size(elements,1);
    [coordinates,elements,boundary] ...
        = QrefineRB(coordinates,elements,boundary,marked);
    if isempty(marked)|| (size(coordinates,1)>N_max)
        break
    end
end

% Coarsen at discrete points
while 1
    mark4 = unique([point2element(coordinates,elements(:,1:3),points);point2element(coordinates,elements(:,[1,3,4]),points)]);
    [coordinates,elements,boundary] = QcoarsenRB(N0,coordinates,elements,boundary,mark4);
    if size(coordinates,1) == c_old
        break
    end
    c_old = size(coordinates,1);
end

% plot results
clf
patch('Faces',elements,'Vertices',coordinates,'Facecolor','none')
axis equal
axis off
pause(1)

%% *** QrefineRG
coordinates = c; elements3 = zeros(0,3); elements4 = e4; boundary = b; c_old = 0;
% Refine uniformly
while 1
    marked3 = 1:size(elements3,1);
    marked4 = 1:size(elements4,1);
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

% Coarsen at discrete points
while 1
    mark4 = unique([point2element(coordinates,elements4(:,1:3),points);point2element(coordinates,elements4(:,[1,3,4]),points)]);
    [coordinates,elements3,elements4,boundary] = ...
        QcoarsenRG(N0,coordinates,elements3,elements4,boundary,mark4);
    if size(coordinates,1) == c_old
        break
    end
    c_old = size(coordinates,1);
end

% plot results
 clf
    patch('Faces',elements4,'Vertices',coordinates,'Facecolor','none')
    hold on
    patch('Faces',elements3,'Vertices',coordinates,'Facecolor','none')
    axis equal
    axis off
pause(1)

%% *** QrefineR
coordinates = c; elements = e4; irregular = zeros(0,3); boundary = b; c_old = 0;
% Refine uniformly
while 1
    marked = 1:size(elements,1);
    [coordinates,elements,irregular,boundary] ...
        = QrefineR(coordinates,elements,irregular,boundary,marked);
    if isempty(marked)|| (size(coordinates,1)>N_max)
        break
    end
end

% Coarsen at discrete points
while 1
    mark4 = unique([point2element(coordinates,elements(:,1:3),points);point2element(coordinates,elements(:,[1,3,4]),points)]);
    [coordinates,elements,irregular,boundary] = ...
        QcoarsenR(N0,coordinates,elements,irregular,boundary,mark4);
    if size(coordinates,1) == c_old
        break
    end
    c_old = size(coordinates,1);
end

% plot results
clf
patch('Faces',elements,'Vertices',coordinates,'Facecolor','none')
axis equal
axis off
pause(1)