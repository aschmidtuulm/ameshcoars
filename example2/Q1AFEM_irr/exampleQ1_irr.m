%example and visualizes the Q1-FEM solution for a quasi-stationary
%Laplace equation via an adaptive algorithm with refinement and coarsening
%
%    example solves Laplace equation
%      - div(grad(u)) = exp(-10*||x-x_0(t)||^2)  in Omega x (0,T)
%                   u = 0                        on the Dirichlet boundary
%    with x_0(t) = [1.5 + cos(t),1.5 + sin(t)] on a geometry described by
%    quadrilaterals.
%
%    Parameters(can be set by direct manipulation in this file):
%       T: real; final time. T=2pi corresponds to one loop of the source
%          term
%       N: integer; number of equidistant time steps
%       rho: in (0,1); controls the amount of refined elements per run of
%          refinement loop. The closer to 1 the more elements are refined
%       sigma: in (0,1); controls the amount of coarsened elements per run
%          of coarsening loop. The closer to 1 the more elements are
%          coarsened
%       tol: real; tolerance to control the stop criterion of the
%          refinement and coarsening loops
%
%    The program contains one global time loop in which time is increased
%    in equidistant time steps. The time t essentially enters the source
%    term when the discrete solution is computed. The time loop includes
%    one refinement sub-loop and one coarsening sub-loop.
%
%    Refinement sub-loop: First, the discrete solution and the residual
%    error indicators are computed. Then, we mark the elements for
%    refinement via the Doerfler criterion and refine the mesh by 
%    QrefineR. The sub-loop is aborted if the error estimator is
%    smaller than the given tolerance.
%
%    Coarsening sub-loop: For the first run, we use the indicators from the
%    final run of the refinement sub-loop. For the other runs, we solve for
%    the actual mesh and compute new error indicators. We mark an element
%    for coarsening if its indicator is smaller than
%    sigma*tol/size(elements,1) and coarsen the mesh via QcoarsenR.
%    The coarsening sub-loop is aborted, if no element was coarsened.
%
%Remark:
%
%    This program is a supplement to the paper 
%    >> CRITERIA FOR NON-RECURSIVE LOCAL COARSENING IN 2D <<
%    by S. Funken, and A. Schmidt. The reader should 
%    consult that paper for more information. 
%
%    based on
%    >> Efficient Implementation of Adaptive P1-FEM in Matlab <<
%    by S. Funken, D. Praetorius, and P. Wissgott.
%
%Authors:
%
%    S. Funken, A. Schmidt  04-11-19
clear all
close all
clc

addpath('../../refinement')
addpath('../../coarsening')

global t

%*** Parameters
T = 2 * pi;     % final time
N = 50;        % number of time steps
rho = 0.25;     % parameter for refinement
sigma = 0.25;   % parameter for coarsening
tol = 0.03;     % tolerance
nEmax = 2000;

%*** Initialization
elements4 = [1,2,6,5;2,3,7,6;3,4,8,7;5,6,10,9;7,8,12,11;9,10,14,13;10,11,15,14;11,12,16,15];
coordinates = [0,0;1,0;2,0;3,0;0,1;1,1;2,1;3,1;0,2;1,2;2,2;3,2;0,3;1,3;2,3;3,3];
dirichlet = [1,2;2,3;3,4;4,8;8,12;12,16;16,15;15,14;14,13;13,9;9,5;5,1;6,7;7,11;11,10;10,6];
neumann = [];
irregular = zeros(0,3);
N0 = size(coordinates,1); % vertices that may not be coarsened

%*** Time loop
for t = linspace(0,T,N+1)
    
    %*** Refinement loop
    k = 0;
    while 1
        k = k + 1;
        %*** Compute discrete solution
        [x,~] = solveLaplace(coordinates,elements4,irregular,dirichlet,neumann,@f,@g,@uD);
        %*** Compute error indicators
        eta = computeEtaR(x,coordinates,elements4,irregular,dirichlet,neumann,@f,@g);
        %*** Stop if solution is sufficiently accurate
        if sum(eta) <= tol^2
            break
        end
        %*** Stopping criteria
        if size(elements4,1) >= nEmax
            break
        end
        %*** Otherwise, mark elements for refinement via Doerfler criterion
        [indicators,idx] = sort(eta,'descend');
        sumeta = cumsum(indicators);
        ell = find(sumeta>=sumeta(end)*rho,1);
        marked = idx(1:ell);
        
        
        %*** Refine mesh by R
        [coordinates,elements4,irregular,dirichlet,neumann] ...
            = QrefineR(coordinates,elements4,irregular,dirichlet,neumann,marked);
    end
    
    %*** Visualization
    figure(1)
    clf
    patch('Faces',elements4,'Vertices',coordinates,'FaceVertexCData',x,'Facecolor','interp');
    title(sprintf('# Elements = %s',int2str(size(elements4,1))),'FontSize',20);
    axis([-0.1 3.1 -0.1 3.1])
    axis equal
    axis off
    view(2)
    pause(1)
    
    %*** Coarsening loop
    while k
        k = k-1;
        nE = size(elements4,1);
        %*** Mark elements for coarsening
        marked = find(eta < sigma* tol^2/nE);
        [coordinates,elements4,irregular,dirichlet,neumann] = ...
            QcoarsenR(N0,coordinates,elements4,irregular,dirichlet,neumann,marked);
        %*** Compute discrete solution
        [x,~] = solveLaplace(coordinates,elements4,irregular,dirichlet,neumann,@f,@g,@uD);
        %*** Compute error indicators
        eta = computeEtaR(x,coordinates,elements4,irregular,dirichlet,neumann,@f,@g);
        %*** Stop if mesh cannot be coarsened anymore
        if size(elements4,1) == nE || k == 0
            break
        end
        
    end
    %*** Visualization
    figure(1)
    clf
    patch('Faces',elements4,'Vertices',coordinates,'FaceVertexCData',x,'Facecolor','interp');
    title(sprintf('# Elements = %s',int2str(size(elements4,1))),'FontSize',20);
    axis([-0.1 3.1 -0.1 3.1])
    axis equal
    axis off
    view(2)
    pause(1)
end