%example and visualizes the P1-FEM solution for a quasi-stationary
%Laplace equation via an adaptive algorithm with refinement and coarsening
%
%    example solves Laplace equation
%      - div(grad(u)) = exp(-10*||x-x_0(t)||^2)  in Omega x (0,T)
%                   u = 0                        on the Dirichlet boundary
%    with x_0(t) = [1.5 + cos(t),1.5 + sin(t)] on a geometry described by
%    triangles.
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
%    TrefineR. The sub-loop is aborted if the error estimator is
%    smaller than the given tolerance.
%
%    Coarsening sub-loop: For the first run, we use the indicators from the
%    final run of the refinement sub-loop. For the other runs, we solve for
%    the actual mesh and compute new error indicators. We mark an element
%    for coarsening if its indicator is smaller than
%    sigma*tol/size(elements,1) and coarsen the mesh via TcoarsenR.
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
T = 2*pi;     % final time
N = 50;        % number of time steps
rho =  0.25;     % parameter for refinement
sigma = 0.25;   % parameter for coarsening
tol = 0.03;     % tolerance
nEmax = 5000;

%*** Initialization
load coordinates.dat
load elements.dat
load dirichlet.dat
neumann = [];
irregular = zeros(0,3);


N0 = size(elements,1); % vertices that may not be coarsened

%*** Time loop
for t = linspace(0,T,N+1)
    
    %*** Refinement loop
    k = 0;
    while 1
        k = k + 1;
        %*** Compute discrete solution
        [x,~] = solveLaplace(coordinates,elements,irregular,dirichlet,neumann,@f,@g,@uD);
        %*** Compute error indicators
        eta = computeEtaR(x,coordinates,elements,irregular,dirichlet,neumann,@f,@g);
        %*** Stop if solution is sufficiently accurate
        if sum(eta) <= tol^2
            break
        end
        % *** Stopping criteria
        if size(elements,1) >= nEmax
            break
        end
        %*** Otherwise, mark elements for refinement via Doerfler criterion
        [indicators,idx] = sort(eta,'descend');
        sumeta = cumsum(indicators);
        ell = find(sumeta>=sumeta(end)*rho,1);
        marked = idx(1:ell);
        %*** Refine mesh by R
        [coordinates,elements,irregular,dirichlet,neumann] ...
            = TrefineR(coordinates,elements,irregular,dirichlet,neumann,marked);
    end
    
    %*** Visualization
    figure(1)
    patch('Faces',elements,'Vertices',coordinates,'FaceVertexCData',x,'Facecolor','interp')
    title(sprintf('# Elements after Refinement = %s',int2str(size(elements,1))),'FontSize',20);
    axis([-0.1 3.1 -0.1 3.1])
    axis equal
    axis off
    view(2)
    pause(1)
    
    %***Coarsening loop
    while k
        k = k-1;
        nE = size(elements,1);
        %*** Mark elements for coarsening
        marked = find(eta < sigma* tol^2/nE);
        [coordinates,elements,irregular,dirichlet,neumann] = ...
            TcoarsenR(N0,coordinates,elements,irregular,dirichlet,neumann,marked);
        %*** Compute discrete solution
        [x,~] = solveLaplace(coordinates,elements,irregular,dirichlet,neumann,@f,@g,@uD);
        %*** Compute error indicators
        eta = computeEtaR(x,coordinates,elements,irregular,dirichlet,neumann,@f,@g);
        %*** Stop if mesh cannot be coarsened anymore
        if size(elements,1) == nE || k == 0
            break
        end
    end
    %*** Visualization
    figure(1)
    patch('Faces',elements,'Vertices',coordinates,'FaceVertexCData',x,'Facecolor','interp')
    title(sprintf('# Elements after Coarsening = %s',int2str(size(elements,1))),'FontSize',20);
    axis([-0.1 3.1 -0.1 3.1])
    axis equal
    axis off
    view(2)
    pause(1)
end