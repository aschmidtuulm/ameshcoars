# ameshcoars
## Efficient Matlab Implementation of Adaptive Mesh Coarsening in 2D

This project provides an efficient implementation of various adaptive mesh coarsening strategies in two dimensions. We coarsen meshes that are adaptively generated using the following refinement strategies. 

![alt text](https://github.com/aschmidtuulm/ameshref/blob/master/RefinementmethodsAMESHREF.png)

## Getting Started

For a use of the coarsening methods download the complete repository _/ameshcoars_ with the test examples and run them on your computer with Matlab. To use it within your own examples you need the repository _/refinement_ and _/coarsening_. Please see the following instructions for a correct use.

### Prerequisites

MATLAB

## Running the test examples

/example1: run _test_moving_circle.m_ (refinement along a moving circle)

/example2: run _example_P1_irr.m_, _exampleQ1_irr,_ and _exampleP1Q1.m_ in the corresponding repositories. Vary the refinement method in _exampleP1Q1.m_ (mesh refinement in the context of the adaptive finite element method for a quasi-stationary partial differential equation)

/example3: run _example3.m_ (triangulation of a GIF)

/example4: run _example4.m_ (local coarsening of a uniformly refined triangulation)

### Deployment for your own examples -  Data structure

To be able to deploy this package within your framework you need the following data strucure of the mesh:

For triangles you need to define _coordinates_, _elements3_, and optionally boundary data _dirichlet_ or _neumann_. Note, that for coarsening of TrefineR it is important to define the initial triangulation such that the smallest index is stored within an element at position one. In general, we number elements in a counterclockwise order.

![alt text](https://github.com/aschmidtuulm/ameshref/blob/master/TriangulationWithTriangles.png)

Similarly, for quadrilaterals you need to define _coordinates, elements4_, and optionally boundary data _dirichlet_ or _neumann_. 

![alt text](https://github.com/aschmidtuulm/ameshref/blob/master/TriangulationWithQuadrilaterals.png)

For meshes with hanging nodes an additional data vector named irregular is needed, where _irregular(l,1)_ and _irregular(l,2_) are the starting and end point of the lth-irregular edge with hanging node stored in _irregular(l,3)_. If there are no irregular edges, a predefinition 
_irregular = zeros(0,3)_ is needed. Also note to clear the variables _nG_ and _nB_ before you start running your code as they are persistent variables!

### How to call the functions

We abbreviate the data structure as _coordinates (C), elements3 (E3), elements4 (E4), irregular (I), dirichlet (D)_, and _neumann (N)_. Furthermore, marked elements are stored with the corresponding index in the variable _marked_. _N0_ is the number of coordinates of the initial mesh.

The different mesh coarsening methods are then called by (remember that _D_ and _N_ are optional arguments)

#### TcoarsenR

```
[C,E3,I,D,N] = TcoarsenR(N0,C,E3,I,D,N,marked)
```
#### QcoarsenR

```
[C,E4,I,D,N] = QcoarsenR(N0,C,E4,I,D,N,marked)
```

#### TcoarsenNVB

```
[C,E3,D,N] = TcoarsenNVB(N0,C,E3,D,N,marked)
```
#### TcoarsenRGB

```
[C,E3,D,N] = TcoarsenRGB(N0,C,E3,D,N,marked)
```
#### TcoarsenRG

```
[C,E3,D,N] = TcoarsenRG(N0,C,E3,D,N,marked)
```
#### QcoarsenRB

```
[C,E4,D,N] = QcoarsenRB(N0,C,E4,D,N,marked)
```
#### QcoarsenRG

```
[C,E3,E4,D,N] = Qcoarsen(N0,C,E3,E4,D,N,marked3,marked4)
```

#### Minimal example
The used functions are provided in the repository _/example1_ and _/refinement_:
```
addpath('../refinement')
addpath('../coarsening')

%% Define initial mesh
coordinates = [0,0;1,0;1,1;0,1;2,0;2,1];
elements = [3,1,2;1,3,4;2,6,3;6,2,5];
boundary = [1,2;2,5;5,6;6,3;3,4;4,1];
N0 = size(coordinates,1);
c_old = 0;

%% Refine uniformly
while 1
    marked = 1:size(elements,1);
    [coordinates,elements,boundary] ...
        = TrefineRGB(coordinates,elements,boundary,marked);
    if isempty(marked)|| (size(coordinates,1)>1e3)
        break
    end
end

% define discrete points (here a disc)
phi = -pi:pi/50:pi;
r = 0.2:1/50:0.4;
[r,phi] = meshgrid(r,phi);
s = r.*cos(phi);
t = r.*sin(phi);
points = [s(:)+1,t(:)+0.5];

%% Coarsen at discrete points
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
```
#### Plot your mesh with Matlab
triangular mesh:
```
patch('Faces',elements3,'Vertices',coordinates,'Facecolor','none')
```
quadrilateral mesh:
```
patch('Faces',elements4,'Vertices',coordinates,'Facecolor','none')
```

## Authors

* **Stefan A. Funken** - **Anja Schmidt** Institute for Numerical Mathematics, Ulm University, Germany

If you use _TcoarsenRGB_ in scientific work, please cite:

* Stefan A. Funken, and Anja Schmidt. **"A Coarsening Algorithm on Adaptive Red-Green-Blue Refined Meshes"**, 2020

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
