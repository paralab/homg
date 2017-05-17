% Parameters:
%   Bo: Bond number. Ratio of surface tension and gravity, characterizes
%       surface tension in single dimensionless scalar
%   k:  Order of the basis elements for FEM. k > 1 drastically increases
%       computation time.
%   sigma:  Shift parameter. Algorithm will compute eigenvalues closest to
%           the given shift.
%   n_evals:    The number of eigenvalues to compute. The given n_evals
%               will be greater than the number of actual computed
%               eigenvalues since spurious eigenvalues and eigenvectors
%               will be eliminatd.
%   d:  The number of mesh subdivisions. The 3D portion of the container
%       will be divided into a (1+kd) by (1+kd) by (1+kd) cube mesh and 
%       the 2D portion will be divided into a (1+kd) by (1+kd) square mesh.
Bo = 10.0;
k = 1;
sigma = 0.9;
n_evals = 30;
d = 10;

% Containers:
%   Box                             -   'box'
%   Bowl                            -   'bowl'
%   Cylinder                        -   'cylinder'
%   Parametric Cylinder (half)      -   'pcylhalf'
%   Parametric Cylinder (double)    -   'pcyldouble'
container_type = 'bowl';

% Once parameters are set, the problem can be initialized using the number
%   of mesh subdivisions 'd', the order of the basis elements 'k', and
%   container type
[S, SF, MF, Mm] = freesurface.initialize_problem(d, k, container_type);

% With the problem initialized, the eigenvalues and eigenvectors of the
%   system can be solved for using the following method.
[evals, evecs] = freesurface.solve(Bo, sigma, n_evals, S, SF, MF, Mm);