%clear all
% Tan Bui and Hari Sundar, Oct 28, 2013
% Testing HDG method for 2D Laplace equation
% modified for multigrid, 02 Jan 2014

% addpath /workspace/tanbui/tanbui/WithHari/homg/

% solution order
order = 4;

% number of elements in x and y directions
nelems = [16, 16];

% generate mesh heirarchy 
grid = create_hdg_grids(2, @homg.xform.identity, [4 2 1], [8 16]);

% generate the hexmesh with identity transform for now
m = grid.Mesh; % homg.hexmesh(nelems,@homg.xform.identity); 

% reference elements
refel = homg.refel(m.dim, order);

% get total number of faces on the skeleton of the mesh
Nsfaces =m.get_num_faces();

% the number of face points
Nfp = refel.Nrp ^ (refel.dim-1);

% the number of volume point
Nv = refel.Nrp ^ (refel.dim);

% number of faces
Nfaces = refel.dim * 2;

% number of elements
K = prod(m.nelems);

% Initialize lam
lam = zeros(Nfp * Nsfaces,1);

% lambda residual
lamRes = zeros(Nfp * Nsfaces, 1);

% forcing
forcing = @(pts) (sin(2.0 * pi * pts(:,1)) .* sin(pi * pts(:,2)));

% exact solution
% u = 0.5/pi^2 * forcing;
uexact = @(pts) 0.2 / pi^2 * forcing(pts);
qxexact = @(pts) 0.4/pi * (cos(2*pi * pts(:,1)) ...
                           .* sin(pi * pts(:,2)));
qyexact = @(pts) 0.2/pi * (sin(2*pi * pts(:,1)) ...
                           .* cos(pi * pts(:,2)));

% Number of volume unknown for a scalar
rhsqx = zeros(Nv,1);
rhsqy = zeros(Nv,1);
rhsu  = zeros(Nv,1);

uu = zeros(Nv,Nv);
uqx = zeros(Nv,Nv);
uqy = zeros(Nv,Nv);

% Also compute uexact for testing
Uexact  = zeros(Nv,K);
Forcing = zeros(Nv,K);
Qxexact = zeros(Nv,K);
Qyexact = zeros(Nv,K);
for k = 1:K
  pts = m.element_nodes(k, refel);
  Uexact(:,k)  = uexact(pts);
  Forcing(:,k) = forcing(pts);
  Qxexact(:,k) = qxexact(pts);
  Qyexact(:,k) = qyexact(pts);
end

% COMPUTE EXACT LAMBDA FOR TESTING
for  sf=1:Nsfaces
    [e1, f1, e2, f2]  =m.get_face_elements(sf);

    %    if (e1 > 0) && (e2 > 0), % interior faces
    if (f1 > 0)
      pts = m.element_nodes(e1, refel);
      idxf = m.get_skeletal_face_indices(refel, e1, f1);      
      idxv = m.get_discontinuous_face_indices(refel, 1, f1);
      lam(idxf) = uexact(pts(idxv,:));
    end
    if (f2 > 0)
      pts = m.element_nodes(e2, refel);
      idxf = m.get_skeletal_face_indices(refel, e2, f2);      
      idxv = m.get_discontinuous_face_indices(refel, 1, f2);
      lam(idxf) = uexact(pts(idxv,:));
    end
end

Bdata = grid.Mesh.get_boundary_data(grid.refel, Uexact);

u = grid.solve_hdg(10, 'jacobi', 3, 3, Forcing(:), zeros(size(Uexact(:))), Bdata);

%% test errors ... 

L2eu = 0;
L2eqx = 0;
L2eqy = 0;

for e = 1:K
  eu = u(:,e) - Uexact(:,e);
  eqx = qx(:,e) - Qxexact(:,e);
  eqy = qy(:,e) - Qyexact(:,e);
  
  pts = m.element_nodes(e, refel);
  [Jv, Dv] = m.geometric_factors(refel, pts);
  eMat = m.element_mass(e, refel, Jv);
  
  L2eu = L2eu +  eu' * eMat * eu;
  L2eqx = L2eqx +  eqx' * eMat * eqx;
  L2eqy = L2eqy +  eqy' * eMat * eqy;
  
end

fprintf('L2 norm error for u  = %1.15e \n',sqrt(L2eu));
fprintf('L2 norm error for qx = %1.15e \n',sqrt(L2eqx));
fprintf('L2 norm error for qy = %1.15e \n',sqrt(L2eqy));



