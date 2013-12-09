clear all
% Tan Bui and Hari Sundar, Oct 28, 2013
% Testing HDG method for 2D Laplace equation

% addpath /workspace/tanbui/tanbui/WithHari/homg/

% solution order
order = 6;

% number of elements in x and y directions
nelems = [6,3];

% generate the hexmesh with identity transform for now
m = homg.hexmesh(nelems,@homg.xform.identity);
m.set_order(order);

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

% predefined normal vector, don't like it but stick with it for now
nx = [-1, 1, 0, 0];
ny = [0, 0, -1, 1];

% stabilization parameter
taur = 1;
uu = zeros(Nv,Nv);
uqx = zeros(Nv,Nv);
uqy = zeros(Nv,Nv);

% Also compute uexact for testing
Uexact = zeros(Nv,K);
Qxexact = zeros(Nv,K);
Qyexact = zeros(Nv,K);
for k = 1:K
  pts = m.element_nodes(k, refel);
  Uexact(:,k)  = uexact(pts);
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

% Construct the Lift and VtoF
LIFT = zeros(Nv, Nfp, Nfaces);
VtoF = zeros(Nfp, Nv, Nfaces);
for f = 1:Nfaces
  idxv = m.get_discontinuous_face_indices(refel, 1, f);
  LIFT(idxv,:,f) = refel.Mr;
  for fp = 1:Nfp
    VtoF(fp,idxv(fp),f) = 1;
  end
end

% loop over all faces of the mesh skeleton
for  sf=1:Nsfaces
    [e1, f1, e2, f2]  =m.get_face_elements(sf);

    %% Task 1
    if (e1 > 0) && (e2 > 0), % interior faces
      
      %---- Solve for u, and q in e1-------
      uu(:) = 0; uqx(:) = 0; uqy(:) = 0;
      rhsqx(:) = 0; rhsqy(:) = 0; rhsu(:) = 0;
      pts = m.element_nodes(e1, refel);
      [Jv, Dv] = m.geometric_factors(refel, pts);
      
      eMat = m.element_mass(e1, refel, Jv);
      eMatInv = inv(eMat);

      % compute the forcing
      rhsu = eMat * forcing(pts);

      % advection stiffness
      [Kex,Key] = m.element_stiffness_advection(e1, refel, Jv, Dv);
      
      uqx = Kex;
      uqy = Key;
      % residual for qx and qy equations
      for f = 1:Nfaces %
        idxf = m.get_skeletal_face_indices(refel, e1, f);      
        % geometrix factors at gll points on face
        Jf = m.geometric_factors_face(refel,e1,f);

        idxv = m.get_discontinuous_face_indices(refel, 1, f);       
        
        % residual due to lambda
        rhsfx = Jf .* (refel.Mr * lam(idxf)) * nx(f);
        rhsfy = Jf .* (refel.Mr * lam(idxf)) * ny(f);
        % lift to volume residual q equation
        rhsqx(idxv) = rhsqx(idxv) + rhsfx;
        rhsqy(idxv) = rhsqy(idxv) + rhsfy;
        % lift to volume residual u equation
        rhsu(idxv)  = rhsu(idxv) - ...
            taur * Jf .* (refel.Mr * lam(idxf));
        
        % lift to volume for uu
        bdry =  LIFT(:,:,f) * (diag(Jf) * VtoF(:,:,f));
        bdryx = LIFT(:,:,f) * (diag(Jf) * VtoF(:,:,f)) * nx(f);
        bdryy = LIFT(:,:,f) * (diag(Jf) * VtoF(:,:,f)) * ny(f);

        uu  = uu  - taur * bdry; 
        uqx = uqx -        bdryx;
        uqy = uqy -        bdryy;
      end
      
      qxMatrix = -eMatInv * Kex;
      qxrhs = eMatInv * rhsqx;
      
      qyMatrix = -eMatInv * Key;
      qyrhs = eMatInv * rhsqy;
      
      F = rhsu - uqx * qxrhs - uqy * qyrhs;
      dF = uqx * qxMatrix + uqy * qyMatrix + uu;
      
      u = dF \ F;
      qx = qxMatrix * u + qxrhs;
      qy = qyMatrix * u + qyrhs;
    
    else
      continue;
    end
      
end
