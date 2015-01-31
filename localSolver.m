function [u,qx,qy] = localSolver(HDG, e1, lam, forcing)
% compute local solution from the trace unknowns

m     = HDG.m;
refel = HDG.refel;
taur  = HDG.taur;
LIFT  = HDG.LIFT;
VtoF  = HDG.VtoF;

% predefined normal vector, don't like it but stick with it for now
nx = [-1, 1, 0, 0];
ny = [0, 0, -1, 1];


% the number of volume point
Nv = refel.Nrp ^ (refel.dim);

uu = zeros(Nv,Nv);
uqx = zeros(Nv,Nv);
uqy = zeros(Nv,Nv);

rhsqx = zeros(Nv,1);
rhsqy = zeros(Nv,1);
rhsu  = zeros(Nv,1);
      
pts = m.element_nodes(e1, refel);
[Jv, Dv] = m.geometric_factors(refel, pts);

eMat = m.element_mass(e1, refel, Jv);
eMatInv = inv(eMat);

% number of faces
Nfaces = refel.dim * 2;

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
  rhsu(idxv)  = rhsu(idxv) + ...
      taur * Jf .* (refel.Mr * lam(idxf));  
  
  % lift to volume for uu
  bdry =  LIFT(:,:,f) * (diag(Jf) * VtoF(:,:,f));
  bdryx = LIFT(:,:,f) * (diag(Jf) * VtoF(:,:,f)) * nx(f);
  bdryy = LIFT(:,:,f) * (diag(Jf) * VtoF(:,:,f)) * ny(f);
  
  uu  = uu  + taur * bdry;     
  uqx = uqx -        bdryx;
  uqy = uqy -        bdryy;
end



uqx=-uqx;  
uqy=-uqy;  

qxMatrix = eMatInv * Kex;  
qxrhs = -eMatInv * rhsqx;  

qyMatrix = eMatInv * Key;  
qyrhs = -eMatInv * rhsqy;  

F = rhsu - uqx * qxrhs - uqy * qyrhs;
dF = uqx * qxMatrix + uqy * qyMatrix + uu;

u = dF \ F;
qx = qxMatrix * u + qxrhs;
qy = qyMatrix * u + qyrhs;
