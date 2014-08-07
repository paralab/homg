function [du_dlam, dqx_dlam, dqy_dlam] = DifflocalSolver(HDG, e1)
% $$$ function [u,qx,qy] = DifflocalSolver(HDG, e1, lam, forcing)
% Tan Bui Dec 11, 2013
% compute derivatives for the local solver

m     = HDG.m;
refel = HDG.refel;
taur  = HDG.taur;
LIFT  = HDG.LIFT;
VtoF  = HDG.VtoF;
Nfp   = HDG.Nfp;

% predefined normal vector, don't like it but stick with it for now
nx = [-1, 1, 0, 0];
ny = [0, 0, -1, 1];

% the number of volume point
Nv = refel.Nrp ^ (refel.dim);

uu = zeros(Nv,Nv);
uqx = zeros(Nv,Nv);
uqy = zeros(Nv,Nv);

% $$$ rhsqx = zeros(Nv,1);
% $$$ rhsqy = zeros(Nv,1);
% $$$ rhsu  = zeros(Nv,1);
      
pts = m.element_nodes(e1, refel);
[Jv, Dv] = m.geometric_factors(refel, pts);

eMat = m.element_mass(e1, refel, Jv);
eMatInv = inv(eMat);

% number of faces
Nfaces = refel.dim * 2;

du_dlam  = zeros(Nv, Nfp * Nfaces);
dqx_dlam = zeros(Nv, Nfp * Nfaces);
dqy_dlam = zeros(Nv, Nfp * Nfaces);

drhsqx_dlam = zeros(Nv, Nfp * Nfaces);
drhsqy_dlam = zeros(Nv, Nfp * Nfaces);
drhu_dlam   = zeros(Nv, Nfp * Nfaces);

dF_dlam     = zeros(Nv, Nfp * Nfaces);

% $$$ % compute the forcing
% $$$ rhsu = eMat * forcing(pts);

% advection stiffness
[Kex,Key] = m.element_stiffness_advection(e1, refel, Jv, Dv);

uqx = Kex;
uqy = Key;
% residual for qx and qy equations
for f = 1:Nfaces %
  index = (f-1)*Nfp+1:f*Nfp;
%  idxf = m.get_skeletal_face_indices(refel, e1, f);      
%  lamlocal = lam(idxf);
  
  % geometrix factors at gll points on face
  Jf = m.geometric_factors_face(refel,e1,f);
  
%  idxv = m.get_discontinuous_face_indices(refel, 1, f);       
  
  %% residual due to lambda
  % rhsfx = Jf .* (refel.Mr * lamlocal) * nx(f);
  %% lift to volume residual q equation
  % rhsqx(idxv) = rhsqx(idxv) + rhsfx;
  drhsqx_dlam(:,index) = LIFT(:,:,f) * (diag(Jf .* nx(f)));

% $$$   %-------- Testing ------------
% $$$   idxf = m.get_skeletal_face_indices(refel, e1, f);      
% $$$   lamlocal = rand(length(idxf),1);;
% $$$   idxv = m.get_discontinuous_face_indices(refel, 1, f);
% $$$   
% $$$   rhsfx = Jf .* (refel.Mr * lamlocal) * nx(f);
% $$$   
% $$$   lamRand = rand(size(lamlocal));
% $$$   lamlocal = lamlocal + i * eps * lamRand;
% $$$   rhsfx = Jf .* (refel.Mr * lamlocal) * nx(f);
% $$$   
% $$$   norm(drhsqx_dlam(idxv,(f-1)*Nfp+1:f*Nfp) * lamRand - imag(rhsfx)/eps)
% $$$   %-------- End testing --------
  
  % rhsfy = Jf .* (refel.Mr * lamlocal) * ny(f);
  % rhsqy(idxv) = rhsqy(idxv) + rhsfy;
  drhsqy_dlam(:,index) = LIFT(:,:,f) * (diag(Jf .* ny(f)));
  
  % lift to volume residual u equation
  % rhsu(idxv)  = rhsu(idxv) - ...
  %    taur * Jf .* (refel.Mr * lamlocal);
  drhsu_dlam(:,index) = -LIFT(:,:,f) * (diag(Jf .* taur));
  
  % lift to volume for uu
  bdry =  LIFT(:,:,f) * (diag(Jf) * VtoF(:,:,f));
  bdryx = LIFT(:,:,f) * (diag(Jf) * VtoF(:,:,f)) * nx(f);
  bdryy = LIFT(:,:,f) * (diag(Jf) * VtoF(:,:,f)) * ny(f);
  
  uu  = uu  - taur * bdry; 
  uqx = uqx -        bdryx;
  uqy = uqy -        bdryy;
end

qxMatrix = -eMatInv * Kex;
% qxrhs = eMatInv * rhsqx;
drhsqx_dlam = eMatInv * drhsqx_dlam;

qyMatrix = -eMatInv * Key;
%qyrhs = eMatInv * rhsqy;
drhsqy_dlam = eMatInv * drhsqy_dlam;

%F = rhsu - uqx * qxrhs - uqy * qyrhs;
dF_dlam = drhsu_dlam - uqx * drhsqx_dlam - uqy * drhsqy_dlam;

dF = uqx * qxMatrix + uqy * qyMatrix + uu;

%u = dF \ F;
% $$$ [L,U] = lu(dF);
% $$$ du_dlam = U\(L \ dF_dlam);
du_dlam = inv(dF) * dF_dlam;

%qx = qxMatrix * u + qxrhs;
%qy = qyMatrix * u + qyrhs;
dqx_dlam = qxMatrix * du_dlam + drhsqx_dlam;
dqy_dlam = qyMatrix * du_dlam + drhsqy_dlam;