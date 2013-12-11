function [res] = residualFast(lamInterior,HDG,forcing, Bdata)
% Dec 10, 2013
% Compute the HDG residual in a faster way

SkelInterior2All = HDG.SkelInterior2All;
SkelAll2Interior = HDG.SkelAll2Interior;
Nifaces = HDG.Nifaces;
Nfp = HDG.Nfp;
Nsfaces = HDG.Nsfaces;
Bmaps = HDG.Bmaps;
InteriorF2AllF = HDG.InteriorF2AllF;
refel = HDG.refel;
m = HDG.m;
LIFT = HDG.LIFT;
taur = HDG.taur;
nx = HDG.nx;
ny = HDG.ny;
K = HDG.K;
Nfaces = HDG.Nfaces;

res = zeros(Nfp * Nifaces, 1);

% lamAll
lamAll = zeros(Nsfaces * Nfp,1);
lamAll(Bmaps) = Bdata;
lamAll(SkelInterior2All) = lamInterior;

for e = 1:K
  % e1 solution
  [u,qx,qy] = localSolver(HDG, e, lamAll, forcing);
  for f = 1:Nfaces
    idxf = m.get_skeletal_face_indices(refel, e, f);
      
    if abs(SkelAll2Interior(idxf(1))) < eps
      continue;
    end
    lam = lamAll(idxf);
    
    Jf = m.geometric_factors_face(refel,e,f);
    idxv = m.get_discontinuous_face_indices(refel, 1, f);
    
    
    % construct the residual
    Fhat = Jf .* (refel.Mr * (qx(idxv) * nx(f) + qy(idxv) * ny(f) + ...
                              taur * (u(idxv) - lam)));
    
    res(SkelAll2Interior(idxf)) = res(SkelAll2Interior(idxf)) + Fhat;
  end
end
