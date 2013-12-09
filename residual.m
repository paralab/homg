function [res] = residual(lamInterior,HDG,forcing, Bdata)
% Dec 9, 2013
% Compute the HDG residual

SkelInterior2All = HDG.SkelInterior2All;
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

res = zeros(Nfp * Nifaces, 1);

% lamAll
lamAll = zeros(Nsfaces * Nfp,1);
lamAll(Bmaps) = Bdata;
lamAll(SkelInterior2All) = lamInterior;

for iface = 1:Nifaces
  index = (iface-1)*Nfp+1:iface*Nfp;
  lam = lamInterior(index);
  
  % global skeleton face
  sf = InteriorF2AllF(iface);
  
  [eL, fL, eR, fR]  =m.get_face_elements(sf);

  JfL = m.geometric_factors_face(refel,eL,fL);
  JfR = m.geometric_factors_face(refel,eR,fR);
  
  idxvL = m.get_discontinuous_face_indices(refel, 1, fL);
  idxvR = m.get_discontinuous_face_indices(refel, 1, fR);
  
  % e1 solution
  [uL,qxL,qyL] = localSolver(HDG, eL, lamAll, forcing);
  
  % e2 solution
  [uR,qxR,qyR] = localSolver(HDG, eR, lamAll , forcing);

  % construct the residual
  FhatL = JfL .* (refel.Mr * (qxL(idxvL) * nx(fL) + qyL(idxvL) * ny(fL) + ...
                             taur * (uL(idxvL) - lam)));
  
  FhatR = JfR .* (refel.Mr * (qxR(idxvR) * nx(fR) + qyR(idxvR) * ny(fR) + ...
                             taur * (uR(idxvR) - lam)));

  res(index) = FhatL + FhatR;
  keyboard
end