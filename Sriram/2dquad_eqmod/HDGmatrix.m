function [A] = HDGmatrix(HDG)
% Dec 11, 2013
% Compute the HDG Jacobian matrix

SkelAll2Interior = HDG.SkelAll2Interior;
Nifaces = HDG.Nifaces;
Nfp = HDG.Nfp;
refel = HDG.refel;
m = HDG.m;
taur = HDG.taur;
nx = HDG.nx;
ny = HDG.ny;
K = HDG.K;
Nfaces = HDG.Nfaces;
%LIFT=HDG.LIFT;

maxnnzeros = Nfaces * 2 * Nfp * Nfp * Nifaces;
II = zeros(maxnnzeros,1);
JJ = zeros(maxnnzeros,1);
SS = zeros(maxnnzeros,1);

Nfp2 = Nfp * Nfp;
nnzeros = 0;
I = eye(Nfp);

for e = 1:K
  % e1 solution
  [du_dlam,dqx_dlam,dqy_dlam] = DifflocalSolver(HDG, e);
  for f = 1:Nfaces
    index = (f-1)*Nfp+1:f*Nfp;  
    idxf = m.get_skeletal_face_indices(refel, e, f);
    idxf_interior = SkelAll2Interior(idxf);

    % if this is boundary face, ignore
    if idxf_interior(1) < eps,
         continue;
    end

    %lam = lamAll(idxf);
    
    Jf = m.geometric_factors_face(refel,e,f);
    
    idxv = m.get_discontinuous_face_indices(refel, 1, f);
    
    
% $$$     % construct the residual
% $$$     Fhat = Jf .* (refel.Mr * (qx(idxv) * nx(f) + qy(idxv) * ny(f) + ...
% $$$                               taur * (u(idxv) - lam)));
    du_dlam_f  = du_dlam(idxv, index);
    dqx_dlam_f = dqx_dlam(idxv,index);
    dqy_dlam_f = dqy_dlam(idxv,index);

    dFhat_dlam = refel.Mr * (diag(Jf .* nx(f)) * dqx_dlam_f +...
                             diag(Jf .* ny(f)) * dqy_dlam_f +...
                             diag(Jf .* taur)  * (du_dlam_f - I));
                      
    %    res(SkelAll2Interior(idxf)) = res(SkelAll2Interior(idxf)) + Fhat;
    % self derivative
    II(nnzeros+1:nnzeros+Nfp2) = repmat(idxf_interior, 1, Nfp);
    JJ(nnzeros+1:nnzeros+Nfp2) = repmat(idxf_interior', Nfp, 1);
    SS(nnzeros+1:nnzeros+Nfp2) = dFhat_dlam;
    
    nnzeros = nnzeros + Nfp2;

    for fn = 1:Nfaces
      indexn = (fn-1)*Nfp+1:fn*Nfp;
      idxfn = m.get_skeletal_face_indices(refel, e, fn);
      idxf_interiorn = SkelAll2Interior(idxfn);
     
      
      % derivative wrt f is already counted in the identity I above
      if (idxf_interiorn(1) < eps) || (fn == f)
        continue;
      end
      
      
      du_dlam_fn  = du_dlam(idxv, indexn);
      dqx_dlam_fn = dqx_dlam(idxv,indexn);
      dqy_dlam_fn = dqy_dlam(idxv,indexn);

      dFhat_dlamn = refel.Mr * (diag(Jf .* nx(f)) * dqx_dlam_fn +...
                                diag(Jf .* ny(f)) * dqy_dlam_fn +...
                                diag(Jf .* taur)  * du_dlam_fn );

      % neighbor derivative
      II(nnzeros+1:nnzeros+Nfp2) = repmat(idxf_interior, 1, Nfp);
      JJ(nnzeros+1:nnzeros+Nfp2) = repmat(idxf_interiorn', Nfp, 1);
      SS(nnzeros+1:nnzeros+Nfp2) = dFhat_dlamn;
      
      nnzeros = nnzeros + Nfp2;
    end

  end
end

A = sparse(II(1:nnzeros),JJ(1:nnzeros),SS(1:nnzeros),...
           Nifaces * Nfp,Nifaces * Nfp);
