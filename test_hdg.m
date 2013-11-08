refel = homg.refel(2, 6);
m = homg.hexmesh([8, 8], @homg.xform.identity);

NV = (refel.Nrp)^(refel.dim);
NP = (refel.Nrp)^(refel.dim - 1);
NF = m.get_num_faces();

lam = zeros(NF*NP, 1);

for gf=1:NF
    [e1, f1, e2, f2] = m.get_face_elements(gf);

    
    % elem1 
    u1    = rand(NV,1);
    q1x    = rand(NV,1);
    q1y    = rand(NV,1);
    res1x = zeros(NV,1);
    res1y = zeros(NV,1);
        
    if (e1 == -1) 
    
    else
      % not boundary
    
      %% face contributions to the residual
      idx_skel1 = m.get_skeletal_face_indices(refel, e1, 1);
      [J,D] = m.geometric_factors_face(refel, e1, 1);
      
      lam_face_1 = refel.Mr * lam(idx_skel1);
      lam_face_1 = lam_face_1 .* J .* refel.wgll;
    
      
      idx_vol_1 = m.get_discontinuous_face_indices(refel, 1, 1);
      res1x(idx_vol_1) = res1x(idx_vol_1) + lam_face_1;
      
      % do same for faces 2,3,4
      
      %% volume contributions to the residual
      pts    = m.element_nodes(e1, refel);
      [J, D] = m.geometric_factors(refel, pts);
      [Jg, Dg] = m.geometric_factors_gll(refel, pts);
      
      % for now the mass matrix term is using gauss quadrature
      Md = refel.W .* J;   
      Me = refel.Q' * diag(Md) * refel.Q;

      % these are at the gll points 
      res1x = res1x + refel.Wgll .* Jg .* homg.tensor.IAX(refel.Dr', u1); 
      
      res1x = res1x + Me*q1x; 
      
      % y term
      res1y(idx_vol_1) = res1y(idx_vol_1) + lam_face_1;
      
      res1y = res1y + homg.tensor.AIX(refel.Dr', u1); 
      res1y = res1y + Me*q1y; 
      
      
    end
    
end
