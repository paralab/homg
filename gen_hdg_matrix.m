function A = gen_hdg_matrix(mesh, refel)
	% Dec 28, 2014
	% Compute the HDG Jacobian matrix, using the mesh

	% variables
	SkelAll2Interior = mesh.SkelAll2Interior; % to be added
	
	Nfp = refel.Nrp ^ (refel.dim - 1);
	Nfaces = refel.dim * 2;
	Nifaces = mesh.Nifaces;                   % to be added
	Nv = refel.Nrp ^ (refel.dim);
	
	K = prod(mesh.nelems);

	% @todo hard-coded for now, get from mesh later
	nx = [-1, 1, 0, 0];
	ny = [0, 0, -1, 1];

	taur = 1; % HDG.taur; pass as parameter to construction ?
	
	% Construct the Lift and VtoF
	LIFT = zeros(Nv, Nfp, Nfaces);
	VtoF = zeros(Nfp, Nv, Nfaces);
	for f = 1:Nfaces
	  idxv = mesh.get_discontinuous_face_indices(refel, 1, f);
	  LIFT(idxv,:,f) = refel.Mr;
	  for fp = 1:Nfp
	    VtoF(fp,idxv(fp),f) = 1;
	  end
	end

	% actual computation ...

	maxnnzeros = Nfaces * 2 * Nfp * Nfp * Nifaces;
	II = zeros(maxnnzeros,1);
	JJ = zeros(maxnnzeros,1);
	SS = zeros(maxnnzeros,1);

	Nfp2 = Nfp * Nfp;
	nnzeros = 0;
	I = eye(Nfp);

	for e = 1:K
		% e1 solution
		[du_dlam, dqx_dlam, dqy_dlam] = DifflocalSolver(mesh, refel, taur, LIFT, VtoF, e);
		for f = 1:Nfaces
			index = (f-1)*Nfp+1:f*Nfp;  
			idxf = mesh.get_skeletal_face_indices(refel, e, f);
			idxf_interior = SkelAll2Interior(idxf);

			% if this is boundary face, ignore
			if idxf_interior(1) < eps,
				continue;
			end

			Jf = mesh.geometric_factors_face(refel,e,f);
			idxv = mesh.get_discontinuous_face_indices(refel, 1, f);
    
    
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
				idxfn = mesh.get_skeletal_face_indices(refel, e, fn);
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


end

function [du_dlam, dqx_dlam, dqy_dlam] = DifflocalSolver(mesh, refel, taur, LIFT, VtoF, e1)
	% compute derivatives for the local solver
	
	% predefined normal vector, don't like it but stick with it for now
	nx = [-1, 1, 0, 0];
	ny = [0, 0, -1, 1];

	Nfp = refel.Nrp ^ (refel.dim - 1);
	Nv = refel.Nrp ^ (refel.dim);

	uu  = zeros(Nv,Nv);
	uqx = zeros(Nv,Nv);
	uqy = zeros(Nv,Nv);
      
	pts = mesh.element_nodes(e1, refel);
	[Jv, Dv] = mesh.geometric_factors(refel, pts);

	eMat = mesh.element_mass(e1, refel, Jv);
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

	% advection stiffness
	[Kex, Key] = mesh.element_stiffness_advection (e1, refel, Jv, Dv);

	uqx = Kex;
	uqy = Key;
	
	% residual for qx and qy equations
	for f = 1:Nfaces %
		index = (f-1)*Nfp+1:f*Nfp;
	
		% geometrix factors at gll points on face
		Jf = mesh.geometric_factors_face(refel,e1,f);
  
	  %% residual due to lambda
	  % rhsfx = Jf .* (refel.Mr * lamlocal) * nx(f);
	  %% lift to volume residual q equation
	  % rhsqx(idxv) = rhsqx(idxv) + rhsfx;
		drhsqx_dlam(:,index) = LIFT(:,:,f) * (diag(Jf .* nx(f)));

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
end
