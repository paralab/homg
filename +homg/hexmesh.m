classdef hexmesh < handle
  %HEXMESH A container class for homg regular grid meshes
  %      RGRID class for homg meshes, made of hexahedral elements.
  
  properties (SetAccess = private)
    dim=2;
    nelems=[8 8];
    
    % coords       % element vertices
    Xf           % transform
    
    % problem specific
    % coeff
    % rhs
  end % properties
  
  methods
    function mesh = hexmesh(nelems, X)
      % dim = 2/3,
      % nelems = number of elements per dimension,
      % X = transform / warping function
      
      mesh.dim    = length(nelems);
      mesh.nelems  = nelems;
      mesh.Xf = X;
    end
    
    function plot(self)
      % create coordinates ...
      if ( self.dim == 2 )
        [x,y] = ndgrid(0:1/self.nelems(1):1.0, 0:1/self.nelems(2):1.0);
        pts = [x(:) y(:)];
      else
        [x,y,z] = ndgrid(0:1/self.nelems(1):1.0, 0:1/self.nelems(2):1.0,0:1/self.nelems(3):1.0);
        pts = [x(:) y(:) z(:)];
      end
      
      coords = self.Xf(pts);
      % display the mesh. Needs X.
      figure(1);
      c = [31/256,171/256,226/256];   % default color of grid
      lw = 1;                         % default line width of grid
      
      if (self.dim == 2 )
        plot(coords(:,1), coords(:,2), 'ko');
        hold on;
        x = reshape(coords(:,1), self.nelems(1)+1, self.nelems(2)+1);
        y = reshape(coords(:,2), self.nelems(1)+1, self.nelems(2)+1);
        plot(x,y, 'Color',c,'LineWidth',lw);
        plot(x',y', 'Color',c,'LineWidth',lw);
        % boundary test
        idx = self.get_boundary_node_indices(1);
        plot(coords(idx,1), coords(idx,2), 'ro');
        axis square
      else
        plot3(coords(:,1), coords(:,2), coords(:,3), 'bo');
        hold on;
        idx = self.get_boundary_node_indices(1);
        plot3(coords(idx,1), coords(idx,2), coords(idx,3), 'ro');
        view(3); axis square
      end
      % title(['Hex Mesh ', num2str(numx,3),'x',num2str(numy,3),'x',num2str(numz,3)])
    end
    
    function u = evaluate(self, fx, order, where)
      % evaluate a function over the domain,
      % where is a string, 'gll', 'gauss' or 'uniform'
      
      if strcmp(where, 'gll')
        pfx = @homg.hexmesh.getGLLcoords;
      elseif strcmp(where, 'gauss')
        pfx = @homg.hexmesh.getGaussCoords;
      else
        pfx = @homg.hexmesh.getUniformCoords;
      end
      
      % create coordinates ...
      if ( self.dim == 2 )
        xp = pfx (order, self.nelems(1));
        yp = pfx (order, self.nelems(2));
        [x,y] = ndgrid(xp, yp);
        pts = [x(:) y(:)];
      else
        xp = pfx (order, self.nelems(1));
        yp = pfx (order, self.nelems(2));
        zp = pfx (order, self.nelems(3));
        [x,y,z] = ndgrid(xp, yp, zp);
        pts = [x(:) y(:) z(:)];
      end
      
      coords = self.Xf(pts);
      
      % fx = @(x,y) or @(x,y,z)
      if (self.dim == 2)
        u = arrayfun( fx, coords(:,1), coords(:,2) );
      else
        u = arrayfun( fx, coords(:,1), coords(:,2), coords(:,3) );
      end
      
      % u = reshape(u, self.nelems + 1);
    end
    
    function set_coeff(self, coeff)
      self.coeff = coeff;
    end
    
    function set_rhs(mesh, rhs)
      mesh.rhs = rhs;
    end
    
    function M = assemble_mass(self, order)
      % assemble the mass matrix
      refel = homg.refel ( self.dim, order );
      dof = prod(self.nelems*order + 1);
      ne  = prod(self.nelems);
      tic;
if (1)
      % storage for indices and values
      NP = (order+1)^self.dim;
      NPNP = NP * NP;
      eM = zeros(NP, NP);
      I = zeros(ne * NP, 1);
      J = zeros(ne * NP, 1);
      val = zeros(ne * NP, 1);
      
      % loop over elements
      for e=1:ne
        detJac = self.geometric_factors(e, refel);
        idx = self.get_node_indices (e, order);
        eM = self.element_mass(e, refel, detJac);
        ind1 = repmat(idx,NP,1);
        ind2 = reshape(repmat(idx',NP,1),NPNP,1);
        st = (e-1)*NPNP+1;
        en = e*NPNP;
        I(st:en) = ind1;
        J(st:en) = ind2;
        val(st:en) = eM(:);
      end
      M = sparse(I,J,val,dof,dof);
else

      num_nz =  dof * ( min(dof, (order+2)^self.dim) );
      M = spalloc(dof, dof, num_nz); 
      % loop over elements
      for e=1:ne
        idx = self.get_node_indices (e, order);
        
        detJac = self.geometric_factors(e, refel);
        M(idx, idx) = M(idx, idx) + self.element_mass(e, refel, detJac);
      end
end
      tspent = toc;
      fprintf('Mass: Assembly time: %g\n', tspent);
      % M = sparse(M);
    end
    
    function K = assemble_stiffness(self, order)
      % assemble the stiffness matrix
      refel = homg.refel ( self.dim, order );
      
      dof = prod(self.nelems*order + 1);
      ne  = prod(self.nelems);
      
      num_nz = dof * ( min(dof, (order+2)^self.dim) ); 
      
      K = spalloc(dof, dof, num_nz); 
      
      % loop over elements
      for e=1:ne
        idx = self.get_node_indices (e, order);
        [detJac, Jac] = self.geometric_factors(e, refel);
        
        K(idx, idx) = K(idx, idx) + self.element_stiffness(e, refel, detJac, Jac);
      end
    end
    
    function [K, M] = assemble_poisson(self, order)
      % assemble the mass matrix
      refel = homg.refel ( self.dim, order );
      dof = prod(self.nelems*order + 1);
      ne  = prod(self.nelems);
      tic;
      
      % storage for indices and values
      NP = (order+1)^self.dim;
      NPNP = NP * NP;
      eMat = zeros(NP, NP);
      
      I = zeros(ne * NP, 1);
      J = zeros(ne * NP, 1);
      mass_val = zeros(ne * NP, 1);
      stiff_val = zeros(ne * NP, 1);
      
      % loop over elements
      for e=1:ne
        idx = self.get_node_indices (e, order);
        
        ind1 = repmat(idx,NP,1);
        ind2 = reshape(repmat(idx',NP,1),NPNP,1);
        st = (e-1)*NPNP+1;
        en = e*NPNP;
        I(st:en) = ind1;
        J(st:en) = ind2;
        
        [detJac, Jac] = self.geometric_factors(e, refel);
        
        eMat = self.element_mass(e, refel, detJac);
        mass_val(st:en) = eMat(:);
        
        eMat = self.element_stiffness(e, refel, detJac, Jac);
        stiff_val(st:en) = eMat(:);
      end
      M = sparse(I,J,mass_val,dof,dof);
      % zero dirichlet bdy conditions
      bdy = self.get_boundary_node_indices(order);

      ii = ismember(I,bdy);
      jj = ismember(J,bdy);
      
      stiff_val = stiff_val.*(~ii).*(~jj);
      I = [I; bdy];
      J = [J; bdy];
      stiff_val = [stiff_val; ones(length(bdy), 1)];
      
      K = sparse(I,J,stiff_val,dof,dof);
    end

% 
%     function [K, M] = assemble_poisson(self, order)
%       % assemble the mass matrix
%       refel = homg.refel ( self.dim, order );
%       dof = prod(self.nelems*order + 1);
%       ne  = prod(self.nelems);
%       
%       num_nz = dof * ( min(dof, (order+2)^self.dim) ); 
%       
%       M = spalloc(dof, dof, num_nz); 
%       K = spalloc(dof, dof, num_nz); 
%       
%       % loop over elements
%       for e=1:ne
%         idx = self.get_node_indices (e, order);
%         [detJac, Jac] = self.geometric_factors(e, refel);
%         
%         M(idx, idx) = M(idx, idx) + self.element_mass(e, refel, detJac);
%         K(idx, idx) = K(idx, idx) + self.element_stiffness(e, refel, detJac, Jac);
%       end
%     end

    function f = assemble_rhs(self, fx, order)
      % assemble the mass matrix
      refel = homg.refel ( self.dim, order );
      
      dof = prod(self.nelems*order + 1);
      ne  = prod(self.nelems);
      
      f = zeros(dof,1);
      fval = self.evaluate(fx, order, 'gauss');
      
      % loop over elements
      for e=1:ne
        idx = self.get_node_indices (e, order);
        J = self.geometric_factors(e, refel);
        Jd = refel.W .* J; 
        fd =  fval(idx); 
      
        f(idx) = f(idx) + refel.Q' * diag(Jd) * fd;
      end
      % M = sparse(M);
    end
    
    function P = assemble_interpolation(self, order)
      % assemble prolongation operator from coarse (self) to fine mesh
      refel = homg.refel ( self.dim, order );
      
      dof_coarse = prod(self.nelems*order + 1);
      dof_fine   = prod(2*self.nelems*order + 1);
      
      ne  = prod(self.nelems);
      
      num_nz = dof_coarse *  (3*order)^self.dim ;
      
      % P = sparse(dof_fine, dof_coarse);
      P = spalloc(dof_fine, dof_coarse, num_nz); 
      
      % loop over elements
      for e=1:ne
        [idx_coarse, idx_fine] = self.get_interpolation_indices (e, order);
        
        P(idx_fine, idx_coarse) = refel.P;
      end
      
      % disp(['factor = ' num2str(nnz(P)/dof_coarse)]); 
    end
    
    function idx = get_node_indices ( self, eid, order )
      % determine global node indices for a given element
      if ( self.dim == 2)
        [i,j] = ind2sub (self.nelems, eid);
        
        i_low   = (i-1)*order + 1;   i_high =  i*order + 1;
        j_low   = (j-1)*order + 1;   j_high =  j*order + 1;
        
        [i,j] = ndgrid(i_low:i_high, j_low:j_high);
        
        idx     = sub2ind (self.nelems*order + 1, i(:), j(:));
      else
        [i,j,k] = ind2sub (self.nelems, eid);
        
        i_low   = (i-1)*order + 1;   i_high =  i*order + 1;
        j_low   = (j-1)*order + 1;   j_high =  j*order + 1;
        k_low   = (k-1)*order + 1;   k_high =  k*order + 1;
        
        [i,j,k] = ndgrid(i_low:i_high, j_low:j_high, k_low:k_high);
        
        idx     = sub2ind (self.nelems*order + 1, i(:), j(:), k(:) );
      end
    end
    
    function [idx_coarse, idx_fine] = get_interpolation_indices ( self, eid, order )
      % determine global node indices for a given element
      if ( self.dim == 2)
        [i,j] = ind2sub (self.nelems, eid);
        
        i_low       = (i-1)*order + 1;   i_high =  i*order + 1;
        j_low       = (j-1)*order + 1;   j_high =  j*order + 1;
        [i,j]       = ndgrid(i_low:i_high, j_low:j_high);
        idx_coarse  = sub2ind (self.nelems*order + 1, i(:), j(:));
        
        [i,j]       = ndgrid(2*i_low-1:2*i_high-1, 2*j_low-1:2*j_high-1);
        idx_fine    = sub2ind (2*self.nelems*order + 1, i(:), j(:));
      else
        [i,j,k] = ind2sub (self.nelems, eid);
        
        i_low       = (i-1)*order + 1;   i_high =  i*order + 1;
        j_low       = (j-1)*order + 1;   j_high =  j*order + 1;
        k_low       = (k-1)*order + 1;   k_high =  k*order + 1;
        [i,j,k]     = ndgrid(i_low:i_high, j_low:j_high, k_low:k_high);
        idx_coarse  = sub2ind (self.nelems*order + 1, i(:), j(:), k(:) );
      
        [i,j,k]     = ndgrid(2*i_low-1:2*i_high-1, 2*j_low-1:2*j_high-1, 2*k_low-1:2*k_high-1);
        idx_fine    = sub2ind (2*self.nelems*order + 1, i(:), j(:), k(:) );
      end
    end
    
		function idx = get_boundary_node_indices(self, order)
			% function idx = get_boundary_node_indices(self, order)
      %    returns indices of boundary nodes, for setting 
      %    boundary conditions       
      if (self.dim == 2)
        [x,y] = ndgrid(1:self.nelems(1)*order+1,1:self.nelems(2)*order+1);
        
        idx = [ find(x == 1);     find(x == (self.nelems(1)*order+1));
                find(y == 1);     find(y == (self.nelems(2)*order+1))  ];
        
        idx = unique(idx);
      else 
         [x,y,z] = ndgrid(1:self.nelems(1)*order+1,1:self.nelems(2)*order+1,1:self.nelems(3)*order+1);
         
         idx = [ find(x == 1);     find(x == (self.nelems(1)*order+1));
                 find(y == 1);     find(y == (self.nelems(2)*order+1));
                 find(z == 1);     find(z == (self.nelems(3)*order+1))  ];
        
        idx = unique(idx);
      end
		end
		
    function Me = element_mass(self, eid, refel, J)
      % element mass matrix
      % J = self.geometric_factors(eid, refel);
      Md = refel.W .* J ; 
      
      Me = refel.Q' * diag(Md) * refel.Q;
    end
    
    function Ke = element_stiffness(self, eid, r, J, D)
      % element mass matrix
      % [J, D] = self.geometric_factors(eid, r);
      
%             | Qx Qy Qz || rx ry rz |     | rx sx tx || Qx |
%    Ke =                 | sx sy sz | J W | ry sy ty || Qy |
%                         | tx ty tz |     | rz sz tz || Qz |
      
      
      factor = zeros(length(J), 6);

      %             1  4  5
      % factor      4  2  6
      %             5  6  3
      
      
      if (self.dim == 2 )
        factor (:,1) = (D.rx.*D.rx + D.ry.*D.ry ) .* J .* r.W ; % d2u/dx^2
        factor (:,2) = (D.sx.*D.sx + D.sy.*D.sy ) .* J .* r.W ; % d2u/dy^2
        factor (:,3) = (D.rx.*D.sx + D.ry.*D.sy ) .* J .* r.W ; % d2u/dxdy
        
        Ke =   r.Qx' * diag(factor(:,1)) * r.Qx ...
             + r.Qy' * diag(factor(:,2)) * r.Qy ...
             + r.Qx' * diag(factor(:,3)) * r.Qy ...
             + r.Qy' * diag(factor(:,3)) * r.Qx ;
      else
        
        % first compute dj.w.J.J'
        factor (:,1) = (D.rx.*D.rx + D.ry.*D.ry + D.rz.*D.rz ) .* J .* r.W ; % d2u/dx^2
        factor (:,2) = (D.sx.*D.sx + D.sy.*D.sy + D.sz.*D.sz ) .* J .* r.W ; % d2u/dy^2
        factor (:,3) = (D.tx.*D.tx + D.ty.*D.ty + D.tz.*D.tz ) .* J .* r.W ; % d2u/dz^2
        factor (:,4) = (D.rx.*D.sx + D.ry.*D.sy + D.rz.*D.sz ) .* J .* r.W ; % d2u/dxdy
        factor (:,5) = (D.rx.*D.tx + D.ry.*D.ty + D.rz.*D.tz ) .* J .* r.W ; % d2u/dxdz
        factor (:,6) = (D.sx.*D.tx + D.sy.*D.ty + D.sz.*D.tz ) .* J .* r.W ; % d2u/dydz
        
        Ke =   r.Qx' * diag(factor(:,1)) * r.Qx ...
             + r.Qy' * diag(factor(:,2)) * r.Qy ...
             + r.Qz' * diag(factor(:,3)) * r.Qz ...
             + r.Qx' * diag(factor(:,4)) * r.Qy ...
             + r.Qy' * diag(factor(:,4)) * r.Qx ...
             + r.Qx' * diag(factor(:,5)) * r.Qz ...
             + r.Qz' * diag(factor(:,5)) * r.Qx ...
             + r.Qz' * diag(factor(:,6)) * r.Qy ...
             + r.Qy' * diag(factor(:,6)) * r.Qz ;
      end
      
    end
    
    
    function [J, D] = geometric_factors( self, elem, refel )
      % Np =  refel.Nrp ^ mesh.dim;
      
      % compute x,y,z for element
      pts = self.element_nodes(elem, refel);
        
      % change to using Qx etc ?
      if (refel.dim == 2)
        [xr, xs] = homg.tensor.grad2 (refel.Dr, pts(:,1));
        [yr, ys] = homg.tensor.grad2 (refel.Dr, pts(:,2));
      
        J = -xs.*yr + xr.*ys;
      else
        [xr, xs, xt] = homg.tensor.grad3 (refel.Dr, pts(:,1));
        [yr, ys, yt] = homg.tensor.grad3 (refel.Dr, pts(:,2));
        [zr, zs, zt] = homg.tensor.grad3 (refel.Dr, pts(:,3));
        
        J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
      end
      
      if (nargout > 1)
        if (refel.dim == 2)
          D.rx =  ys./J;
          D.sx = -yr./J;
          D.ry = -xs./J;
          D.sy =  xr./J;
        else
          D.rx =  (ys.*zt - zs.*yt)./J;
          D.ry = -(xs.*zt - zs.*xt)./J;
          D.rz =  (xs.*yt - ys.*xt)./J;
          
          D.sx = -(yr.*zt - zr.*yt)./J;
          D.sy =  (xr.*zt - zr.*xt)./J;
          D.sz = -(xr.*yt - yr.*xt)./J;
          
          D.tx =  (yr.*zs - zr.*ys)./J;
          D.ty = -(xr.*zs - zr.*xs)./J;
          D.tz =  (xr.*ys - yr.*xs)./J;
        end
        
      end
    end
    
    function coords = element_nodes(self, elem, refel)
      h = 1./self.nelems;
      
      if ( self.dim == 2)
        [i,j] = ind2sub (self.nelems, elem);
        idx = [i j];
      else
        [i,j,k] = ind2sub (self.nelems, elem);
        idx = [i j k];
      end
      
      p_mid = (idx + 0.5) .* h;
      p_gll = refel.r * 0.5 * h;
      nodes = bsxfun(@plus, p_mid, p_gll) ;
      
      if ( self.dim == 2)
        [x, y] = ndgrid(nodes(:,1), nodes(:,2));
        pts = [x(:) y(:)];
      else
        [x, y, z] = ndgrid(nodes(:,1), nodes(:,2), nodes(:,3));
        pts = [x(:) y(:) z(:)];
      end
      
      coords = self.Xf(pts);
    end
		
  end % methods
  
  methods(Static) 
    function coords = getGLLcoords(order, elems)
      % function coords=getGLLcoords(order, elems)
      % returns location of gll coordinates of order
      % for elements in [0,1]
      
      fac = 1.0/(2*elems);
      
      % gll coordinates in [-1,1]
      x = homg.basis.gll (0,0,order)';
      
      x = (x + 1)*fac;
      
      
      coords = [];
      for i=1:elems
        y = x + (i-1)/elems;
        coords = [coords y(1:end-1)];
      end
      
      coords = [coords 1.0];
    end
    
    function coords = getGaussCoords(order, elems)
      % function coords=getGaussCoords(order, elems)
      % returns location of gauss coordinates of order
      % for elements in [0,1]
      
      fac = 1.0/(2*elems);
      
      % gll coordinates in [-1,1]
      x = homg.basis.gauss (0,0,order)';
      
      x = (x + 1)*fac;
      
      
      coords = [];
      for i=1:elems
        y = x + (i-1)/elems;
        coords = [coords y(1:end-1)];
      end
      
      coords = [coords 1.0];
    end
      
    function coords = getUniformCoords(order, elems)
      coords = linspace(0, 1, order*elems+1);
    end
      
    function C = stats(nelems, order)
      % function Ch = stats(nelems, order)
      %   given number of elements and the order,
      %   this function calculates different node 
      %   stats for the mesh
      d               = length(nelems);
      C.num_nodes     = prod(nelems*order + 1);
      C.num_elements  = prod(nelems);
      C.num_bdy_nodes = C.num_nodes - prod(nelems*order - 1);
      
      C.num_int_elements = prod(nelems - 1);
      C.num_bdy_elements = C.num_elements - C.num_int_elements;
      
      C.nnz = (order+2)^d*C.num_nodes;
%       if (d == 2)
%         
%       else
%         
%       end
    end
  end  % static methods
  
  
end % class
