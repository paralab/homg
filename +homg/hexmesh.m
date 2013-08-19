classdef hexmesh < handle
  %HEXMESH A container class for homg regular grid meshes
  %      RGRID class for homg meshes, made of hexahedral elements.
  
  properties (SetAccess = private)
    dim=2;
    nelems=[8 8];
    
    coords       % element vertices
    Xf           % transform
    
    % problem specific
    coeff
    rhs
  end % properties
  
  methods
    function mesh = hexmesh(nelems, X)
      % dim = 2/3,
      % nelems = number of elements per dimension,
      % X = transform / warping function
      
      mesh.dim    = length(nelems);
      mesh.nelems  = nelems;
      mesh.Xf = X;
      
      % create coordinates ...
      if ( mesh.dim == 2 )
        [x,y] = ndgrid(0:1/nelems(1):1.0, 0:1/nelems(2):1.0);
        pts = [x(:) y(:)];
      else
        [x,y,z] = ndgrid(0:1/nelems(1):1.0, 0:1/nelems(2):1.0,0:1/nelems(3):1.0);
        pts = [x(:) y(:) z(:)];
      end
      
      mesh.coords = X(pts);
    end
    
    function plot(self)
      % display the mesh. Needs X.
      figure(1);
      c = [31/256,171/256,226/256];   % default color of grid
      lw = 1;                         % default line width of grid
      
      if (self.dim == 2 )
        plot(self.coords(:,1), self.coords(:,2), 'ko');
        hold on;
        x = reshape(self.coords(:,1), self.nelems(1)+1, self.nelems(2)+1);
        y = reshape(self.coords(:,2), self.nelems(1)+1, self.nelems(2)+1);
        plot(x,y, 'Color',c,'LineWidth',lw);
        plot(x',y', 'Color',c,'LineWidth',lw);
        axis square
      else
        plot3(self.coords(:,1), self.coords(:,2), self.coords(:,3), 'bo');
        
        view(3); axis square
      end
      % title(['Hex Mesh ', num2str(numx,3),'x',num2str(numy,3),'x',num2str(numz,3)])
    end
    
    function u = evaluate(self, fx)
      % evaluate a function over the domain,
      % fx = @(x,y) or @(x,y,z)
      if (self.dim == 2)
        u = arrayfun( fx, self.coords(:,1), self.coords(:,2) );
      else
        u = arrayfun( fx, self.coords(:,1), self.coords(:,2), self.coords(:,3) );
      end
      
      u = reshape(u, self.nelems + 1);
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
      
      M = zeros(dof, dof); % maybe sparse ?
      
      % loop over elements
      for e=1:ne
        idx = self.get_node_indices (e, order);
        
        M(idx, idx) = M(idx, idx) + self.element_mass(e, refel);
      end
      M = sparse(M);
    end
    
    function K = assemble_stiffness(self, order)
      % assemble the mass matrix
      refel = homg.refel ( self.dim, order );
      
      dof = prod(self.nelems*order + 1);
      ne  = prod(self.nelems);
      
      K = zeros(dof, dof); % maybe sparse ?
      
      % loop over elements
      for e=1:ne
        idx = self.get_node_indices (e, order);
        
        K(idx, idx) = K(idx, idx) + self.element_stiffness(e, refel);
      end
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
    
    function Me = element_mass(self, eid, refel)
      % element mass matrix
      J = self.geometric_factors(eid, refel);
      Md = refel.W .* J ; 
      
      Me = refel.Q' * diag(Md) * refel.Q;
    end
    
    function Ke = element_stiffness(self, eid, r)
      % element mass matrix
      [J, D] = self.geometric_factors(eid, r);
      
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
    
    function coords = element_nodes(self,  elem, refel)
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
    function Xout = identity (Xin)
      Xout = Xin;
    end
    
    function Xout = shell (Xin)
      d = size(Xin, 2);
      R1 = 1.0; % hard coded for now.
      R2 = 0.55; % hard coded for now.
      R2byR1 = R2 / R1;
      R1sqrbyR2 = R1 * R1 / R2;
      
      if (d == 2)
        x = zeros( size(Xin(:,1)) );
        y = tan ( Xin(:,2)  * pi/4 );
        R = R1sqrbyR2 * ( R2byR1.^(Xin(:,1) + 1) ) ;
      else
        x = tan ( Xin(:,1)  * pi/4 );
        y = tan ( Xin(:,2)  * pi/4 );
        R = R1sqrbyR2 * ( R2byR1.^(Xin(:,3) + 1) );
      end
      
      q = R ./ sqrt (x.*x + y.*y + 1);
      
      if (d == 3)
        Xout(:,1) =  q.* y;
        Xout(:,2) = -q.* x;
        Xout(:,3) = q;
      else
        Xout(:,1) =   q.* y;
        Xout(:,2) =   q;
      end
    end
  end  % static methods
  
  
end % class