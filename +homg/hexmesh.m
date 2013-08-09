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
        [x,y] = meshgrid(0:1/nelems(1):1.0, 0:1/nelems(2):1.0);
        pts = [x(:) y(:)];
      else
        [x,y,z] = meshgrid(0:1/nelems(1):1.0, 0:1/nelems(2):1.0,0:1/nelems(3):1.0);
        pts = [x(:) y(:) z(:)];
      end
      
      mesh.coords = X(pts);
    end
    
    function plot(mesh)
        % display the mesh. Needs X.
        figure(1);
        if (mesh.dim == 2 )
            plot(mesh.coords(:,1), mesh.coords(:,2), 'bo');
            axis square
        else
            plot3(mesh.coords(:,1), mesh.coords(:,2), mesh.coords(:,3), 'bo');
            view(3); axis square
        end
        % title(['Hex Mesh ', num2str(numx,3),'x',num2str(numy,3),'x',num2str(numz,3)])
    end
    
    function set_coeff(mesh, coeff)
      mesh.coeff = coeff;
    end

    function set_rhs(mesh, rhs)
      mesh.rhs = rhs;
    end
    
    function M = assemble_mass(mesh, order)
       % assemble the mass matrix
       refel = homg.refel ( mesh.dim, order );
       
       dof = prod(mesh.nelems*order + 1);
       ne  = prod(mesh.nelems);
       
       M = zeros(dof, dof); % maybe sparse ?
       
       % loop over elements
       for e=1:ne
         idx = mesh.get_node_indices (e, order);
         
         M(idx, idx) = mesh.element_mass(order, refel);
       end
    end
    
    
    function idx = get_node_indices ( mesh, elem, order )
       % determine global node indices for a given element
       if ( mesh.dim == 2) 
        [i,j] = ind2sub (mesh.nelems, elem);
        
        i_low   = (i-1)*order + 1;   i_high =  i*order + 1;
        j_low   = (j-1)*order + 1;   j_high =  j*order + 1;
        
        [i,j] = ndgrid(i_low:i_high, j_low:j_high);
       
        idx     = sub2ind (mesh.nelems*order + 1, i(:), j(:));
       else
        [i,j,k] = ind2sub (mesh.nelems, elem);
       
        i_low   = (i-1)*order + 1;   i_high =  i*order + 1;
        j_low   = (j-1)*order + 1;   j_high =  j*order + 1;
        k_low   = (k-1)*order + 1;   k_high =  k*order + 1;
       
        [i,j,k] = ndgrid(i_low:i_high, j_low:j_high, k_low:k_high);
       
        idx     = sub2ind (mesh.nelems*order + 1, i(:), j(:), k(:) );
       end
    end
    
    function Me = element_mass(mesh, elem, refel)
       % element mass matrix 
       Np =  refel.Nrp ^ mesh.dim;
       Me = zeros(Np);
    end
    
    function [J, rst] = geometric_factors( mesh, elem, refel ) 
        Np =  refel.Nrp ^ mesh.dim;
        
        J   = zeros(Np,1);
        rst = zeros(Np,mesh.dim^2);
        
        
        
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