classdef homg_mesh
  %HOMG_MESH A container class for homg Meshes
  %   Mesh class for homg meshes. Is a wrapper around comsol 
  
  properties (SetAccess = private)
    dim=2;
    nelem=8; 
    fem
    % permutations
    perm_full 
    perm_rest
    coords
  end % properties 
  
  methods
    function mesh = homg_mesh(dim,nelem)
      mesh.dim    = dim;
      mesh.nelem  = nelem;
      mesh.perm_full = [];
      
      if (mesh.dim == 2)
        mesh.fem.geom = rect2(0,1,0,1);
        mesh.fem.mesh = meshmap(mesh.fem, 'Edgelem', {1,nelem,2,nelem});
      elseif (mesh.dim == 3)
        tmp_fem.geom = rect2(0,1,0,1);
        mesh.fem = meshextrude(tmp_fem, 'distance', 1, 'elextlayers', {nelem});
        clear tmp_fem;
      else
        disp(['Error: mesh is not supported for ',num2str(mesh.dim),' dimensions.'])
      end
    end % constructor
    
    function M = assemble_mass(mesh, order)
      mesh.fem.dim = {'u'};
      mesh.fem.shape = order;
      
      mesh.fem.equ.weak = '(u*u_test)';
      M = assemble(mesh.fem);
    end
    
    function [K,L,M,N] = assemble_poisson(mesh, order)
      mesh.fem.dim   = {'u'};
      mesh.fem.shape = order;
      
      if (mesh.dim == 2)
        mesh.fem.equ.weak = '-(ux * ux_test + uy * uy_test + 1*u_test)';
      else
        mesh.fem.equ.weak = '-(ux * ux_test + uy * uy_test + uz * uz_test + 1*u_test)';
      end
      
      mesh.fem.bnd.r = {'u-0'};
      mesh.fem.xmesh = meshextend(mesh.fem);
      
      % K system matrix, L rhs, boundary conditions are to be incorporated
      % by imposing N*U = M
      [K,L,M,N] = assemble(mesh.fem);
      
      % [Kc,Lc,Null,Ud] = femlin(fem);
      
      nodes = xmeshinfo(mesh.fem ,'out', 'nodes');
      dofs = nodes.dofs;
      crds = nodes.coords';

      [~,idof] = sort(dofs);

      mesh.coords = crds(idof,:);
      fac = 2*mesh.nelem;
      
      if (mesh.dim == 2)
        sortval = mesh.coords(:,2)*fac + mesh.coords(:,1);
      else
        sortval = mesh.coords(:,3)*fac*fac + mesh.coords(:,2)*fac + mesh.coords(:,1);
      end
      
      [~, p] = sort(sortval);
      mesh.perm_full = p; 
      %crds2 = Null' * crds;
      %sortval = crds2(:,3)*1e6 + crds2(:,2)*1e3 + crds2(:,1);
      %[~,mesh.p2] = sort(sortval);
    end
    
    function P = assemble_interpolation(mesh, order, pts)
      % todo: check if matrices assembled before this ...
      % vector of FEM coefficients and its length
      if (order ~= mesh.fem.shape)
        disp('The order specified in the interpolation does not match the one used for assembly');
      end
      
      X = mesh.fem.sol.u; 
      no_dofs = length(X);
      
      % allocate storage for interpolation operator
      P = zeros(numel(pts),no_dofs);

      % build interpolation operator
      for i = 1:no_dofs
          X(:) = 0;
          X(i) = 1;
          P(:,i) = postinterp(mesh.fem, 'u', pts, 'U', X)';
      end

    end
    
  end % methods
  
end % classdef

