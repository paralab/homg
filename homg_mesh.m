classdef homg_mesh < handle
  %HOMG_MESH A container class for homg Meshes
  %   Mesh class for homg meshes. Is a wrapper around comsol
  
  properties (SetAccess = private)
    dim=2;
    nelem=8;
    fem
    rhs
    % permutations
    perm_full
    perm_rest
    coords
  end % properties
  
  methods
    function mesh = homg_mesh(dim,nelem)
      % dim = 2/3, nelem = number of elements per dimension.
      if nargin > 0
        if nargin > 1
          mesh.dim    = dim;
        end
        mesh.nelem  = nelem;
      end
      
      if (mesh.dim == 2)
        mesh.fem.geom = rect2(0,1,0,1);
        mesh.fem.mesh = meshmap(mesh.fem, 'Edgelem', {1, mesh.nelem,2, mesh.nelem}, 'report', 'off');
      elseif (mesh.dim == 3)
        tmp_fem.geom = rect2(0,1,0,1);
        tmp_fem.mesh = meshmap(tmp_fem, 'Edgelem', {1, mesh.nelem,2, mesh.nelem}, 'report', 'off');
        mesh.fem = meshextrude(tmp_fem, 'distance', 1, 'elextlayers', {mesh.nelem});
        clear tmp_fem;
      else
        error(['Error: mesh is not supported for ',num2str(mesh.dim),'D.'])
      end
      % default rhs
      mesh.rhs = '1';
    end % constructor
    
    function show(mesh)
      % display the mesh. Needs X.
      meshplot(mesh.fem);
    end
    
    function set_rhs(mesh, rhs)
      mesh.rhs = rhs;
    end
    
    function M = assemble_mass(mesh, order)
      % Assembles the mass matrix
      % Check that the Stiffness matrix is assembled first
      mesh.fem.dim = {'u'};
      mesh.fem.shape = order;
      
      mesh.fem.equ.weak = '(-u*u_test)';
      mesh.fem.xmesh = meshextend(mesh.fem,  'report', 'off');
      
      M = assemble(mesh.fem, 'Out', {'K'}, 'report', 'off');
      
      if ( isempty(mesh.perm_full) )
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
      end
      M = M(mesh.perm_full, mesh.perm_full);
    end
    
    function [K,L,Null,Ud] = assemble_poisson(mesh, order)
      % Assembles the Stiffness matrix and RHS for the poisson equation
      mesh.fem.dim   = {'u'};
      mesh.fem.shape = order;
      
      mesh.fem.equ.expr = {'f' mesh.rhs};
      if (mesh.dim == 2)
        mesh.fem.equ.weak = '-(ux * ux_test + uy * uy_test + f * u_test)';
      else
        mesh.fem.equ.weak = '-(ux * ux_test + uy * uy_test + uz * uz_test + f * u_test)';
      end
      
      mesh.fem.bnd.r = {'u-0'};
      mesh.fem.xmesh = meshextend(mesh.fem,  'report', 'off');
      
      % K system matrix, L rhs, boundary conditions are to be incorporated
      % by imposing N*U = M
      [K,L] = assemble(mesh.fem,'Out',{'K','L'}, 'report', 'off');
      
      [~,~,Null,Ud] = femlin(mesh.fem,  'report', 'off');
      
      nodes = xmeshinfo(mesh.fem ,'out', 'nodes');
      dofs = nodes.dofs;
      crds = nodes.coords';
      
      % mesh.fem.sol = femlin(mesh.fem, 'report', 'off');
      
      [~,idof] = sort(dofs);
      
      mesh.coords = crds(idof,:);
      fac = 10*mesh.nelem;
      
      if (mesh.dim == 2)
        sortval = mesh.coords(:,2)*fac + mesh.coords(:,1);
      else
        sortval = mesh.coords(:,3)*fac*fac + mesh.coords(:,2)*fac + mesh.coords(:,1);
      end
      
      [~, p] = sort(sortval);
      mesh.perm_full = p;
      
      crds2 = Null' * mesh.coords;
      if (mesh.dim == 2)
        sortval = crds2(:,2)*fac + crds2(:,1);
      else
        sortval = crds2(:,3)*fac*fac + crds2(:,2)*fac + crds2(:,1);
      end
      
      [~, pc] = sort(sortval);
      mesh.perm_rest = pc;
      
      % re-arrange matrices and coords ...
      K     = K(p, p);
      L     = L(p);
      Null  = Null(p, pc);
      Ud    = Ud(p);
      
      mesh.coords = mesh.coords(p,:);
      
    end
    
    function P = assemble_interpolation(mesh, order, pts)
      % todo: check if matrices assembled before this ...
      % vector of FEM coefficients and its length
      if (order ~= mesh.fem.shape)
        disp('The order specified in the interpolation does not match the one used for assembly');
      end
      
      % X = mesh.fem.sol.u;
      % no_dofs = length(X);
      no_dofs =  (order*mesh.nelem+1)^mesh.dim;
      X = zeros(no_dofs,1);
      
      % allocate storage for interpolation operator
      P = zeros(size(pts,2),no_dofs);
      
      % build interpolation operator
      for i = 1:no_dofs
        X(:) = 0;
        X(i) = 1;
        P(:,i) = postinterp(mesh.fem, 'u', pts, 'U', X)';
      end
      
    end
    
  end % methods
  
end % classdef

