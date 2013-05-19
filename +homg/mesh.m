classdef mesh < handle
  %MESH A container class for homg Meshes
  %   Mesh class for homg meshes. Is a wrapper around comsol
  
  properties (SetAccess = private)
    dim=2;
    nelem=8;
    fem
    rhs
    coeff
    % permutations
    perm_full
    perm_rest
    coords
    geom_shape
  end % properties
  
  methods
    function mesh = mesh(dim, geom, nelem, sp)
      % dim = 2/3, nelem = number of elements per dimension.
      if nargin > 0
        if nargin > 1
          mesh.dim    = dim;
        end
        mesh.nelem  = nelem;
      end
      
      mesh.geom_shape = geom;

      if (mesh.dim == 2)
        if (geom == 'box')
          mesh.fem.geom = rect2(0,sp,0,1);
          mesh.fem.mesh = meshmap(mesh.fem, 'Edgelem', {1, mesh.nelem,2, mesh.nelem}, 'report', 'off');
        elseif (geom == 'fan')
          % Geometry
          g1 = ellip2(1.0, 1.0, 'base','center','pos',[0.,0.0]);
          g2=ellip2(0.8,0.8,'base','center','pos',[0.0,0.0]);
          % g2=ellip2(0.55,0.55,'base','center','pos',[0.0,0.0]);
          g3=geomcomp({g1,g2},'ns',{'g1','g2'},'sf','g1-g2','edge','none');
          g4=rect2(1.0,1.0,'base','corner','pos',[0.0,0.0]);
          g5=geomcomp({g3,g4},'ns',{'g3','g4'},'sf','g3*g4','edge','none');

          % Analyzed geometry
          s.objs={g5};
          s.name={'Sector'};
          % s.tags={'g5'};

          mesh.fem.draw=struct('s',s);
          mesh.fem.geom = geomcsg(mesh.fem);

          % Create mapped quad mesh
          mesh.fem.mesh=meshmap(mesh.fem, ...
                           'edgegroups',{{[4],[2],[3],[1]}}, ...
                           'Edgelem', {1,mesh.nelem,2,mesh.nelem,3,3*mesh.nelem,4,3*mesh.nelem}, 'report', 'off'); 
          
          clear g* s
        else
          error(['Error: Unknown mesh shape' geom])
        end

      elseif (mesh.dim == 3)
        if (geom == 'box')
          tmp_fem.geom = rect2(0,1,0,1);
          tmp_fem.mesh = meshmap(tmp_fem, 'Edgelem', {1, mesh.nelem,2, mesh.nelem}, 'report', 'off');
          mesh.fem = meshextrude(tmp_fem, 'distance', 1, 'elextlayers', {mesh.nelem});
          clear tmp_fem;
        else
           % Geometry
          g1 = ellip2(1.0, 1.0, 'base','center','pos',[0.,0.0]);
          % g2=ellip2(0.8,0.8,'base','center','pos',[0.0,0.0]);
          g2=ellip2(0.55,0.55,'base','center','pos',[0.0,0.0]);
          g3=geomcomp({g1,g2},'ns',{'g1','g2'},'sf','g1-g2','edge','none');
          g4=rect2(1.0,1.0,'base','corner','pos',[0.0,0.0]);
          g5=geomcomp({g3,g4},'ns',{'g3','g4'},'sf','g3*g4','edge','none');

          % Analyzed geometry
          s.objs={g5};
          s.name={'Sector'};
          % s.tags={'g5'};

          tmp_fem.draw=struct('s',s);
          tmp_fem.geom = geomcsg(tmp_fem);

          % Create mapped quad mesh
          tmp_fem.mesh=meshmap(tmp_fem, ...
                           'edgegroups',{{[4],[2],[3],[1]}}, ...
                           'Edgelem', {1,mesh.nelem,2,mesh.nelem,3,3*mesh.nelem,4,3*mesh.nelem}, 'report', 'off'); 
          
          mesh.fem = meshextrude(tmp_fem, 'distance', 1, 'elextlayers', {mesh.nelem});
          clear g* s tmp_fem;
          % error(['Error: Unknown mesh shape' geom])
        end
      else
        error(['Error: mesh is not supported for ',num2str(mesh.dim),'D.'])
      end
      % default rhs
      mesh.rhs   = '1';
      mesh.coeff = '1';
    end % constructor
    
    function plot(mesh)
      % display the mesh. Needs X.
      meshplot(mesh.fem);
    end
    
    function show(mesh)
      % display the mesh. Needs X.
      meshplot(mesh.fem);
    end
    
    function set_coeff(mesh, coeff)
      mesh.coeff = coeff;
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
      
      mesh.fem.equ.expr = {'f' mesh.rhs 'q' mesh.coeff};
      % mesh.fem.equ.expr = {'mu' mesh.coeff};
      if (mesh.dim == 2)
        mesh.fem.equ.weak = '-(u*u_test + q*ux*ux_test + q*uy*uy_test + f * u_test)';
      else
        mesh.fem.equ.weak = '-(u*u_test + q*ux*ux_test + q*uy*uy_test + q*uz*uz_test + f*u_test)';
      end
      
      % mesh.fem.bnd.r = {'u-0'};
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
    
    function P = assemble_interpolation(mesh, pts)
      % todo: check if matrices assembled before this ...
      % vector of FEM coefficients and its length

      % making change to allow order to be different ...
      % if (order ~= mesh.fem.shape)
        % if ( order ~= 1 )
        %  disp('The order specified in the interpolation does not match the one used for assembly');
        % else
          % Create a new fem that is linear ...
         % num_dups = mesh.fem.shape / order;
         % lin_nelem = num_dups * mesh.nelem;
         % meshlin = homg.mesh(mesh.dim, lin_nelem);
         %  meshlin.assemble_poisson(order);
         %
         % no_dofs =  (order*lin_nelem+1)^meshlin.dim;
         %  X = zeros(no_dofs,1);
         %
         % % allocate storage for interpolation operator
         %  P = zeros(size(pts,1), no_dofs);
         %
         % % build interpolation operator
         %  for i = 1:no_dofs
         %   X(:) = 0;
         %   X(i) = 1;
         %    P(:,i) = postinterp(meshlin.fem, 'u', pts', 'U', X)';
         % end
         % P = P(:, meshlin.perm_full);
         %  P = sparse(P);

        % end % if order != 1
      % else 
        
       %if (mesh.geom_shape == 'box')
       %   no_dofs =  (mesh.fem.shape * mesh.nelem + 1)^mesh.dim;
       % else
       no_dofs = size(mesh.coords, 1);
       % end
       elsize = 1.5 / mesh.nelem;

        Xi = zeros(no_dofs,1);

        % allocate storage for interpolation operator
        P = zeros(size(pts,1), no_dofs);
        prog_step = ceil(no_dofs/100);
        % build interpolation operator
	for i = 1:no_dofs
	    crds = mesh.coords(i,:);
	    if (mesh.dim == 2)
        ind = find(((pts(:,1)-crds(1)).^2 + (pts(:,2)-crds(2)).^2)<=elsize^2);
      else 
        ind = find(((pts(:,1)-crds(1)).^2 + (pts(:,2)-crds(2)).^2 + (pts(:,3)-crds(3)).^2)<=elsize^2);
      end

      % if ( mod(i, prog_step) == 0 )
	    %   disp(['Assembling P ' num2str(i/prog_step) ' %']);
	    % end
	    Xi(:) = 0;
	    %% Xi = zeros(no_dofs,1);
	    Xi(mesh.perm_full(i)) = 1;
	    P(ind,mesh.perm_full(i)) = postinterp(mesh.fem, 'u', pts(ind,:)', 'U', Xi)';
	   % plot(pts(ind,1),pts(ind,2)','ro'); hold on;
	   % plot(crds(1),crds(2),'bo','MarkerFaceColor', 'b'); hold off;
	   % axis([0,1,0,1]);
	   % pause(.4);
       end
	P = P(:, mesh.perm_full);
        P = sparse(P);
      % end 
    end 
end % methods

end % classdef

