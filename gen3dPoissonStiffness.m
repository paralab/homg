function genPoissonStiffness(p, nelem, typ)
% function genPoissonStiffness(N, nelem, type)
%
% generate the stiffness matrix using a discretization of order p
% on a 2D mesh with nelem elements in each direction.
% The type determines the type of matrix generated,
%     'ho'   : generates using high-order shape functions
%     'lin'  : approximates ho using N*elem linear elements
%     'gll'  : same as linear, except that the nodes are spaced
%              at GLL coordinates

order = p;

% geometry and mesh
% fem.geom = block3(1,1,1,'base','corner','pos',[0,0,0]); 
fem1.geom = rect2(0,1,0,1);


if ( strcmpi(typ, 'lin') )
  coords = 0.0:1/(p*nelem):1.0; 
  fem1.mesh = meshmap(fem1, 'Edgelem', {1,coords,2,coords,3,coords,4,coords});
  fem = meshextrude(fem1, 'distance', 1, 'elextlayers', {coords});
  fem.shape = 1;
elseif ( strcmpi( typ, 'gll') )
  % sp = [-1.0 -0.6 0 0.6 1.0];
  coords = getGLLcoords(p, nelem);
  fem1.mesh = meshmap(fem1, 'Edgelem', {1,coords,2,coords,3,coords,4,coords});
  fem = meshextrude(fem1, 'distance', 1, 'elextlayers', {coords});
  fem.shape = 1;
elseif ( strcmpi( typ, 'ho') )
    fem1.mesh = meshmap(fem1, 'Edgelem', {1,nelem,2,nelem});
    fem = meshextrude(fem1, 'distance', 1, 'elextlayers', {nelem});
    fem.shape = p;
else 
  fprintf('Unknown type: %s', typ);
end

% meshplot(fem);

% finite element functions
fem.dim = {'u'};

fem.equ.weak = '-(ux * ux_test + uy * uy_test + uz * uz_test + 1*u_test)';

fem.bnd.r = {'u-0'};
fem.xmesh = meshextend(fem);


% K system matrix, L rhs, boundary conditions are to be incorporated
% by imposing N*U = M
[K,L,M,N] = assemble(fem);


[Kc,Lc,Null,Ud] = femlin(fem);

% U = Null * (Kc \ Lc) + Ud;

% get indices for dofs corresponding to FE functions
nodes = xmeshinfo(fem ,'out', 'nodes');
dofs = nodes.dofs';
coords = nodes.coords';

[~,idof] = sort(dofs);

crds = coords(idof,:);
sortval = crds(:,3)*1e6 + crds(:,2)*1e3 + crds(:,1);
[~,p] = sort(sortval);

crds2 = Null' * crds;
sortval = crds2(:,3)*1e6 + crds2(:,2)*1e3 + crds2(:,1);
[~,p2] = sort(sortval);

K_new = K(p,p);
N2 = Null(p,p2);
Kc_new = N2' * K_new * N2;

% write out the matrix ...
K = Kc_new;
size(K)
fname = sprintf('K_%s_%d_%d.mat', typ, order, nelem);
save(fname, 'K');

