% clear up
clear all; close all;

% polynomial order
N = 6;
% number of elements per space dimension
no_elem = 4;

% get GLL coordinates
GLLcoords = getGLLcoords(N, no_elem);

% high order mesh
fem.geom = rect2(0,1,0,1);
fem.mesh = meshmap(fem, 'Edgelem', {1,no_elem,2,no_elem});
figure; meshplot(fem);
fem.dim = {'u'};
fem.shape = N;
fem.equ.weak = '-(ux * ux_test + uy * uy_test + 1*u_test)';
fem.bnd.r = {'u - 0'};
fem.xmesh = meshextend(fem);
[K,L,M,N] = assemble(fem);

% low oder GLL mesh
fem_low.geom = rect2(0,1,0,1);
fem_low.mesh = meshmap(fem_low, 'Edgelem', {1,GLLcoords,2,GLLcoords,...
    3,GLLcoords,4,GLLcoords});
figure; meshplot(fem_low);
fem_low.dim = {'u'};
fem_low.shape = 1;
fem_low.equ.weak = '-(ux * ux_test + uy * uy_test + 1*u_test)';
fem_low.bnd.r = {'u - 0'};
fem_low.xmesh = meshextend(fem_low);
[K_low,L_low,M_low,N_low] = assemble(fem_low);

% get indices for GLL dofs corresponding to FE functions
nodes = xmeshinfo(fem_low, 'out', 'nodes');
dofs = nodes.dofs';
crds = nodes.coords';
[~,idof] = sort(dofs);
coords = crds(idof,:);

fac = 10000*no_elem;
% this is 2D
sortval = coords(:,2)*fac + coords(:,1);
[~, p] = sort(sortval);


% low order permuted system
Kl = K_low(p,p);
Ll = L_low(p,1);
Nl = N_low(:,p);

% vector of FEM coefficients and its length
fem.sol = femlin(fem);
X = fem.sol.u;
no_dofs = length(X);
no_pts = size(coords,1);

% coordinate transformation matrix uniform->GLL
R = zeros(no_pts,no_dofs);
for i = 1:no_dofs
    X(:) = 0;
    X(i) = 1;
    R(:,i) = postinterp(fem,'u', coords(p,:)', 'U', X)';
end
% Rinv transforms from GLL to uniform Lagrange basis
Rinv = inv(R);

% high-order stuff in GLL basis
KK = Rinv' * K * Rinv;
NN = N * Rinv;
LL = Rinv' * L;

% new basis
[KC,LC,NULL,UD] = femlin('in', {'K', KK, 'L', LL, 'M', M, 'N', NN});

% old basis
[K1C,L1C, NULL1, U1D] = femlin('in', {'K', K, 'L', L, 'M', M, 'N', N});

% low order basis
[KCl,LCl, NULLl, UCl] = femlin('in', {'K', Kl, 'L', Ll, 'M', M_low, 'N', Nl});

% check reordering of solution
%sol = Rinv * (NULL * (KC\LC));
%sol1 = NULL1 * (K1C\L1C);
%norm(sol-sol1)/norm(sol)

soll = (NULLl * (KCl \ LCl));
soll1(p,1) = soll;
postsurf(fem_low,'u','U',soll1);

condKC = cond(full(KC));
condKCl = cond(full(KCl));
condKCinvKCl = cond(full(KC*inv(KCl)));

fprintf('%f, %f, %f \n', condKC, condKCl, condKCinvKCl);


return;
