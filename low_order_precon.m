% Example usage: low_order_precon([16 16], @homg.xform.identity, 3)
function [it_gll, it_smooth_1, it_smooth_3] = low_order_precon(nelems, xform, order)

mesh = homg.hexmesh(nelems, xform); 

if (mesh.dim == 2)
  mu = @(xx,yy)(1 + 1000000*( (cos(2*pi*xx))^2 + (cos(2*pi*yy))^2 ));
else
  mu = @(xx,yy,zz)(1 + 1000000*( (cos(2*pi*xx))^2 + (cos(2*pi*yy))^2 + (cos(2*pi*zz))^2 ));
end

mesh.set_coeff(mu);

[K, M]          =  mesh.assemble_poisson (order);
[K_lin, M_lin]  =  mesh.assemble_poisson_linearized (order);

% grid = homg.grid(mesh, order);
% grid.assemble_poisson(mu);

% muni = homg.hexmesh(nelems*order, xform); 
% muni.set_coeff(mu);
% 
% [K_uni, M_uni]  =  muni.assemble_poisson_linearized (1);

grid_lin = homg.grid(mesh, order);
grid_lin.assemble_poisson(mu);
grid_lin.use_linearized_smoothers();
grid_lin.is_finest = true;

bdy     = mesh.get_boundary_node_indices(order);

% syms x,y,z;
% if ( mesh.dim==2 )
%   fx = matlabFunction(-8*pi^2*(sin(2*pi*x) * sin(2*pi*y)));
% else
%   fx = matlabFunction(-12*pi^2*(sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) ));
% end

% rhs       = mesh.assemble_rhs(fx, order);

n = size(K, 1);
x_gt = rand(n,1);
rhs = K*x_gt;

rhs (bdy) = 0;

% now solve and test 
maxit = min(350, size(K,1));


%tic
%[x0,fl0,rr0,it0,rv0] = gmres(K, rhs, [], 1e-8, maxit);
%toc

%figure;
%semilogy(rv0/norm(rhs),'-o');
%xlabel('Iteration number');
%ylabel('Relative residual'); hold on;


%tic
% permutation to reduce fill in when computing factorization
per = symamd(K_lin);
K_lin_chol = chol(K_lin(per,per));

% solve reordered system and revert ordering in solution
[x1,fl1,rr1,it1,rv1] = pcg(K(per,per), rhs(per), 1e-8, maxit, K_lin_chol', K_lin_chol,[]);
% x1(per) = x1;
%toc

%tic
% per_uni = symamd(K_uni);
% K_uni_chol = chol(K_uni(per_uni,per_uni));
% [x2,fl2,rr2,it2,rv2] = pcg(K(per_uni,per_uni), rhs(per_uni), 1e-8, maxit, K_uni_chol', K_uni_chol,[]);
% x2(per_uni) = x2;
%toc


num_smooth = 1;
[x2,fl2,rr2,it2,rv2] = pcg(K, rhs, 1e-8, maxit, @smooth_gll);

function yo = smooth_gll(xi) 
	% first smooth 
	xs = grid_lin.smoother_chebyshev (num_smooth, xi, zeros(size(xi)));
	
  ys = xs - (K_lin \ grid_lin.residual(xi, xs));
	
	yo = grid_lin.smoother_chebyshev (num_smooth, xi, ys);
end

num_smooth = 3;
[x3,fl3,rr3,it3,rv3] = pcg(K, rhs, 1e-8, maxit, @smooth_gll);


%fprintf('Difference between solutions: %g\n', norm(x1-x0,'fro')/norm(x0,'fro'));

% disp(['order: ' num2str(order) ' -- iterations: ' num2str(it1(2)) ' --- uniform : ' num2str(it2(2))]);

it_gll 				= it1;
it_smooth_1 	= it2;
it_smooth_3 	= it3;

%semilogy(rv1/norm(rhs),'r-o');
%hold off;

end
