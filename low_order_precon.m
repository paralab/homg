% Example usage: low_order_precon([16 16], @homg.xform.identity, 3)

function low_order_precon(nelems, xform, order)

mesh = homg.hexmesh(nelems, xform); 

[K, M]          =  mesh.assemble_poisson (order);
[K_lin, M_lin]  =  mesh.assemble_poisson_linearized (order);

% permutation to reduce fill in when computing factorization
per = symamd(K_lin);

K_lin_chol = chol(K_lin(per,per));


bdy = mesh.get_boundary_node_indices(order);

syms x y z
if ( mesh.dim==2 )
  fx = matlabFunction(-8*pi^2*(sin(2*pi*x) * sin(2*pi*y)));
else
  fx =matlabFunction(-12*pi^2*(sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) ));
end

rhs       = mesh.assemble_rhs(fx, order);

n = size(K, 1);
x_gt = rand(n,1);
rhs = K*x_gt;

rhs (bdy) = 0;

% now solve and test 
maxit = 200;


tic
[x0,fl0,rr0,it0,rv0] = gmres(K, rhs, [], 1e-8, maxit);
toc

%figure;
semilogy(rv0/norm(rhs),'-o');
xlabel('Iteration number');
ylabel('Relative residual'); hold on;


tic
% solve reordered system and revert ordering in solution
[x1,fl1,rr1,it1,rv1] = gmres(K(per,per), rhs(per), [], 1e-8, maxit, K_lin_chol', K_lin_chol);
x1(per) = x1;
toc

fprintf('Difference between solutions: %g\n', norm(x1-x0,'fro')/norm(x0,'fro'));

semilogy(rv1/norm(rhs),'r-o');

end