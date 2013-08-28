function low_order_precon(nelems, xform, order)

mesh = homg.hexmesh(nelems, xform); 

[K, M]          =  mesh.assemble_poisson (order);
[K_lin, M_lin]  =  mesh.assemble_poisson_linearized (order);

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
[x1,fl1,rr1,it1,rv1] = gmres(K, rhs, [], 1e-8, maxit, K_lin);
toc

semilogy(rv1/norm(rhs),'r-o');

end