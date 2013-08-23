function test_hex_conv(dim, order, nelem)
% function test_convergence(dim, order, nelem)
% nelem can be an array of grid sizes.
% 
% e.g.  test_convergence(2, 3, [8 16 32 64]); 

syms x y z

% disp('grid -> err');
fprintf('.-------------------------------------------.\n');
fprintf('| nelem\t|    err (M-norm)\t|    fac    |\n');
fprintf('|-------------------------------------------|\n');
for i=1:length(nelem)
  m = homg.hexmesh(repmat(nelem(i), 1, dim), @homg.xform.identity);
    
  if ( dim==2 )
    fx = matlabFunction(-8*pi^2*(sin(2*pi*x) * sin(2*pi*y)));
  else  
    fx =matlabFunction(-12*pi^2*(sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) ));
  end
  
  bdy = m.get_boundary_node_indices(order);
  
  K   = m.assemble_stiffness(order);
  M   = m.assemble_mass(order);
  f   = m.evaluate(fx, order);
  f(bdy) = 0;
  rhs = - M * f;
  % rhs = -( m.assemble_rhs(fx, order) );
  
  N = size(K,1);
  
  % zero-Dirichlet boundary conditions 
  K(bdy,:)    = 0;
  K(:,bdy)    = 0;
  K((N+1)*(bdy-1)+1) = 1;
  
  % rhs(bdy) = 0;
  % disp(['Solving: ' num2str(i)])
  U = K \ rhs;
  
  % check error
  if (dim == 2)
    u_fx = matlabFunction( sin(2*pi*x) * sin(2*pi*y) );
  else
    u_fx = matlabFunction( sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) );
  end % dim == 2
  u_exact = m.evaluate(u_fx, order);
  err = u_exact - U;
  e = sqrt(err' * M * err );
  if (i==1)
    fprintf('| %d\t|\t%g\t|     -     |\n', nelem(i), e);
  else
    fprintf('| %d\t|\t%g\t|  %g  |\n', nelem(i), e, prev_e/e);
  end
  prev_e = e;
  % clean up
  clear m K L Null Ud Kc Lc U u_sol u_exact err;
end
fprintf('`-------------------------------------------.\n');
