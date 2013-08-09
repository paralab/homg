function test_convergence(dim, order, nelem)
% function test_convergence(dim, order, nelem)
% nelem can be an array of grid sizes.
% 
% e.g.  test_convergence(2, 3, [8 16 32 64]); 

% disp('grid -> err');
fprintf('.-------------------------------------------.\n');
fprintf('| nelem\t|    err (M-norm)\t|    fac    |\n');
fprintf('|-------------------------------------------|\n');
for i=1:length(nelem)
  m = homg.mesh(dim, nelem(i),1);
  if ( dim==2 )
    m.set_rhs('-8*pi^2*(sin(2*pi*x) * sin(2*pi*y))');
  else
    m.set_rhs('-12*pi^2*(sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) )');
  end
  [K,L,Null,Ud] = m.assemble_poisson(order);
  M = m.assemble_mass(order);
  % solve
  Kc = Null'*K*Null;
  Lc = Null'*L;
  U = Null *(Kc \ Lc) + Ud;
  % check error
  nn1d = nelem(i)*order + 1;
  if (dim == 2)
    u_sol = reshape(U, nn1d, nn1d);
    u_exact = zeros(size(u_sol));
    for k=1:nn1d
      for j=1:nn1d
        u_exact(j,k) = sin(2*pi*(k-1)/(nn1d-1)) * sin(2*pi*(j-1)/(nn1d-1));
      end
    end
  else
    u_sol = reshape(U, nn1d, nn1d, nn1d);
    u_exact = zeros(size(u_sol));
    for l=1:nn1d
      for k=1:nn1d
        for j=1:nn1d
          u_exact(j,k,l) = sin(2*pi*(k-1)/(nn1d-1)) * sin(2*pi*(j-1)/(nn1d-1)) * sin(2*pi*(l-1)/(nn1d-1));
        end
      end
    end
  end % dim == 2
  err = u_exact(:) - u_sol(:);
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
