function test_convergence(dim, order, nelem)
% function test_convergence(dim, order, nelem)

disp('grid -> err');
for i=1:length(nelem)
  m = homg_mesh(dim, nelem(i));
  m.set_rhs('-8*pi^2*(sin(2*pi*x) * sin(2*pi*y))');
  [K,L,Null,Ud] = m.assemble_poisson(order);
  % solve 
  Kc = Null'*K*Null;
  Lc = Null'*L;
  U = Null *(Kc \ Lc) + Ud;
  % check error
  nn1d = nelem(i)*order + 1;
  u_sol = reshape(U, nn1d, nn1d); 
  u_exact = zeros(size(u_sol)); 
  for k=1:nn1d
    for j=1:nn1d
      u_exact(j,k) = sin(2*pi*(k-1)/(nn1d-1)) * sin(2*pi*(j-1)/(nn1d-1));
    end
  end
  e = norm(u_exact - u_sol);
  disp([num2str(nelem(i)) ' -> ' num2str(e)]);
  % clean up
  clear m K L Null Ud Kc Lc U u_sol u_exact;
end