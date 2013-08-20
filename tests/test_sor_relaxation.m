function [rres, iterations] = test_sor_relaxation(dim, order, nelems)
% function [rres, iterations] = test_sor_relaxation(dim, order, nelems)

omegas = 0.5:0.01:1.5;
iterations = zeros(size(omegas));
rres = zeros(size(omegas));

clear g;
g = create_grid_hierarchy(dim, order, nelems, 1);
  
for i=1:length(omegas)
  g.set_sor_omega(omegas(i));
  [u, rr, iter] = g.solve_pcg(100, 'ssor', 2, g.L, g.get_u0() );
  iterations(i) = iter;
  rres(i) = rr;
end
clf;
plot(omegas, iterations);
tt = [num2str(dim) 'D - order ' num2str(order)]; 
title(tt);
figure,
plot(omegas, rres);
title(tt);
