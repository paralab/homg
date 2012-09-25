function [rres, iterations] = test_sor_relaxation(dim, order, nelems)
% function [rres, iterations] = test_sor_relaxation(dim, order, nelems)

omegas = 0.9:0.02:1.2;
iterations = zeros(size(omegas));
rres = zeros(size(omegas));

clear g;
g = create_grid_hierarchy(dim, order, nelems);
  
for i=1:length(omegas)
  g.set_sor_omega(omegas(i));
  [u, rr, iter] = g.solve(10, 'sor', 3, g.L, g.get_u0() );
  iterations(i) = iter;
  rres(i) = rr;
end
clf;
%plot(omegas, iterations);
tt = [num2str(dim) 'D - order ' num2str(order)]; 
%title(tt);
% figure,
plot(omegas, rres);
title(tt);
