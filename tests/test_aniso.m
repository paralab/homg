fp = fopen('aniso-1-128.txt', 'w');

aniso = [1 2 4 8 16 32 64 128 256 512 1024];

fprintf(fp, 'aniso\tjac\tchebyshev\tssor\n');
for i=1:length(aniso)
  clear g;
  g = create_grid_hierarchy(2, 1, [8 16 32 64 128], aniso(i));
  [u, rr, it_jac]  = g.solve_pcg(500, 'jacobi', 3, g.L, g.get_u0() );
  [u, rr, it_cheb] = g.solve_pcg(500, 'chebyshev', 3, g.L, g.get_u0() );
  [u, rr, it_ssor] = g.solve_pcg(500, 'ssor', 2, g.L, g.get_u0() );
  fprintf(fp, '%d\t%d\t%d\t%d\n', aniso(i), it_jac, it_cheb, it_ssor);
end
fclose(fp);

