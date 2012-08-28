function grid = create_grid_hierarchy(dim, nelem_coarse, order, num_grids)
% function grid = create_grid_hierarchy(dim, nelem_coarse, order, num_grids)

nelems = nelem_coarse * 2.^(0:num_grids-1);

disp(['Creating grid: ' num2str(1) ' of ' num2str(num_grids) ', nelem = ' num2str(nelem_coarse)]);
m = homg.mesh(dim, nelem_coarse);
coarse = homg.grid(m, order);

for i=2:num_grids
  disp(['Creating grid: ' num2str(i) ' of ' num2str(num_grids) ', nelem = ' num2str(nelems(i))]);
  m = homg.mesh(dim, nelems(i));
  grid = homg.grid(m, order, coarse);
  coarse = grid;
end

