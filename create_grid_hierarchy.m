function grid = create_grid_hierarchy(dim, order, nelems, sp)
% function grid = create_grid_hierarchy(dim, order, nelems)
%   specify nelems as an array of sizes, preferably factors of 2
%     from coarse to fine.
%   returns handle to the finest grid.
%
% example:
%          fine_grid = create_grid_heirarchy(2, 3, [4 8 16]);
%
% @author: Hari Sundar - hari@ices.utexas.edu 
% @date  : 28-Aug-2012


% nelems = nelem_coarse * 2.^(0:num_grids-1);
num_grids = length(nelems);

disp(['Creating grid: ' num2str(1) ' of ' num2str(num_grids) ', nelem = ' num2str(nelems(1))]);
m = homg.mesh(dim, nelems(1), sp);
coarse = homg.grid(m, order);

for i=2:num_grids
  disp(['Creating grid: ' num2str(i) ' of ' num2str(num_grids) ', nelem = ' num2str(nelems(i))]);
  m = homg.mesh(dim, nelems(i), sp);
  grid = homg.grid(m, order, coarse);
  % grid.debug = 1;
  % evc = grid.get_eigenvectors();
  if ( dim==2 )
    m.set_rhs('-8*pi^2*(sin(2*pi*x) * sin(2*pi*y))');
  else
    m.set_rhs('-12*pi^2*(sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) )');
  end
  coarse = grid;
end

