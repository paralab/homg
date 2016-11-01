function grid = create_hdg_grids(dim, xform, orders, nelems)
% function grid = create_hdg_grids(dim, xform, orders, nelems)
%
%   dim  is the dimension of the mesh, either 2 or 3
%   geom is the domain geometry, currently 'box', 'fan'
%   specify nelems as an array of sizes, preferably factors of 2
%     from coarse to fine.
%   specify order as the desired hdg order at the finest level  
%   returns handle to the finest grid.
%
% example:
%          fine_grid = create_hdg_grids(2, @homg.xform.identity, [1 2 4], [8 16]);
%
% @author: Hari Sundar - hari@ices.utexas.edu 
% @date  : 28-Aug-2012, 31-Apr-2013, 7/20/2014

% h grids are CG
% p grids are hDG

num_hgrids = length(nelems); % also number of CG grids
num_pgrids = length(orders); % also number of hDG grids

num_grids = num_hgrids + num_pgrids;

% disp('Creating linear CG h-grids first');

disp(['Creating CG-grid: ' num2str(1) ' of ' num2str(num_grids) ', order = ' num2str(orders(1)) ', nelem = ' num2str(nelems(1))]);

m = homg.hexmesh(repmat(nelems(1), 1, dim), xform);
coarse = homg.grid(m, orders(1));
grid = coarse;
% disp('---- created grid ----')

for i=2:num_hgrids
  disp(['Creating CG-grid: ' num2str(i) ' of ' num2str(num_grids) ', order = ' num2str(orders(1)) ', nelem = ' num2str(nelems(i))]);
  m = homg.hexmesh(repmat(nelems(i), 1, dim), xform);
  grid = homg.grid(m, orders(1), coarse);
  coarse = grid;
end

% assemble poisson for the CG grids 
mu = @(x,y)(1);
grid.assemble_poisson(mu);

hfine = nelems(num_hgrids);

% disp('Creating p-grids now');
for i=1:num_pgrids
  disp(['Creating hDG-grid: ' num2str(i+num_hgrids) ' of ' num2str(num_grids) ', order = ' num2str(orders(i)) ', nelem = ' num2str(hfine)]);
  m = homg.hexmesh(repmat(hfine, 1, dim), xform);
  grid = homg.grid(m, orders(i), coarse);
  grid.gen_hdg_matrix();
%  grid.skel_to_cg_matrix();
%  grid.cg_to_skel_matrix();
  coarse = grid;
end

grid.is_finest = true;

