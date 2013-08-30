function grid = create_hexmesh_grids(dim, mu, xform, orders, nelems)
% function grid = create_hexmesh_grids(dim, mu, xform, orders, nelems)
%
%   dim  is the dimension of the mesh, either 2 or 3
%   geom is the domain geometry, currently 'box', 'fan'
%   specify nelems as an array of sizes, preferably factors of 2
%     from coarse to fine.
%   specify orders as an array of sizes, as multiples of 2, and 
%     going down to 1
%   sp is the aspect ratio of the element, for stretching the mesh
%    
%   returns handle to the finest grid.
%
% example:
%          fine_grid = create_grid_heirarchy(2, 'box', [1 2 4], [4 8 16], 1);
%
% @author: Hari Sundar - hari@ices.utexas.edu 
% @date  : 28-Aug-2012, 31-Apr-2013

% nelems = nelem_coarse * 2.^(0:num_grids-1);
num_hgrids = length(nelems);
num_pgrids = length(orders);

num_grids = num_hgrids + num_pgrids - 1;

% disp('Creating linear h-grids first');

%disp(['Creating h-grid: ' num2str(1) ' of ' num2str(num_grids) ', order = ' num2str(orders(1)) ', nelem = ' num2str(nelems(1))]);

m = homg.hexmesh(repmat(nelems(1), 1, dim), xform);
coarse = homg.grid(m, orders(1), mu );
%disp('---- created grid ----')

for i=2:num_hgrids
  %disp(['Creating h-grid: ' num2str(i) ' of ' num2str(num_grids) ', order = ' num2str(orders(1)) ', nelem = ' num2str(nelems(i))]);
  m = homg.hexmesh(repmat(nelems(i), 1, dim), xform);
  %%% disp('---- created mesh')
  % grid.debug = 1;
  % evc = grid.get_eigenvectors();
  %if ( dim==2 )
    %m.set_rhs('0');
    % m.set_rhs('(1 - 8*pi^2)*(sin(2*pi*x) * sin(2*pi*y))');
    %m.set_coeff('1');
    % m.set_coeff('1 + 1000000*((cos(2*pi*x))^2 + (cos(2*pi*y))^2 )');
  %else
    % m.set_rhs('(1 - 12*pi^2)*(sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) )');
    % m.set_coeff('1');
    % m.set_coeff('1 + 1000000*((cos(2*pi*x))^2 + (cos(2*pi*y))^2 + (cos(2*pi*z))^2)');
  %end
  grid = homg.grid(m, orders(1), mu, coarse);
  %disp('---- created grid ----')
  
  coarse = grid;
end

hfine = nelems(num_hgrids);

% disp('Creating p-grids now');
for i=2:num_pgrids
  %disp(['Creating p-grid: ' num2str(i+num_hgrids-1) ' of ' num2str(num_grids) ', order = ' num2str(orders(i)) ', nelem = ' num2str(hfine)]);
  % m = homg.mesh(dim, geom, hfine, sp);
  m = homg.hexmesh(repmat(hfine, 1, dim), xform);
  grid = homg.grid(m, orders(i), mu, coarse);
%   if ( dim == 2 )
%     m.set_rhs('(1 - 8*pi^2)*(sin(2*pi*x) * sin(2*pi*y))');
%     m.set_coeff('1');
%   else
%     m.set_rhs('(1 - 12*pi^2)*(sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) )');
%     m.set_coeff('1');
%   end
  coarse = grid;
end

