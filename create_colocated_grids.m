function grid = create_colocated_grids (dim, geom, order, nelems, sp)
% function grid = create_colocated_grids (dim, geom, order, nelem, sp)


%% create linear grid heirarchy
lin_elems = order*nelems;

lin_grid = create_grid_hierarchy(dim, geom, 1, lin_elems, sp)

%% create ho grid,
m = homg.mesh(dim, geom, nelems(end), sp);
grid = homg.grid(m, order, lin_grid);


