classdef grid < handle
  %GRID A single grid in a multigrid heirarchy
  % compare with mgm_grid structure in mgm_multigrid.h 
  
  properties
    level
    eig_max
    eig_min
    K
    L
    Null
    Ud
    M
    R
    P
    Coarse
  end % properties
  
  methods
    function grid = grid(mesh, order, coarse) 
      if ((nargin < 3) || isempty(coarse))
        grid.level = 0;
        grid.coarse = [];
      else
        grid.level = coarse.level + 1;
        grid.coarse = coarse;
      end
      
      [grid.K, grid.L, grid.Null, grid.Ud] = mesh.assemble_poisson(order);
      grid.M = mesh.assemble_mass(order);
      
      if (~ isempty(grid.coarse) )
        grid.P = mesh.assemble_interpolation(order, mesh.coords);
      end
    end
    
  end %methods

end %classdef

