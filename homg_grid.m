classdef homg_grid
  %HOMG_GRID grid class for multigrid
  % compare with mgm_grid structure in mgm_multigrid.h 
  
  properties
    level
    eig_max
    eig_min
    K
    M
    R
    P
    Coarse
  end % properties
  
  methods
    function grid = homg_grid(mesh, coarse) 
      if ((nargin == 1) || isempty(coarse))
        grid.level = 0;
        grid.coarse = [];
      else
        grid.level = coarse.level + 1;
        grid.coarse = coarse;
      end
      
      % grid.K = mesh.assemble_poisson();
      % grid.M = M;
      % grid.R = R;
      % grid.P = P;
    end
    
  end %methods

end %classdef

