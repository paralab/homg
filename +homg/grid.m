classdef grid < handle
  %GRID A single grid in a multigrid heirarchy
  % compare with mgm_grid structure in mgm_multigrid.h
  
  properties
    level
    eig_max
    eig_min
    jacobi_omega
    jacobi_invdiag
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
        grid.Coarse = [];
      else
        grid.level = coarse.level + 1;
        grid.Coarse = coarse;
      end
      
      [grid.K, grid.L, grid.Null, grid.Ud] = mesh.assemble_poisson(order);
      grid.M = mesh.assemble_mass(order);
      
      % fixme
      if (~ isempty(grid.Coarse) )
        grid.P = mesh.assemble_interpolation(order, mesh.coords);
      end
    end
 
    % compute the residual
    function r = residual(grid, rhs, u)
    % function r = residual(grid, u, rhs)
      r = rhs - grid.K*u;  
    
      Kc = grid.Null' * grid.K * grid.Null;
      Lc = grid.Null' * rhs;
      
      % r = ƒ - Au 
      r = Lc - (Kc * grid.Null'(u - grid.Ud));
    end

    % main v-cycle
    function u = vcycle(grid, v1, v2, rhs, u)
      % function u = vcycle(grid, v1, v2, rhs, u)
      % solve system using initial guess u, given rhs
      % with v1 pre and v2 post-smoothing steps
      
      % handle for the coarsest level
      if ( isempty( grid.Coarse ) )
        Kc = grid.Null' * grid.K * grid.Null;
        Lc = grid.Null' * rhs;
        u = grid.Null *(Kc \ Lc) + grid.Ud;
        return;
      end
      
      Kc = grid.Null' * grid.K * grid.Null;
      Lc = grid.Null' * rhs;
      
      % 1. pre-smooth
      u = grid.smooth ( v1, rhs, u );
      
      % 2. compute residual
      res = Lc - Kc * grid.Null' * (u - grid.Ud);
      
      % 3. restrict
      res_coarse = grid.R * res;
      
      % 4. recurse
      u_corr_coarse = grid.coarse.vcycle(v1, v2, res_coarse, zeros(size(res_coarse)));
      
      % 5. prolong and correct
      u = u - grid.P * u_corr_coarse;
      
      % 6. post-smooth
      u = grid.smooth ( v2, rhs, u );
      
    end % v-cycle
    
    % smoothers
    
    function u = smoother_jacobi (grid, v, rhs, u)
      % standard jacobi smoother
      if ( isempty(grid.jacobi_invdiag) )
        D = diag(grid.K);
        grid.jacobi_invdiag = 1./D;
      end
      
      Kc = grid.Null' * grid.K * grid.Null;
      Lc = grid.Null' * rhs;
      
      for i=1:v
        r  = grid.jacobi_invdiag .* (Lc - Kc * u);
        u = u + omega*r;
      end
    end
    
    function u = smoother_2sr (grid, v, rhs, u)
      % 2-step stationary iterative smoother
      % factors
      l_max = grid.eig_max;
      l_min = grid.eig_min;
      rho       = (1 - l_min/l_max)/(1 + l_min/l_max);
      alpha     = 2/( 1 + sqrt(1-rho*rho));
      epsilon   = 2/(l_min + l_max);
      epsalpha  = epsilon * alpha;
      
      Kc = grid.Null' * grid.K * grid.Null;
      Lc = grid.Null' * rhs;
      
      % variables
      r = -Lc; % or rhs ?
      d0 = zeros(size(u));
      d1 = epsilon * r ;
      
      for iter = 1:v,                            % begin iteration
        r = Kc*d1 - Lc;
        d = alpha*d1 + (1-alpha)*d0 - epsalpha*r;
        d0 = d1; d1 = d;
      end % end iteration
      
    end
    
  end %methods
  
end %classdef

