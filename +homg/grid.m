classdef grid < handle
  %GRID A single grid in a multigrid heirarchy
  % compare with mgm_grid structure in mgm_multigrid.h
  
  properties
    level
    eig_max
    eig_min
    jacobi_omega
    jacobi_invdiag
    smoother
    K
    L
    Null
    ZeroBoundary
    Ud
    M
    R
    P
    Mesh
    Coarse  % handle to coarse grid 
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
      
      grid.Mesh = mesh;
      [grid.K, grid.L, grid.Null, grid.Ud] = mesh.assemble_poisson(order);
      grid.M = mesh.assemble_mass(order);
      grid.ZeroBoundary = grid.Null * grid.Null';
      grid.smoother = 'jacobi';
      % fixme
      if (~ isempty(grid.Coarse) )
         grid.P = grid.Coarse.Mesh.assemble_interpolation(order, mesh.coords);
      end
    end
    
    % compute the residual
    function r = residual(grid, rhs, u)
    % function r = residual(grid, u, rhs)
      if ( nargin < 2 )
        rhs = grid.L;
      end
      if ( nargin < 3 )
        u = zeros(size(rhs));
      end

      % r = Au - f 
      r = grid.ZeroBoundary * (grid.K*(grid.ZeroBoundary*u+grid.Ud) - grid.L) + (u - grid.ZeroBoundary*u - grid.Ud);
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
        u = grid.Null * (Kc \ Lc) + grid.Ud;
        return;
      end
      
      % 1. pre-smooth
      u = grid.smooth ( v1, rhs, u );
      
      % 2. compute residual
      res = grid.residual(rhs, u);
      
      % 3. restrict
      res_coarse = grid.R * res;
      
      % 4. recurse
      u_corr_coarse = grid.coarse.vcycle(v1, v2, res_coarse, zeros(size(res_coarse)));
      
      % 5. prolong and correct
      u = u + grid.P * u_corr_coarse;
      
      % 6. post-smooth
      u = grid.smooth ( v2, rhs, u );
      
    end % v-cycle
    
    % smoothers
    function u = smooth (grid, v, rhs, u)
      switch(grid.smoother)
        case 'jacobi', return grid.smoother_jacobi(v, rhs, u);
        case 'chebyshev', return grid.smoother_chebyshev(v, rhs, u);
        case '2sr', return grid.smoother_2sr(v, rhs, u);
        otherwise disp('ERROR: Unrecognized smoother type'), return;
      end
    end

    function set_smoother(grid, sm)
      grid.smoother = sm;
    end

    function u = smoother_jacobi (grid, v, rhs, u)
      % standard jacobi smoother
      if ( isempty(grid.jacobi_invdiag) )
        D = diag(grid.K);
        grid.jacobi_invdiag = 1./D;
      end
      
      for i=1:v
        r  = grid.jacobi_invdiag .* grid.residual(rhs, u);
        u = u + omega*r;
      end
    end % jacobi
    
    function u = smoother_2sr (grid, v, rhs, u)
      % 2-step stationary iterative smoother
      % factors
      if ( isempty ( grid.eig_max ) )
        Kc = grid.Null' * grid.K * grid.Null;
        grid.eig_max = eigs(Kc,1, 'LM');  
        grid.eig_min = eigs(Kc,1, 'SM');  
      end

      l_max = grid.eig_max;
      l_min = grid.eig_min;
      
      rho       = (1 - l_min/l_max)/(1 + l_min/l_max);
      alpha     = 2/( 1 + sqrt(1-rho*rho));
      epsilon   = 2/(l_min + l_max);
      epsalpha  = epsilon * alpha;
      
      % variables
      res  = -rhs ?
      u0 = zeros(size(u));
      u1 = epsilon * r ;
      
      for iter = 1:v,                            % begin iteration
        res = grid.residual ( rhs, u );
 
        u = alpha*u1 + (1-alpha)*u0 - epsalpha*res;
        u0 = u1; u1 = u;
      end % end iteration
      
    end % 2sr
    
    function u = smoother_chebyshev (grid, v, rhs, u)
      if ( isempty ( grid.eig_max ) )
        Kc = grid.Null' * grid.K * grid.Null;
        grid.eig_max = eigs(Kc,1, 'LM');  
        grid.eig_min = eigs(Kc,1, 'SM');  
      end
      
      % FIXME adjust the eigenvalues to hit the upper spectrum
      l_max = grid.eig_max;
      l_min = grid.eig_min;
      
      c = (l_min - l_max)/2;
      d = (l_min + l_max)/2;

      p = zeros(size(u));

      for iter = 1:v
        res = grid.residual ( rhs, u ); 

        if ( iter == 1 )
          alpha = 1.0/d;
        else if (iter == 2)
          alpha = 2*d / (2*d*d - c*c);
        else
          alpha = 1.0/(d - alpha*c*c*0.25);
        end

        beta = alpha * d - 1.0;

        p = alpha * res + beta * p;
        u = u + p;  
      end
    end % chebyshev

  end %methods
  
end %classdef

