classdef grid < handle
  %GRID A single grid in a multigrid heirarchy
  % compare with mgm_grid structure in mgm_multigrid.h
  
  properties
    level
    eig_max
    eig_min
    k_evec
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
    dbg_spaces
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
      
      grid.dbg_spaces = '      ';
      grid.dbg_spaces = grid.dbg_spaces(1:end-4*grid.level);

      % for i=1:(grid.level)
      %  grid.dbg_spaces = [grid.dbg_spaces  '  '];
      % end

      grid.Mesh = mesh;
      [grid.K, grid.L, grid.Null, grid.Ud] = mesh.assemble_poisson(order);
      grid.M = mesh.assemble_mass(order);
      grid.ZeroBoundary = grid.Null * grid.Null';
      grid.smoother = 'jacobi';
      grid.jacobi_omega = 2/3;
      
      if (~ isempty(grid.Coarse) )
         grid.P = grid.Coarse.Mesh.assemble_interpolation(order, mesh.coords);
         % grid.R = inv(grid.Coarse.M) * grid.P' * grid.M ; 
         grid.R = grid.P';
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
      r = grid.ZeroBoundary * (grid.K*(grid.ZeroBoundary*u+grid.Ud) - rhs) + (u - grid.ZeroBoundary*u - grid.Ud);
    end

    function u = solve(grid, num_vcyc, smoother, smooth_steps, rhs, u)
      grid.set_smoother(smoother);
      % bdy conditions ...
      u = grid.ZeroBoundary*u;
      
      r = grid.residual(rhs, u);
      disp(['Initial residual is ' num2str(norm(r))]);
      for i=1:num_vcyc
        u = grid.vcycle(smooth_steps, smooth_steps, rhs, u);
        r = grid.residual(rhs, u);
        disp('------------------------------------------');
        disp([num2str(i) ' residual is ' num2str(norm(r))]);
        disp('------------------------------------------');
      end
    end

    % main v-cycle
    function u = vcycle(grid, v1, v2, rhs, u)
      % function u = vcycle(grid, v1, v2, rhs, u)
      % solve system using initial guess u, given rhs
      % with v1 pre and v2 post-smoothing steps
      
      disp( [grid.dbg_spaces 'v-cycle at level ' num2str(grid.level)]);
      % handle for the coarsest level
      if ( isempty( grid.Coarse ) )
        % disp('------- coarse solve -------');
        % disp( [' grid level is ' num2str(grid.level)]);
        Kc = grid.Null' * grid.K * grid.Null;
        Lc = grid.Null' * rhs;
        u = grid.Null * (Kc \ Lc) + grid.Ud;
        return;
      end
      
      % 1. pre-smooth
      clf;
      grid.plot_spectrum(u, 'k', rhs); hold on;
      u = grid.smooth ( v1, rhs, u );
      grid.plot_spectrum(u, 'b', rhs);
      
      % 2. compute residual
      res = grid.residual(rhs, u);
      
      % 3. restrict
      res_coarse = grid.R * res;
      
      % 4. recurse
      u_corr_coarse = grid.Coarse.vcycle(v1, v2, res_coarse, zeros(size(res_coarse)));
      
      % 5. prolong and correct
      u = u - grid.P * u_corr_coarse;
      grid.plot_spectrum(u, 'r', rhs);
      % 6. post-smooth
      u = grid.smooth ( v2, rhs, u );
      grid.plot_spectrum(u, 'g', rhs);
      
    end % v-cycle
    
    % smoothers
    function u = smooth (grid, v, rhs, u)
      switch(grid.smoother)
        case 'jacobi', 
          u = grid.smoother_jacobi(v, rhs, u); 
          return;
        case 'chebyshev', 
          u = grid.smoother_chebyshev(v, rhs, u); 
          return;
        case '2sr', 
          u = grid.smoother_2sr(v, rhs, u); 
          return;
        case 'hybrid',
          u = grid.smoother_hybrid(v, rhs, u);
          return;
        otherwise
          disp('ERROR: Unrecognized smoother type'); 
          return;
      end
    end

    function set_smoother(grid, sm)
      grid.smoother = sm;
      if (~ isempty(grid.Coarse) )
        grid.Coarse.set_smoother(sm);
      end
    end

    function u = smoother_jacobi (grid, v, rhs, u)
      % standard jacobi smoother
      if ( isempty(grid.jacobi_invdiag) )
        Kc = (eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
        D = diag(Kc);
        grid.jacobi_invdiag = 1./D;
      end
      
      for i=1:v
        res  = grid.jacobi_invdiag .* grid.residual(rhs, u);
        u = u - grid.jacobi_omega.*res;
        r = norm(res);
        % disp([grid.dbg_spaces 'residual: ' num2str(r)]); 
        % norm(r)
      end
    end % jacobi
    
    function u = smoother_hybrid (grid, v, rhs, u)
      u = grid.smoother_chebyshev(v, rhs, u);
    end
    
    function u = smoother_2sr (grid, v, rhs, u)
      % 2-step stationary iterative smoother
      % factors
      if ( isempty ( grid.eig_max ) )
        % Kc = grid.Null' * grid.K * grid.Null;
        Kc = (eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
        grid.eig_max = eigs(Kc,1, 'LM');  
        grid.eig_min = eigs(Kc,1, 'SM');  
      end

      l_max = grid.eig_max;
      % l_min = grid.eig_max*.9; 
      l_min = (grid.eig_min + grid.eig_max)/2;
      
      rho       = (1 - l_min/l_max)/(1 + l_min/l_max);
      alpha     = 2/( 1 + sqrt(1-rho*rho));
      epsilon   = 2/(l_min + l_max);
      epsalpha  = epsilon * alpha;
      
      % variables
      u0  = u; %zeros(size(u));
      res = grid.residual(rhs, u);
      u1 = epsilon * res ;
      
      for iter = 1:v,                            % begin iteration
        res = grid.residual ( rhs, u );
 
        u = alpha*u1 + (1-alpha)*u0 - epsalpha*res;
        u0 = u1; u1 = u;
        r = norm(res);
        % n0 = norm(u0);
        % n1 = norm(u1);
        disp([grid.dbg_spaces 'residual: ' num2str(r)]); % ' u0: ' num2str(n0) ' u1: ' num2str(n1)]);
      end % end iteration
      
    end % 2sr
    
    function u = smoother_chebyshev (grid, v, rhs, u)
      if ( isempty ( grid.eig_max ) )
        % Kc = grid.Null' * grid.K * grid.Null;
        Kc = (eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
        grid.eig_max = eigs(Kc, 1, 'LM');  
        grid.eig_min = eigs(Kc, 1, 'SM');  
      end
      
      % adjust the eigenvalues to hit the upper spectrum
      l_max = grid.eig_max*1.1;
      l_min = grid.eig_max*0.35; 
      %l_min =  (grid.eig_min + grid.eig_max)/2;
      
      c = (l_min - l_max)/2;
      d = (l_min + l_max)/2;

      p = zeros(size(u));

      for iter = 1:v
        res = grid.residual ( rhs, u ); 
        r = norm(res);
        % disp([grid.dbg_spaces 'residual: ' num2str(r)]); 
        if ( iter == 1 )
          alpha = 1.0/d;
        elseif (iter == 2)
          alpha = 2*d / (2*d*d - c*c);
        else
          alpha = 1.0/(d - alpha*c*c*0.25);
        end

        beta = alpha * d - 1.0;

        p = -alpha * res + beta * p;
        u = u + p;  
      end
    end % chebyshev

    function evec = get_eigenvectors(grid)
      % generate the correct matrix 
      Kc = (eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
      [evec, ~] = eig(full(Kc));
      grid.k_evec = evec;
    end
    
    function plot_spectrum(grid, u, clr, rhs)
      subplot(1,2,1);
      q = repmat(u,size(u'));
      b = abs(dot (grid.k_evec, q));
      % plot eigenvalues 
      plot(b, clr); hold on;
      subplot(1,2,2); 
      rr = grid.residual(rhs, u);
      n = sqrt(length(rr));
      imagesc(reshape(rr, n, n)); colorbar; hold off;
    end

    function u0 = get_u0(grid)
      n = size(grid.k_evec, 1);
      lam = ones(n,1);
      u0 = grid.k_evec*lam;
    end
    
  end %methods
  
end %classdef

