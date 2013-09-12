classdef grid < handle
  %GRID A single grid in a multigrid heirarchy
  % compare with mgm_grid structure in mgm_multigrid.h
  
  properties
    level
    eig_max
    eig_min
    k_evec
    k_lam
    jacobi_omega
    jacobi_invdiag
    gs_G
    gs_c
    ssor_M
    ssor_N
    sor_G
    sor_c
    sor_omega
    smoother
    K
    L
    % Null
    % Zero
    Boundary
    Ud
    M
    R
    P
    Mesh
    Coarse  % handle to coarse grid 
    debug
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
      grid.debug = 0;
      
      grid.Mesh = mesh;
      
			mesh.set_order(order); 
      grid.sor_omega = 1;
      if (~ isempty(grid.Coarse) )
         grid.P = grid.Coarse.Mesh.assemble_interpolation(order);
         grid.R = grid.P';
      end
			
			% for Dirichlet boundary conditions
			grid.Boundary = mesh.get_boundary_node_indices(order);

			%% defaults ...
		  grid.smoother = 'sor';
      grid.jacobi_omega = 2/3;
    end
    
		function assemble_poisson(grid, mu)
			% fine grid material props ...
			if isnumeric(mu)
				grid.Mesh.set_muvec (mu) ;
			else
				grid.Mesh.set_coeff (mu) ;
			end
			% assemble for this level ... 
      [grid.K, grid.M] = grid.Mesh.assemble_poisson(grid.Mesh.order);
			syms x y z
			if ( grid.Mesh.dim == 2 )
				fx = matlabFunction(-8*pi^2*(sin(2*pi*x) * sin(2*pi*y)));
			else
				fx =matlabFunction(-12*pi^2*(sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) ));
			end
			
			grid.L = grid.Mesh.assemble_rhs(fx, grid.Mesh.order);
			grid.L(grid.Boundary) = 0;
      
      % propagate to lower grids
      if (~ isempty(grid.Coarse) )
        if isnumeric(mu)
          % M = grid.Coarse.Mesh.assemble_mass(1);
          % mu_coarse =   M \ ( grid.R * grid.M * mu );
          % max(mu)
          % max(mu_coarse)
          harmonic = 0;
          if (grid.Mesh.dim == 2)
            mu2 = reshape(mu, grid.Mesh.nelems(1), grid.Mesh.nelems(2));
            if (harmonic)
              mu_coarse = 4 ./ ( 1./mu2(1:2:end, 1:2:end) + 1./mu2(2:2:end, 1:2:end) + 1./mu2(1:2:end, 2:2:end) + 1./mu2(2:2:end, 2:2:end) );
            else
              mu_coarse = 0.25*(mu2(1:2:end, 1:2:end) + mu2(2:2:end, 1:2:end) + mu2(1:2:end, 2:2:end) + mu2(2:2:end, 2:2:end));
            end
          else
            mu3 = reshape(mu, grid.Mesh.nelems(1), grid.Mesh.nelems(2), grid.Mesh.nelems(3));
            if (harmonic)
              mu_coarse = 8 ./ ( 1./mu3(1:2:end, 1:2:end, 1:2:end) + 1./mu3(2:2:end, 1:2:end, 1:2:end) ...
                               + 1./mu3(1:2:end, 2:2:end, 1:2:end) + 1./mu3(2:2:end, 2:2:end, 1:2:end) ...
                               + 1./mu3(1:2:end, 1:2:end, 2:2:end) + 1./mu3(2:2:end, 1:2:end, 2:2:end) ...
                               + 1./mu3(1:2:end, 2:2:end, 2:2:end) + 1./mu3(2:2:end, 2:2:end, 2:2:end) );
            else
              mu_coarse = 0.125*(mu3(1:2:end, 1:2:end, 1:2:end) + mu3(2:2:end, 1:2:end, 1:2:end) + mu3(1:2:end, 2:2:end, 1:2:end) + mu3(2:2:end, 2:2:end, 1:2:end) + ...
                mu3(1:2:end, 1:2:end, 2:2:end) + mu3(2:2:end, 1:2:end, 2:2:end) + mu3(1:2:end, 2:2:end, 2:2:end) + mu3(2:2:end, 2:2:end, 2:2:end) );
            end
          end
          grid.Coarse.assemble_poisson (mu_coarse(:)) ;
        else
          grid.Coarse.assemble_poisson (mu) ;
        end
      end
      
    end
    
		function set_stiffness(grid, K)
			grid.K = K;
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
      % r = grid.ZeroBoundary * (grid.K*(grid.ZeroBoundary*u+grid.Ud) - rhs) + (u - grid.ZeroBoundary*u - grid.Ud);
      % u(grid.Boundary) = 0;
      r = grid.K*u - rhs;
      % r(grid.Boundary) = 0;
    end


    function [u, rr, iter] = solve_lin_pcg(grid, num_vcyc, smoother, smooth_steps, rhs, u)
      % disp('setting smoother');
      grid.set_smoother(smoother);
      
      % disp('computing initial residual');
      r = grid.residual(rhs, u);
      rho = zeros(size(u));
      % disp('outer v-cycle');
      rho = grid.Coarse.vcycle(smooth_steps, smooth_steps, r, rho);
      p = rho;
      disp(['Initial residual is ' num2str(norm(r))]);
      disp('------------------------------------------');
      r0 = norm(r);
      for i=1:num_vcyc
        % disp(['inner v-cycle: ' num2str(i)]);
        h = grid.K * p;
        rho_res = dot (rho, r);
        alpha = rho_res / dot ( p, h );
        u = u + alpha*p;
        r = r - alpha*h;

        % rho_res_prev = rho_res;
        
        disp([num2str(i, '%03d\t') ': |res| = ' num2str(norm(r),'\t%8.4e')]);
        if (norm(r)/r0 < 1e-8)
          iter = i;
          rr = norm(r)/r0;
          return;
        end
        
        % precondition ..
        rho = zeros(size(u)); % needed ?
        rho = grid.Coarse.vcycle(smooth_steps, smooth_steps, r, rho);

        beta = dot(rho, r) / rho_res ;
        p = rho + beta*p;
      end
      disp('------------------------------------------');
      iter = num_vcyc;
      rr = norm(r)/r0;
    end


    function [u, rr, iter] = solve_pcg(grid, num_vcyc, smoother, smooth_steps, rhs, u)
      % disp('setting smoother');
      grid.set_smoother(smoother);
      
      % disp('computing initial residual');
      r = grid.residual(rhs, u);
      rho = zeros(size(u));
      % disp('outer v-cycle');
      rho = grid.vcycle(smooth_steps, smooth_steps, r, rho);
      p = rho;
      % disp(['Initial residual is ' num2str(norm(r))]);
      % disp('------------------------------------------');
      r0 = norm(r);
      for i=1:num_vcyc
        % disp(['inner v-cycle: ' num2str(i)]);
        h = grid.K * p;
        rho_res = dot (rho, r);
        alpha = rho_res / dot ( p, h );
        u = u + alpha*p;
        r = r - alpha*h;

        % rho_res_prev = rho_res;
        
        % disp([num2str(i, '%03d\t') ': |res| = ' num2str(norm(r),'\t%8.4e')]);
        if (norm(r)/r0 < 1e-8)
          iter = i;
          rr = norm(r)/r0;
          return;
        end
        
        % precondition ..
        rho = zeros(size(u)); % needed ?
        rho = grid.vcycle(smooth_steps, smooth_steps, r, rho);
        
        beta = dot(rho, r) / rho_res ;
        p = rho + beta*p;
      end
      % disp('------------------------------------------');
      iter = num_vcyc;
      rr = norm(r)/r0;
    end
    
    function [u, rr, iter] = solve(grid, num_vcyc, smoother, smooth_steps, rhs, u)
      grid.set_smoother(smoother);
      % bdy conditions ...
      % u = grid.ZeroBoundary*u;
      
      r = grid.residual(rhs, u);
      % disp(['Initial residual is ' num2str(norm(r))]);
      % disp('------------------------------------------');
      r0 = norm(r);
      for i=1:num_vcyc
        u = grid.vcycle(smooth_steps, smooth_steps, rhs, u);
        r = grid.residual(rhs, u);
        % disp([num2str(i) ': |res| = ' num2str(norm(r))]);
        if (norm(r)/r0 < 1e-8)
          iter = i;
          rr = norm(r)/r0;
          return;
        end
      end
      % disp('------------------------------------------');
      iter = num_vcyc;
      rr = norm(r)/r0;
    end

    % main v-cycle
    function u = vcycle(grid, v1, v2, rhs, u)
      % function u = vcycle(grid, v1, v2, rhs, u)
      % solve system using initial guess u, given rhs
      % with v1 pre and v2 post-smoothing steps
      
      % disp( [grid.dbg_spaces 'v-cycle at level ' num2str(grid.level)]);
      % handle for the coarsest level
      if ( isempty( grid.Coarse ) )
        % disp('------- coarse solve -------');
        % disp( [' grid level is ' num2str(grid.level)]);
        % Kc = grid.Null' * grid.K * grid.Null;
        % Lc = grid.Null' * rhs;
        % u = grid.Null * (Kc \ Lc) + grid.Ud;
        u = grid.K \ rhs;
        return;
      end
      
      % 1. pre-smooth
      if (grid.debug) clf; end
      grid.plot_spectrum(u, 'k', rhs); 
      if (grid.debug) hold on; end
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
        case 'gs',
          grid.sor_omega = 1.0;
          u = grid.smoother_gauss_seidel(v, rhs, u); 
          return;
        case 'chebssor'
          u = grid.smoother_chebyshev_ssor (v, rhs, u);
          return;
        case 'chebyshev2'
          u = grid.smoother_chebyshev2 (v, rhs, u);
          return;
        case 'chebyshev', 
          u = grid.smoother_chebyshev(v, rhs, u); 
          return;
        case 'sor',
          u = grid.smoother_sor(v, rhs, u);
          return;
        case 'ssor',
          u = grid.smoother_sym_sor(v, rhs, u);
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

    function set_coeff(grid, mu)
      grid.Mesh.set_coeff (mu) ;
      if (~ isempty(grid.Coarse) )
        grid.Coarse.Mesh.set_coeff (mu) ;
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
        % Kc = (eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
        D = diag(grid.K);
        grid.jacobi_invdiag = 1./D;
      end
      
      for i=1:v
        res  = grid.jacobi_invdiag .* grid.residual(rhs, u);
        u = u - grid.jacobi_omega.*res;
        % r = norm(res);
        % disp([grid.dbg_spaces num2str(r)]); 
        % norm(r)
      end  
    end % jacobi
    
    
    function u = smoother_sor (grid, v, rhs, u)
      if ( isempty ( grid.sor_G ) )
        % Kc = (eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
        kL = tril(grid.K, -1);
        kD = diag(diag(Kc));
        kU = triu(Kc, 1);
        grid.sor_G = - (kD + grid.sor_omega*kL) \ (grid.sor_omega*kU + (grid.sor_omega - 1.0)*kD );
        grid.sor_c = grid.ZeroBoundary *( (kD + grid.sor_omega*kL) \ (grid.sor_omega*rhs) );
      end

      for i=1:v
         u = grid.sor_G*u + grid.sor_c;
      end
    end
    
    function u = smoother_sym_sor (grid, v, rhs, u)
      if ( isempty ( grid.ssor_M ) )
        w = grid.sor_omega;
        n = length(u);
        grid.ssor_M = spdiags( (1/w)*diag(grid.K), 0, n, n) + tril(grid.K,-1);
        grid.ssor_N = spdiags(((1-w)/w)*diag(grid.K), 0, n, n) - triu(grid.K,1);
      end

      for i=1:v
        r = grid.residual(rhs, u);
        u = u - grid.ssor_M \ r;
        u = grid.ssor_M' \ (grid.ssor_N'*u + rhs);
      end
    end
   
    function u = smoother_chebyshev_jacobi (grid, v, rhs, u)
      if ( isempty ( grid.eig_max ) )
        % Kc = grid.Null' * grid.K * grid.Null;
        Kc = (eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
        D = diag(Kc);
        grid.jacobi_invdiag = 1./D;
        Kc = (eye(size(Kc)) - (diag(D) \ Kc) );
        grid.eig_max = eigs(Kc, 1, 'LM');  
        grid.eig_min = eigs(Kc, 1, 'SM');  
      end
      
      l_max = grid.eig_max;
      l_min =  (grid.eig_min + grid.eig_max)/2;
      
      rho = 2/(l_min + l_max);
      
      mu_0 = 1;
      mu_1 = rho;
      y_0 = u;
      y_1 = grid.smoother_jacobi(1, rhs, u); 
      for i=2:v
        mu_2 = 1.0 / ( 2.0/(rho*mu_1) - 1.0/mu_0);
        u = (2.0*mu_2)/(rho*mu_1)* grid.smoother_jacobi(1, rhs, y_1) - (mu_2/mu_1)*y_0;
        y_0 = y_1; mu_0 = mu_1;
        y_1 =  u ; mu_1 = mu_2; 
      end
    end

    function u = smoother_chebyshev_ssor (grid, v, rhs, u)
      if ( isempty ( grid.eig_max ) )
        % Kc = grid.Null' * grid.K * grid.Null;
        Kc = (eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
        grid.eig_max = eigs(Kc, 1, 'LM');  
        grid.eig_min = eigs(Kc, 1, 'SM');  
      end
      if ( isempty ( grid.ssor_M ) )
        w = grid.sor_omega;
        n = length(u);
        grid.ssor_M = spdiags( (1/w)*diag(grid.K), 0, n, n) + tril(grid.K,-1);
        grid.ssor_N = spdiags(((1-w)/w)*diag(grid.K), 0, n, n) - triu(grid.K,1);
      end
      
      l_max = grid.eig_max;
      l_min =  (grid.eig_min + grid.eig_max)/2;
      
      rho = 2/(l_min + l_max);
      
      mu_0 = 1;
      mu_1 = rho;
      y_0 = u;
      y_1 = grid.smoother_sym_sor(1, rhs, u); 
      for i=2:v
        mu_2 = 1.0 / ( 2.0/(rho*mu_1) - 1.0/mu_0);
        u = (2.0*mu_2)/(rho*mu_1)* grid.smoother_sym_sor(1, rhs, y_1) - (mu_2/mu_1)*y_0;
        y_0 = y_1; mu_0 = mu_1;
        y_1 =  u ; mu_1 = mu_2; 
      end

    end

    function set_sor_omega(grid, w)
      grid.sor_omega = w;
      grid.sor_G = [];
      grid.sor_c = [];
      grid.ssor_M = [];
      grid.ssor_N = [];
    end
    
    function u = smoother_gauss_seidel (grid, v, rhs, u)
      if ( isempty ( grid.gs_G ) )
        Kc = (eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
        LD = tril(Kc);
        grid.gs_G = -LD \ triu(Kc, 1);
        grid.gs_c = grid.ZeroBoundary *( LD \ rhs );
      end

      for i=1:v
         u = grid.gs_G*u + grid.gs_c;
      end
    end

    function u = smoother_hybrid (grid, v, rhs, u)
      % u = grid.smoother_gauss_seidel (v, rhs, u);
      u = grid.smoother_chebyshev (v, rhs, u);
      % u = grid.smoother_chebyshev (v, rhs, u);
      u = grid.smoother_sym_sor(2, rhs, u);
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
      u0  = zeros(size(u));
      res = -rhs; % grid.residual(rhs, u);
      u1 = epsilon * res ;
      
      for iter = 1:v,                            % begin iteration
        res = grid.residual ( rhs, u );  
        u = alpha*u1 + (1-alpha)*u0 - epsalpha*res;
        u0 = u1; u1 = u;
        % r = norm(res);
        % n0 = norm(u0);
        % n1 = norm(u1);
        % disp([grid.dbg_spaces 'residual: ' num2str(r)]); % ' u0: ' num2str(n0) ' u1: ' num2str(n1)]);
      end % end iteration
      
    end % 2sr
    
    function u = smoother_chebyshev2 (grid, v, rhs, u)
      if ( isempty ( grid.eig_max ) )
        disp('computing eigenvalues');
        tic;
        % Kc = grid.Null' * grid.K * grid.Null;
        Kc = grid.K; %(eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
        % d = eigs(Kc, 2, 'be');
	grid.eig_max = eigs(Kc, 1, 'lm');  
        grid.eig_min = eigs(Kc, 1, 'sm');  
        toc;
      end
      
      % adjust the eigenvalues to hit the upper spectrum
      l_max = grid.eig_max;
      l_min = (grid.eig_min + grid.eig_max)/2;
      
      c = (l_min - l_max)/2;
      d = (l_min + l_max)/2;

      p = zeros(size(u));

      for iter = 1:v
        res = grid.residual ( rhs, u ); 
        %r = norm(res)
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

    function u = smoother_chebyshev (grid, v, rhs, u)
      if ( isempty ( grid.eig_max ) )
        % disp('computing eigenvalues');
        % Kc = grid.Null' * grid.K * grid.Null;
        D = diag(grid.K);
        grid.jacobi_invdiag = 1./D;
        Kc = spdiags(grid.jacobi_invdiag,0,length(D), length(D)) * grid.K;
        % d = eigs(Kc, 2, 'be');
        opts.tol = 0.01;
        grid.eig_max = eigs(Kc, 1, 'lm', opts);  
        % grid.eig_min = eigs(Kc, 1, 'sm');  
      end

      % adjust the eigenvalues to hit the upper spectrum
      beta = grid.eig_max;
      alpha = 0.25*grid.eig_max;% (grid.eig_min + grid.eig_max)/2;

      delta = (beta - alpha)/2;
      theta = (beta + alpha)/2;
      s1 = theta/delta;
      rhok = 1./s1;

      d = zeros(size(u));

      % first loop
      res = -grid.residual ( rhs, u );
      d = res/theta.* grid.jacobi_invdiag;
      u = u + d;

      for iter = 2:v
	  rhokp1 = 1/ (2*s1 - rhok);
	  d1 = rhokp1 * rhok;
	  d2 = 2*rhokp1 / delta;
	  rhok = rhokp1;
	  res = -grid.residual ( rhs, u ); 
	  %norm(res)
	  d = d1 * d + d2 * res.*grid.jacobi_invdiag;
	  u = u + d;
      end
    end % chebyshev


    function evec = get_eigenvectors(grid)
      % generate the correct matrix 
      Kc = grid.K; %(eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
      [evec, eval] = svd(full(Kc)); %eig(full(Kc), full(grid.M));
      [eval,per] = sort(diag(eval),'ascend');
      evec = evec(:,per);
      grid.k_evec = evec;
      grid.k_lam = eval;
    end
    
    function plot_spectrum(g, u, clr, rhs)
      if (g.debug)  
        subplot(1,2,1);
        a = g.M * u;
        q = repmat(a, 1, 80);
        b = abs(dot (g.k_evec, q));
        % plot eigenvalues 
        plot(b, clr); hold on;
        subplot(1,2,2); 
        rr = g.residual(rhs, u);
        n = sqrt(length(rr));
        imagesc(reshape(rr, n, n)); colorbar; hold off;
        grid on;
        odr = g.Mesh.fem.shape;
        set(gca, 'xtick', odr+0.5:odr:odr*g.Mesh.nelem);
        set(gca, 'ytick', odr+0.5:odr:odr*g.Mesh.nelem);
      end
    end

    function u0 = get_u0(grid)
      if (grid.debug)
        if ( isempty( grid.k_evec ) )
          [grid.k_evec, ~] = eigs(grid.K, grid.M, 80, 'BE');
        end
        n = size(grid.k_evec, 2);
        lam = ones(n,1);
        % lam(1:n/4) = 1;
        u0 = grid.k_evec*lam;
      else
        u0 = rand(size(grid.L()));
        u0(grid.Boundary) = 0;
      end
    end
  end %methods
  
end %classdef

