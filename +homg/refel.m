classdef refel < handle
    %REFEL hexahedral reference element
    %    - Equivalent to mangll_refel_t
    
    properties
        dim
        N      % polynomial order
        Nrp    % number of 1D interpolation points
        
        r      % 1D reference coordinates of the interpolation nodes ( 1 x Nrp )
        g      % 1D reference coordinates of the gauss quadrature points
        w      % 1d weights for gauss quadrature
        
        W      % dim-dimensional weights 
        
        Vr     % 1D Vandermonde matrix of Legendre polynomials at r 
        gradVr %     and their derivative (Nrp x Nrp)
        
        Vg     % 1D Vandermonde matrix of Legendre polynomials at g 
        gradVg %     and their derivative (Nrp x Nrp)
        
        Dr     % Derivative of Lagrange interpolants at the interpolation nodes.
               % ( Nrp x Nrp )
               % Dr(i,j) = lagrange_i' (r_j)
        Dg     % Derivative of Lagrange interpolants at the gauss quad points
               % ( Nrp x Nrp )
               % Dr(i,j) = lagrange_i' (r_j)
        
        Q      % map to gauss points 
    
        Qx
        Qy
        Qz
        
        Mr     % exact 1D Mass matrix (Nrp x Nrp)
        invMr  % and its inverse
    end
    
    methods
        function elem = refel(d, order)
            % Setup a d-dimensional reference element 
            % order = polynomial order of the reference element
            elem.dim    = d;
            elem.N      = order;
            elem.Nrp    = order + 1;
            
            elem.r      = homg.jacobi.gll (0, 0, elem.N);
            
            [elem.g, elem.w] = homg.jacobi.gauss(0, 0, elem.N);
            
            elem.Vr     = zeros (order+1, order+1);
            elem.gradVr = zeros (order+1, order+1);
            
            elem.Vg     = zeros (order+1, order+1);
            elem.gradVg = zeros (order+1, order+1);
            
            for i=1:elem.Nrp
                elem.Vr(i,:)     = homg.jacobi.polynomial (elem.r, 0, 0, i-1);
                elem.gradVr(i,:) = homg.jacobi.gradient (elem.r, 0, 0, i-1);
                
                elem.Vg(i,:)     = homg.jacobi.polynomial (elem.g, 0, 0, i-1);
                elem.gradVg(i,:) = homg.jacobi.gradient (elem.g, 0, 0, i-1);
            end
        
            elem.Dr     = transpose(elem.Vr \ elem.gradVr);
            
            elem.Dg     = transpose(elem.Vr \ elem.gradVg);
            
            iVr         = elem.Vr \ eye(order+1);
            
            q1d         = transpose (elem.Vr \ elem.Vg);  
            
            elem.W           = zeros(elem.Nrp^elem.dim, 1);
            
            if (d == 2)
              elem.Q  = kron(q1d, q1d) ;
              
              elem.Qx = kron(q1d, elem.Dg);
              elem.Qy = kron(elem.Dg, q1d);
              
              sk = 1;
              for i=1:elem.Nrp
                for j=1:elem.Nrp
                  elem.W(sk) = elem.w(i) * elem.w(j);
                  sk = sk + 1;
                end
              end
              
            else
              elem.Q  = kron(kron(q1d, q1d), q1d);
              
              elem.Qx = kron(kron(q1d, q1d), elem.Dg);
              elem.Qy = kron(kron(q1d, elem.Dg), q1d);
              elem.Qz = kron(kron(elem.Dg, q1d), q1d);
              
              sk = 1;
              for i=1:elem.Nrp
                for j=1:elem.Nrp
                  for k=1:elem.Nrp
                    elem.W(sk) = elem.w(i) * elem.w(j) * elem.w(k);
                    sk = sk + 1;
                  end
                end
              end
            end
            
            elem.Mr     = iVr * iVr';
            elem.invMr  = elem.Mr \ eye(order+1);
            
        end
    end
    
end

