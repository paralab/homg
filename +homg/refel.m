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
        
        % Prolongation 
        P      % interpolation from this element to its 4/8 children
        
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
            
            elem.r      = homg.basis.gll (0, 0, elem.N);
            
            rp          = [0.5*(elem.r - 1); 0.5*(elem.r(2:end) + 1)];
            
            [elem.g, elem.w] = homg.basis.gauss(0, 0, elem.N);
            
            elem.Vr     = zeros (order+1, order+1);
            elem.gradVr = zeros (order+1, order+1);
            
            elem.Vg     = zeros (order+1, order+1);
            elem.gradVg = zeros (order+1, order+1);
            
            Vp     = zeros (order+1, 2*order+1);
            
            for i=1:elem.Nrp
                elem.Vr(i,:)     = homg.basis.polynomial (elem.r, 0, 0, i-1);
                elem.gradVr(i,:) = homg.basis.gradient (elem.r, 0, 0, i-1);
                
                elem.Vg(i,:)     = homg.basis.polynomial (elem.g, 0, 0, i-1);
                elem.gradVg(i,:) = homg.basis.gradient (elem.g, 0, 0, i-1);
                
                Vp(i,:)          = homg.basis.polynomial (rp, 0, 0, i-1);
            end
        
            elem.Dr     = transpose(elem.Vr \ elem.gradVr);
            
            elem.Dg     = transpose(elem.Vr \ elem.gradVg);
            
            iVr         = elem.Vr \ eye(order+1);
            
            q1d         = transpose (elem.Vr \ elem.Vg);  
            p1d         = transpose (elem.Vr \ Vp);  
            
            elem.W      = zeros(elem.Nrp^elem.dim, 1);
            
            if (d == 2)
              elem.Q  = kron(q1d, q1d) ;
              elem.P  = kron(p1d, p1d) ;
              
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
              elem.P  = kron(kron(p1d, p1d), p1d);
              
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

