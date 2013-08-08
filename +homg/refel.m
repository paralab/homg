classdef refel < handle
    %REFEL hexahedral reference element
    %    - Equivalent to mangll_refel_t
    
    properties
        dim
        N      % polynomial order
        Nrp    % number of 1D interpolation points
        
        r      % 1D reference coordinates of the interpolation nodes ( 1 x Nrp )
        
        Vr     % 1D Vandermonde matrix of Legendre polynomials at r 
        gradVr %     and their derivative (Nrp x Nrp)
        
        Dr     % Derivative of Lagrange interpolants at the interpolation nodes.
               % ( Nrp x Nrp )
               % Dr(i,j) = lagrange_i' (r_j)
        
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
            
            elem.r      = homg.JacobiGL(0, 0, elem.N);
            
            elem.Vr     = zeros (order+1, order+1);
            elem.gradVr = zeros (order+1, order+1);
            for i=1:elem.Nrp
                elem.Vr(i,:)     = homg.JacobiP(elem.r, 0, 0, i-1);
                elem.gradVr(i,:) = homg.GradJacobiP(elem.r, 0, 0, i-1);
            end
        
            elem.Dr     = elem.Vr \ elem.gradVr;
            
            iVr         = elem.Vr \ eye(order+1);
            elem.Mr     = iVr * iVr';
            elem.invMr  = elem.Mr \ eye(order+1);
        end
    end
    
end

