classdef tensor
    %TENSOR Basic tensor routines
    
    methods(Static)
        % 2D routines
        function y = IAX (A, x)
            N = size (A, 1);
            y = A * reshape(x, N, N);
            y = y(:);
        end
        
        function y = AIX (A, x)
            N = size (A, 1);
            y = A * reshape(x, N, N)';
            y = y'; 
            y = y(:);
        end
        
        % 3D routines
        function y = IIAX (A, x)
            N = size (A, 1);
            y = reshape(x, N*N, N) * A';
            y = y(:);
        end
        
        function y = IAIX (A, x)
            N = size (A, 1);
            q = reshape(x, N, N, N);
            y = zeros(N,N,N);
            for i=1:N
                y(i,:,:) = A * squeeze( q(i,:,:) );
            end
            y = y(:);
        end
        
        function y = AIIX (A, x)
            N = size (A, 1);
            y = A * reshape(x, N, N*N);
            y = y(:);
        end
        
    end
    
end

