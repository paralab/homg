classdef freesurface
    %FREESURFACE: The freesurface class contains the code that intializes
    %and solves the free surface sloshing problem with surface tension.
    %   initialize_problem:
    %       This method uses the homg code to construct the stiffness and
    %       mass matrices for the 3D and 2D portions of the free surface
    %       problem. Once these matrices have been constructed,
    %       assemble_freesurface_mixed is called to construct the matrix
    %       corresponding to the coupling terms.
    %   solve:
    %       With the problem initialized for a given container type and 
    %       problem size, this method solves the system for a number of
    %       eigenvalues and eigenvectors that will be supplied by the user.
    %       The first (1+kd)^3 elements of the eigenvectors correspond to
    %       the velocity potential and the final (1+kd)^2 elements
    %       correspond to the free surface displacement.
    %   assemble_freesurface_mixed:
    %       This method then constructs the matrix that couples the 3D
    %       portion of the problem with the 2D portion.
    methods (Static)
        function [ S, SF, MF, Mm ] = initialize_problem( d, k, container_type)
            if strcmpi(container_type, 'box')
                ctype3 = @homg.xform.identity;
                ctype2 = @homg.xform.identity;
            elseif strcmpi(container_type, 'bowl')
                ctype3 = @homg.xform.bowl;
                ctype2 = @homg.xform.disc3;
            elseif strcmpi(container_type, 'cylinder')
                ctype3 = @homg.xform.cylinder;
                ctype2 = @homg.xform.disc3;
            elseif strcmpi(container_type, 'pcylhalf')
                ctype3 = @homg.xform.pcylhalf;
                ctype2 = @homg.xform.disc3;
            elseif strcmpi(container_type, 'pcyldouble')
                ctype3 = @homg.xform.pcyldouble;
                ctype2 = @homg.xform.disc3;
            else
                ctype3 = @homg.xform.identity;
                ctype2 = @homg.xform.identity;
            end
            m1 = homg.hexmesh([d d d],ctype3);
            g1 = homg.grid(m1,k);
            S = g1.Mesh.assemble_stiffness(k);

            m2 = homg.hexmesh([d d],ctype2);
            g2 = homg.grid(m2,k);
            SF = g2.Mesh.assemble_stiffness(k);
            MF = g2.Mesh.assemble_mass(k);
            Mm = freesurface.assemble_freesurface_mixed(MF,d,k,3);
        end
        function [ Mm ] = assemble_freesurface_mixed( MF, d, k, dim )
            % MM and MMt matrices have the same values as MF but are each shifted to
            % their respective positions
            [I_MM, J_MM, V_MM] = find(MF);
            I_MM = I_MM + (1+k*d)^dim - (1+k*d)^(dim-1);
            Mm = sparse(I_MM,J_MM,V_MM);
        end
        function [ evals, evecs ] = solve(Bo, sigma, n_evals, S, MF, SF, Mm)
            n1 = size(S,1); n2 = size(MF,1);
            [i,j,s] = find(S);
            A = sparse(i,j,s,n1+n2,n1+n2);
            [i,j,s] = find(MF + (1/Bo)*SF);
            A = A + sparse(i+n1,j+n1,s);

            B = sparse([],[],[],n1+n2,n1+n2);
            [i,j,s] = find(Mm);
            B = B + sparse(i,j+n1,s,n1+n2,n1+n2) + sparse(j+n1,i,s,n1+n2,n1+n2);

            [evecs, evals] = eigs(A,B,n_evals,sigma);
            evals = sqrt(diag(evals));

            I = find (real(evals) > 0.1);
            evecs = evecs(:,I);
            evals = real(evals(I));

            s = size(evecs);
            for i = 1:s(2)
               alpha = transpose(evecs(:,i))*B*evecs(:,i);
               alpha = 1.0 / sqrt(alpha);
               evecs(:,i) = alpha * evecs(:,i);
               if evecs(1,i) < 0
                   evecs(:,i) = (-1.0) * evecs(:,i);
               end
            end

        end
    end    
end

