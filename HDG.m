clear all
% Tan Bui and Hari Sundar, Oct 28, 2013
% Testing HDG method for 2D Laplace equation

% addpath /workspace/tanbui/tanbui/WithHari/homg/

% HDG structure
HDGdata = [];

% solution order
order = 6;

% number of elements in x and y directions
nelems = [6,3];

% generate the hexmesh with identity transform for now
m = homg.hexmesh(nelems,@homg.xform.identity);
m.set_order(order);
HDGdata.m = m;

% reference elements
refel = homg.refel(m.dim, order);
HDGdata.refel = refel;

% get total number of faces on the skeleton of the mesh
Nsfaces =m.get_num_faces();
HDGdata.Nsfaces = Nsfaces;

% the number of face points
Nfp = refel.Nrp ^ (refel.dim-1);
HDGdata.Nfp = Nfp;

% the number of volume point
Nv = refel.Nrp ^ (refel.dim);
HDGdata.Nv = Nv;

% number of faces
Nfaces = refel.dim * 2;
HDGdata.Nfaces = Nfaces;

% number of elements
K = prod(m.nelems);
HDGdata.K = K;

% Initialize lam
lam = zeros(Nfp * Nsfaces,1);

% lambda residual
lamRes = zeros(Nfp * Nsfaces, 1);

% forcing
forcing = @(pts) (sin(2.0 * pi * pts(:,1)) .* sin(pi * pts(:,2)));

% exact solution
% u = 0.5/pi^2 * forcing;
uexact = @(pts) 0.2 / pi^2 * forcing(pts);
qxexact = @(pts) 0.4/pi * (cos(2*pi * pts(:,1)) ...
                           .* sin(pi * pts(:,2)));
qyexact = @(pts) 0.2/pi * (sin(2*pi * pts(:,1)) ...
                           .* cos(pi * pts(:,2)));

% Number of volume unknown for a scalar
rhsqx = zeros(Nv,1);
rhsqy = zeros(Nv,1);
rhsu  = zeros(Nv,1);

% predefined normal vector, don't like it but stick with it for now
nx = [-1, 1, 0, 0];
ny = [0, 0, -1, 1];

% stabilization parameter
taur = 1;
HDGdata.taur = taur;
uu = zeros(Nv,Nv);
uqx = zeros(Nv,Nv);
uqy = zeros(Nv,Nv);

% Also compute uexact for testing
Uexact = zeros(Nv,K);
Qxexact = zeros(Nv,K);
Qyexact = zeros(Nv,K);
for k = 1:K
  pts = m.element_nodes(k, refel);
  Uexact(:,k)  = uexact(pts);
  Qxexact(:,k) = qxexact(pts);
  Qyexact(:,k) = qyexact(pts);
end

% COMPUTE EXACT LAMBDA FOR TESTING
for  sf=1:Nsfaces
    [e1, f1, e2, f2]  =m.get_face_elements(sf);

    %    if (e1 > 0) && (e2 > 0), % interior faces
    if (f1 > 0)
      pts = m.element_nodes(e1, refel);
      idxf = m.get_skeletal_face_indices(refel, e1, f1);      
      idxv = m.get_discontinuous_face_indices(refel, 1, f1);
      lam(idxf) = uexact(pts(idxv,:));
    end
    if (f2 > 0)
      pts = m.element_nodes(e2, refel);
      idxf = m.get_skeletal_face_indices(refel, e2, f2);      
      idxv = m.get_discontinuous_face_indices(refel, 1, f2);
      lam(idxf) = uexact(pts(idxv,:));
    end
end

% Construct the Lift and VtoF
LIFT = zeros(Nv, Nfp, Nfaces);
VtoF = zeros(Nfp, Nv, Nfaces);
for f = 1:Nfaces
  idxv = m.get_discontinuous_face_indices(refel, 1, f);
  LIFT(idxv,:,f) = refel.Mr;
  for fp = 1:Nfp
    VtoF(fp,idxv(fp),f) = 1;
  end
end
HDGdata.LIFT = LIFT;
HDGdata.VtoF = VtoF;

% find boundary faces and indices
Nbfaces = 0; Bmaps = zeros(Nfp,1); Bdata = zeros(Nfp,1);
iindex = 0;
for  sf=1:Nsfaces
    [e1, f1, e2, f2]  =m.get_face_elements(sf);

    if (e1 < 0) || (e2 < 0), % boundary faces
      if e1 < 0, 
        e = e2; f = f2;
      else
        e = e1; f = f1;
      end
      Nbfaces = Nbfaces + 1;
      idxf = m.get_skeletal_face_indices(refel, e, f);
      idxv = m.get_discontinuous_face_indices(refel, 1, f);
      
      Bmaps(iindex+1:iindex+Nfp) = idxf;
      Bdata(iindex+1:iindex+Nfp) = Uexact(idxv,e);
      iindex = iindex + Nfp;
    end
end
HDGdata.Nbfaces = Nbfaces;
HDGdata.Bmaps = Bmaps;

Nifaces = Nsfaces - Nbfaces;
SkelInterior2All = zeros(Nifaces * Nfp,1);
SkelAll2Interior = zeros(Nsfaces * Nfp,1);
InteriorF2AllF   = zeros(Nifaces,1);

HDGdata.Nifaces = Nifaces;

% Construct the skeleton maps
iindex = 0; iface = 0;
for  sf=1:Nsfaces
    [e1, f1, e2, f2]  =m.get_face_elements(sf);
    
    % interior faces
    if (e1 > 0) && (e2 > 0), 
      % global trace index for f1
      idxf = m.get_skeletal_face_indices(refel, e1, f1);
      SkelInterior2All(iindex+1:iindex+Nfp) = idxf;
      SkelAll2Interior(idxf) = iindex+1:iindex+Nfp;
      iindex = iindex + Nfp;
      iface = iface + 1;
      InteriorF2AllF(iface) = sf;
    end
end
HDGdata.SkelInterior2All = SkelInterior2All;
HDGdata.SkelAll2Interior = SkelAll2Interior;
HDGdata.InteriorF2AllF = InteriorF2AllF;

% Form the HDG matrix and RHS
% Quick, dirty, and expensive way

% loop over all faces of the mesh skeleton
for  sf=1:Nsfaces
    [e1, f1, e2, f2]  =m.get_face_elements(sf);

    %% Task 1
    if (e1 > 0) && (e2 > 0), % interior faces
      
      % e1 solution
      [u1,qx1,qy1] = localSolver(m, refel, e1, lam, taur,...
                                 forcing, LIFT, VtoF);

      % e2 solution
      [u2,qx2,qy2] = localSolver(m, refel, e2, lam, taur,...
                                 forcing, LIFT, VtoF);

      keyboard
    else
      continue;
    end
      
end
