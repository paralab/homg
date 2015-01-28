function[u_l2err, qx_l2err, qy_l2err] = HDG(nElem,nPoly)
%clc;
%clear all
% Tan Bui and Hari Sundar, Oct 28, 2013
% Testing HDG method for 2D Laplace equation

% addpath /workspace/tanbui/tanbui/WithHari/homg/

% HDG structure
HDGdata = [];

% solution order
order = nPoly;

% number of elements in x and y directions
nelems = [nElem,nElem];

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
qxexact = @(pts) -0.4/pi * (cos(2*pi * pts(:,1)) ...  %since q=-grad u there is negative sign
                           .* sin(pi * pts(:,2)));
qyexact = @(pts) -0.2/pi * (sin(2*pi * pts(:,1)) ...
                           .* cos(pi * pts(:,2)));

% Number of volume unknown for a scalar
rhsqx = zeros(Nv,1);
rhsqy = zeros(Nv,1);
rhsu  = zeros(Nv,1);

% predefined normal vector, don't like it but stick with it for now
nx = [-1, 1, 0, 0];
ny = [0, 0, -1, 1];
HDGdata.nx = nx;
HDGdata.ny = ny;

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
    if ((e1 > 0) && (e2 > 0)) 
      % global trace index for f1h = [4,8,16,32];

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

%-------- form the RHS------------------
lamInterior = zeros(size(HDGdata.SkelInterior2All));
%rhs = -residual(lamInterior,HDGdata,forcing, Bdata);
rhs = -residualFast(lamInterior,HDGdata,forcing, Bdata);
%--------- end form the RHS------------

% form the HDG matrix
A = HDGmatrix(HDGdata);


% $$$ % Quick, dirty, and expensive way
%--------- Construct the HDG matrix-------
% forcingn = @(pts) zeros(size(pts,1),1);
% Bdatan = zeros(size(Bdata));
% Nh = Nfp*Nifaces; 
% maxnnzeros = 5 * Nh;
% II = zeros(maxnnzeros,1);
% JJ = zeros(maxnnzeros,1);
% SS = zeros(maxnnzeros,1);
% $$$ 
% nnzeros = 0;
% $$$ 
% for n = 1:Nh
%   lamInterior(n) = 1;
% $$$   
%     Aj = residual(lamInterior,HDGdata,forcingn, Bdatan);    
     %Aj = residualFast(lamInterior,HDGdata,forcingn, Bdatan);
% $$$   
%   [i,j,s] = find(Aj); j(:) = n;
% $$$   
%ni = length(i);
% $$$   
%   if (nnzeros + ni) > maxnnzeros,
%     maxnnzeros = 2*maxnnzeros;
%     II(maxnnzeros) = 0;
%     JJ(maxnnzeros) = 0;
%     SS(maxnnzeros) = 0;
%   end
% $$$   
%   maxnnzeros = 2*maxnnzeros;
% $$$   
%   II(nnzeros+1:nnzeros+ni) = i;
%   JJ(nnzeros+1:nnzeros+ni) = j;
%   SS(nnzeros+1:nnzeros+ni) = s;
% $$$   
%   nnzeros = nnzeros + ni; 
% $$$   
%   lamInterior(n) = 0;
% $$$   
% end
% Aa = sparse(II(1:nnzeros),JJ(1:nnzeros),SS(1:nnzeros),Nh,Nh);
% $$$ 
% $$$ keyboard

% Now solve for lam
lamInterior = A \ rhs;

eig_A=eig((A+A')/2);

eig_A


lamAll = zeros(Nsfaces * Nfp,1);
lamAll(Bmaps) = Bdata;
lamAll(SkelInterior2All) = lamInterior;

u = zeros(Nv,K);
qx = zeros(Nv,K);
qy = zeros(Nv,K);

L2eu = 0;
L2eqx = 0;
L2eqy = 0;

for e = 1:K
  [u(:,e),qx(:,e),qy(:,e)] = localSolver(HDGdata, e, lamAll, ...
                                         forcing);
  
  eu = u(:,e) - Uexact(:,e);
  eqx = qx(:,e) - Qxexact(:,e);
  eqy = qy(:,e) - Qyexact(:,e);
  
  pts = m.element_nodes(e, refel);
  [Jv, Dv] = m.geometric_factors(refel, pts);
  eMat = m.element_mass(e, refel, Jv);
  
  L2eu = L2eu +  eu' * eMat * eu;
  L2eqx = L2eqx +  eqx' * eMat * eqx;
  L2eqy = L2eqy +  eqy' * eMat * eqy;
  
end
u_l2err=0;
qx_l2err=0;
qy_l2err=0;
%fprintf('L2 norm error for u  = %1.15e \n',sqrt(L2eu));
%fprintf('L2 norm error for qx = %1.15e \n',sqrt(L2eqx));
%fprintf('L2 norm error for qy = %1.15e \n',sqrt(L2eqy));
u_l2err=sqrt(L2eu);
qx_l2err=sqrt(L2eqx);
qy_l2err=sqrt(L2eqy);
fprintf('L2 norm error for u  = %1.15e \n',u_l2err);
fprintf('L2 norm error for qx = %1.15e \n',qx_l2err);
fprintf('L2 norm error for qy = %1.15e \n',qy_l2err);


end



