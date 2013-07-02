function test_smoother (dim, geom, order, nelem, smoother) 
% function test_smoother (dim, order, nelem, smoother) 
% 
% valid smoothers are,
%       'jacobi'
%       'chebyshev'
%       'gs'
%       '2sr'

% generate matrix 
m = homg.mesh(dim, geom, nelem, 1);
m.set_rhs('0');
% m.set_rhs('-8*pi^2*(sin(2*pi*x) * sin(2*pi*y))');
g = homg.grid(m, order);

% get eigenvectors 
evec = g.get_eigenvectors();
n = size(evec, 1);

% generate u
lam = ones(n,1);


% smooth
g.set_smoother(smoother);
g.jacobi_omega = 2.0/3.0;
clrs = 'bkrgcym';
stps = [1 2 3 4 5 6 8 12 16 20];
% stps = [order-1 order order+3];

clf;
u0 = evec*lam;
q = repmat(u0,size(u0'));
b = dot (evec, q);
plot(b, 'k'); hold on;

for i=1
  u0 = evec*lam; % compare with sum(evec, 2);
  u = g.smooth(stps(i), g.L, u0);
  % compute projections ...
  q = repmat(u,size(u'));
  b = abs(dot (evec, q));
  % plot eigenvalues 
  plot(b, clrs(i)); % hold on;
end

grid on;
