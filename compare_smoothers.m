function compare_smoothers (dim, geom, order, nelem) 
% function compare_smoothers (dim, geom, order, nelem) 
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

clrs = 'bkrgcym';


% smooth
if (dim == 2)
  g.jacobi_omega = 2.0/3.0;
else
  g.jacobi_omega = 6.0/7.0;
end


r0 = evec*lam;
u0 = g.K \ (r0);

close all;
% figure('DefaultAxesFontSize', 20);
% figure (1);
hFig = figure(1);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [200 200 800 800])


% clf; % hold on;
% u0 = evec*lam;
r = g.residual(g.L, u0);
% q = repmat(u0,size(u0'));
b = r' * evec; % dot (evec, q);
semilogy(b, 'k'); hold on;

% jacobi 
g.set_smoother('jacobi');
u1 = u0; %evec*lam;
u = g.smooth(3, g.L, u1);
% compute projections ...
% q = repmat(u,size(u'));
r = g.residual(g.L, u);
b = r' * evec; % abs(dot (evec, q));
% plot eigenvalues
semilogy(abs(b), 'b.', 'MarkerSize', 5); hold on;

% chebyshev 
g.set_smoother('chebyshev');
u1 = u0; % evec*lam;
u = g.smooth(3, g.L, u1);
% compute projections ...
% q = repmat(u,size(u'));
r = g.residual(g.L, u);
b = r' * evec; % abs(dot (evec, q));
% plot eigenvalues
semilogy(abs(b), 'm', 'LineWidth', 3); 

% ssor 
g.set_smoother('ssor');
u1 = u0; % evec*lam;
u = g.smooth(2, g.L, u1);
% compute projections ...
% q = repmat(u,size(u'));
r = g.residual(g.L, u);
b = r' * evec; %abs(dot (evec, q));
% plot eigenvalues
semilogy(abs(b), 'g.', 'MarkerSize', 5); 

title(['N = ' num2str(order*nelem+1) '^2' ' , p = ' num2str(order)]);

grid on;

ylim([0.00001, 2]);

print ('-depsc2', ['smoothers-order' num2str(order) '.eps']);
matlab2tikz (['smoothers-order' num2str(order) '.tikz']);
