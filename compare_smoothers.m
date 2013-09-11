function [res, evec] = compare_smoothers (dim, xform, mu, order, nelem, sfix) 
% function compare_smoothers (dim, geom, order, nelem) 
% 
% valid smoothers are,
%       'jacobi'
%       'chebyshev'
%       'gs'
%       '2sr'

% generate matrix 
% m = homg.mesh(dim, geom, nelem, 1);
% m.set_rhs('0');
% m.set_rhs('-8*pi^2*(sin(2*pi*x) * sin(2*pi*y))');
% g = homg.grid(m, order);

% g = create_grid_hierarchy(dim, geom, order, [nelem/4 nelem/2 nelem], 1);
g = create_hexmesh_grids(dim, mu, xform, order, [nelem/4 nelem/2 nelem]);

fname1 = ['smoother-' sfix '.dat'];
fname2 = ['vcycle-' sfix '.dat'];

nun = nelem*order + 1;

% get eigenvectors 
evec = g.get_eigenvectors();
n = size(evec, 1);

% generate u
lam = ones(n,1);

clrs = 'bkrgcym';

res = zeros(nun, nun, 4);

% smooth
if (dim == 2)
  g.jacobi_omega = 2.0/3.0;
else
  g.jacobi_omega = 6.0/7.0;
end


u0 = evec*lam;
% u0 = g.K \ (r0 + g.L);

close all;
% figure('DefaultAxesFontSize', 20);
% figure (1);
hFig = figure(1);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [200 200 800 800])


% clf; % hold on;
% u0 = evec*lam;
% r = g.residual(g.L, u0);
% q = repmat(u0,size(u0'));
b = u0' * evec; % dot (evec, q);
semilogy(abs(b), 'k'); hold on;

%~~~~~~~~~~~~~~~~~~~~~~~~~ 
if (order == 1)
  dlmwrite(fname1,abs(b),'delimiter','\t','precision',6);
  dlmwrite(fname2,abs(b),'delimiter','\t','precision',6);
end


res(:,:,1) = reshape(u0, nun, nun);

% jacobi 
g.set_smoother('jacobi');
u1 = u0; %evec*lam;
u = g.smooth(4, g.L, u1);
% compute projections ...
% q = repmat(u,size(u'));
% r = g.residual(g.L, u);
b = u' * evec; % abs(dot (evec, q));
% plot eigenvalues
semilogy(abs(b), 'b'); %, 'MarkerSize', 5); hold on;

dlmwrite(fname1,abs(b),'delimiter','\t','precision',6, '-append');

res(:,:,2) = reshape(u, nun, nun);

% chebyshev 
g.set_smoother('chebyshev');
u1 = u0; % evec*lam;
u = g.smooth(4, g.L, u1);
%g.set_smoother('ssor');
%%u = g.smooth(2, g.L, u);
% compute projections ...
% q = repmat(u,size(u'));
r = g.residual(g.L, u);
b = u' * evec; % abs(dot (evec, q));
% plot eigenvalues
semilogy(abs(b), 'm'); 
dlmwrite(fname1,abs(b),'delimiter','\t','precision',6, '-append');

res(:,:,3) = reshape(u, nun, nun);

% ssor 
g.set_smoother('ssor');
u1 = u0; % evec*lam;
u = g.smooth(2, g.L, u1);
% compute projections ...
% q = repmat(u,size(u'));
r = g.residual(g.L, u);
b = u' * evec; %abs(dot (evec, q));
% plot eigenvalues
semilogy(abs(b), 'g'); %, 'MarkerSize', 5); 
dlmwrite(fname1,abs(b),'delimiter','\t','precision',6, '-append');

res(:,:,4) = reshape(u, nun, nun);

title(['N = ' num2str(order*nelem+1) '^2' ' , p = ' num2str(order)]);

grid on;

% ylim([0.001, 10]);

% print ('-depsc2', ['smoothers-order' num2str(order) sfix '.eps']);
% matlab2tikz (['smoothers-order' num2str(order) sfix '.tikz'], 'checkForUpdates', false, 'showInfo', false);

% Now for one v-cycle ...
% close all;
% figure('DefaultAxesFontSize', 20);
% figure (1);
hFig = figure(2);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [200 200 800 800])


r = g.residual(g.L, u0);
b = u0' * evec; % dot (evec, q);
semilogy(abs(b), 'k'); hold on;

% jacobi 
u1 = u0; %evec*lam;
[u, rr, iter3] = g.solve(1, 'jacobi', 4, g.L, u1); 
r = g.residual(g.L, u);
b = u' * evec; % abs(dot (evec, q));
% plot eigenvalues
semilogy(abs(b), 'b'); %, 'LineWidth', 3); hold on;
dlmwrite(fname2,abs(b),'delimiter','\t','precision',6, '-append');

% chebyshev 
u1 = u0; % evec*lam;
[u, rr, iter3] = g.solve(1, 'chebyshev', 4, g.L, u1); 
%[u, rr, iter3] = g.solve(1, 'ssor', 2, g.L, u); 
% compute projections ...
% q = repmat(u,size(u'));
r = g.residual(g.L, u);
b = u' * evec; % abs(dot (evec, q));
% plot eigenvalues
semilogy(abs(b), 'm'); 
dlmwrite(fname2,abs(b),'delimiter','\t','precision',6, '-append');
% ssor 
u1 = u0; % evec*lam;
[u, rr, iter3] = g.solve(1, 'ssor', 2, g.L, u1); 
% compute projections ...
r = g.residual(g.L, u);
b = u' * evec; %abs(dot (evec, q));
% plot eigenvalues
semilogy(abs(b), 'g'); %, 'LineWidth', 3); 
dlmwrite(fname2,abs(b),'delimiter','\t','precision',6, '-append');

title(['N = ' num2str(order*nelem+1) '^2' ' , p = ' num2str(order)]);

grid on;

% ylim([0.00001, 1.2]);

% print ('-depsc2', ['vcycle-order' num2str(order) sfix '.eps']);
% matlab2tikz (['vcycle-order' num2str(order) sfix '.tikz'], 'checkForUpdates', false, 'showInfo', false);

% close all;
% hFig = figure(1);
% set(gcf,'PaperPositionMode','auto')
% set(hFig, 'Position', [200 200 800 800])
% 
% subplot(2,2,1);
% surf(res(:,:,1)); axis([0,nun,0,nun,-10,10]); %axis equal;
% subplot(2,2,2);
% surf(res(:,:,2)); axis([0,nun,0,nun,-10,10]); %axis equal;
% subplot(2,2,3);
% surf(res(:,:,3)); axis([0,nun,0,nun,-10,10]); %axis equal;
% subplot(2,2,4);
% surf(res(:,:,4)); axis([0,nun,0,nun,-10,10]); %axis equal;
