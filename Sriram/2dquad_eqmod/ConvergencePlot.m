clear all
% Tan Bui, July 12, 2014
% For Shinhoo and Sriram
% Convergence plot
% This file plots h-convergence
% Cell A{i}: is for the ith order solution
% For each cell A{i}: First column is solution order (p), second
% column is mesh size (h), third column is the L2 error
% corresponding to each mesh size. Top row is for the coarsest mesh
% and last row is for the finest mesh (you can have as many rows as the
% number of meshes).
% You can have as many A{i} as the number of solution orders you test

%------------------------------------------------------------------
% h: # of elements in 1 direction in the  domain
% p: the degree of polynomials


iVar = 1; % 1 for u, 2 for qx, 3 for qy

%h = [4,8,16,32];
h=[2];
%p = [1,2,3,4];
p=[2];
nh = length(h);
np = length(p);
nVar = 3; %u, qx, qy
l2err = zeros(nh,np,nVar);
u_l2err=0;
qx_l2err=0;
qy_l2err=0;

fid=fopen('L2_Error_HDG_tau_eqmodfinal.txt','w');
fprintf(fid,'K \t N \t U_L2_Error \t qx_L2_Error \t qy_L2_Error\r\n');


%---------------------
% calculate l2 errors
%---------------------
for ip=1:np
for ih=1:nh

    nElem = h(ih);
    nPoly = p(ip);
    
    [u_l2err, qx_l2err, qy_l2err] = HDG(nElem,nPoly);
      
   
    l2err(ih,ip,1) = u_l2err;
    l2err(ih,ip,2) = qx_l2err;
    l2err(ih,ip,3) = qy_l2err;
    
  fprintf(fid,'%g \t %g \t %g \t %g \t %g\r\n',nElem,nPoly,u_l2err,qx_l2err,qy_l2err);

end
fprintf(fid,'\r\n');
end

fclose(fid);




%---------------------
% Build A
%---------------------
A = cell(np,1);

for ip=1:np
    A{ip} = zeros(nh,3);
    for ih=1:nh
        A{ip}(ih,1) = p(ip);
        A{ip}(ih,2) = 1/h(ih);
        A{ip}(ih,3) = l2err(ih,ip,iVar);
    end
end




%% For your particular problem, replace h and the L2 error norm correspondingly.

% % First order solution
% A{1} = [1 0.015625000000000 0.009092198566947
% 1 0.007812500000000   0.001941095878590
% 1 0.003906250000000   0.000472314346570
% 1 0.001953125000000   0.000115962509411];
% 
% % Second order solution
% A{2} = [2 0.015625000000000 0.001529334669987
% 2 0.007812500000000   0.000358791388418
% 2 0.003906250000000   0.000042735296594
% 2 0.001953125000000   0.000005086547954];
% 
% % Third order solution 
% A{3} = [3 0.015625000000000 6.590037046239085e-04
% 3 0.007812500000000   0.000049080023736
% 3 0.003906250000000   0.000002413575104
% 3 0.001953125000000   0.000000150814174];
% 
% % Fourth order solution
% A{4} = [4 0.015625000000000 1.326513627798934e-04
% 4 0.007812500000000   0.000002402383781
% 4 0.003906250000000   0.000000169556135
% 4 0.001953125000000   0.000000005378503];
% 
% % Fifth order solution
% A{5} = [5 0.015625000000000 1.465040597678935e-05
% 5 0.007812500000000   0.000000956316574
% 5 0.003906250000000   0.000000013931646
% 5 0.001953125000000   0.000000000196333];
% 
% % You can add 6th order and so on.....

%%%%%%%%%%%%%%%------------h convergence----------------------------
N = length(A);
order = zeros(N,1);
figure
axes('fontsize',14);
hold on;

for n=1:N
  h = A{n}(:,2); err = sqrt(sum(A{n}(:,3:end).^2,2));
  % slope
  temp = polyfit(log(h),log(err),1);
  order(n) = temp(1);
  %order(n) = (log(err(end))-log(err(end-1)))/(log(h(end))-log(h(end-1)));
  loglog(h,err,'linewidth',1.5)
  
  % Draw the solution order
  text(h(1)+5.e-4,err(1),strcat('p = ',num2str(A{n}(1,1))))
  
  % drawing the slope, working between log and linear scale
  x1 = h(end); x2 = h(end-1)-0.5*(h(end-1)-h(end));
  y1 = err(end) + err(end)*0.2;
  b = log(y1) - order(n)*log(x1);
  y2 = exp(order(n)*log(x2) + b);
  x3 = x1;
  y3 = y2;
  px = [x1,x2,x3,x1];
  py = [y1,y2,y3,y1];
  loglog(px,py,'r','linewidth',1.5)
  
  % plot the order of convergene
  xm = x1; xm = xm - 0.2*xm;
  ym = y1 + 0.3*(y3-y1); 
  text(xm,ym,num2str(order(n),3))
  
  xm = x1 + 0.3*(x2-x1);
  ym = y2 + 0.3*y2; 
  text(xm,ym,'1')
  
end
%xlim([0.04 0.5])
%axis([0.02 0.45 1.e-10 1])
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('h'); ylabel('error in L^2-norm')
