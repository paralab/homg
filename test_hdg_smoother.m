function test_hdg_smoother(A, rhs)

% generate random vector ...
u = rand(size(rhs)); 

% diagonal of operator
D    = diag(A);
invD = 1./D;


plot(u); hold on;

% smooth using Jacobi & Chebyshev
us = smoother_jacobi (A, invD, rhs, u, 100);

plot(us, 'r');

end


function r = hdg_residual(A, rhs, u)
  r = A*u - rhs;
end 

function u = smoother_jacobi (K, invD, rhs, u, v)
  omega = 2/3;
  % standard jacobi smoother
  for i=1:v
    res  = invD .* hdg_residual(K, rhs, u);
    
    u = u - omega .* res;
    r = norm(res);
    disp([num2str(i) ' : ' num2str(r)]);
  end
end % jacobi