function test_gradient( order )
% function test_gradient(nelems, transform, order)

%% 2D Case
quad = homg.refel(2, order);


f2 = @(x,y)(x*y^order + y - 3);
g2x = @(x,y)(y^(order));  
g2y = @(x,y)(x*order*y^(order-1) + 1);

[X, Y] = ndgrid(quad.r, quad.r);

u    = arrayfun(f2,   X(:), Y(:)); 
dudx = arrayfun(g2x,  X(:), Y(:)); 
dudy = arrayfun(g2y,  X(:), Y(:)); 

% numerical approximation ...
dudx_t = homg.tensor.IAX(quad.Dr, u);
dudy_t = homg.tensor.AIX(quad.Dr, u);

% error ...
err_x = norm(dudx - dudx_t);
err_y = norm(dudy - dudy_t);
disp(['2D errors: (' num2str(err_x) ', ' num2str(err_y) ')']);

%% 3D Case
oct  = homg.refel(3, order);

f3 = @(x,y,z)(x*y^order + z*x^order + y*z^order);
g3x = @(x,y,z)(y^order + order*z*x^(order-1));
g3y = @(x,y,z)(order*x*y^(order-1) + z^order);
g3z = @(x,y,z)(x^order + order*y*z^(order-1));

[X, Y, Z] = ndgrid(oct.r, oct.r, oct.r);

u    = arrayfun(f3,   X(:), Y(:), Z(:)); 
dudx = arrayfun(g3x,  X(:), Y(:), Z(:)); 
dudy = arrayfun(g3y,  X(:), Y(:), Z(:));
dudz = arrayfun(g3z,  X(:), Y(:), Z(:));

% numerical approximation ...
dudx_t = homg.tensor.IIAX(oct.Dr, u);
dudy_t = homg.tensor.IAIX(oct.Dr, u);
dudz_t = homg.tensor.AIIX(oct.Dr, u);

err_x = norm(dudx - dudx_t);
err_y = norm(dudy - dudy_t);
err_z = norm(dudz - dudz_t);

disp(['3D errors: (' num2str(err_x) ', ' num2str(err_y) ', ' num2str(err_z) ')']);