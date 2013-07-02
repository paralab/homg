function coords=getGLLcoords(order, elems)
% function coords=getGLLcoords(order, elems)
% returns location of gll coordinates of order
% for elements in [0,1]

fac = 1.0/(2*elems);

% gll coordinates in [-1,1]
x = JacobiGL(0,0,order)';

x = (x + 1)*fac;


coords = [];
for i=1:elems
    y = x + (i-1)/elems;
    coords = [coords y(1:end-1)];
end

coords = [coords 1.0];
