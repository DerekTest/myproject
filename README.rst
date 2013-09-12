My First Projection
===================

Step 1.
-------

Solving the Poission equation of :math:`u^*` 

Step 2.
-------

Project the inter velocity into divergence space

function my_test


nx = 5;
ny = 5;
hx = 1/nx;
hy = 1/ny;

hxc = 0.5*(hx(2:end)+hx(1:end-1));
hyc = 0.5*(hy+hy);
A = getLaplcain(nx,ny,hx,hyc);

full(A)*max(hx).^2


% u: d2u/dx2 ; v: d2v/dy2
% Dp = spdiags(hx(2:end)',0,nx-1,nx-1)\spdiags([-1 0;ones(nx-3,1)*[-1,1];0 1],[0,1],nx-1,nx-1);
%
% Dn = spdiags(hx(1:end-1)',0,nx-1,nx-1)\spdiags([-1,1;ones(nx-3,1)*[-1,1];0,1],[-1,0],nx-1,nx-1);
%
% Ku = spdiags(hx(2:end)'+hx(1:end-1)',0,nx-1,nx-1)\ (2*(Dp-Dn));
% full(Ku*max(hx).^2)


% vertKu = spdiags([1 -1 0;ones(ny-2,1)*[1 -2 1];0 -1 1],[-1,0,1],ny,ny) ./ hy.^2;
% Lu = kron(speye(ny),horzKu) + kron(vertKu,speye(nx-1));

% figure(1)
% mesh(X,Y,zeros(size(X))); hold on
% plot(D/2*cos(linspace(0,2*pi,101)),D/2*sin(linspace(0,2*pi,101)),'b-');
% plot(xy(fd(xy)<=eps,1),xy(fd(xy)<=eps,2),'r.')
% hold off
% axis equal
% axis(domain)
% view(2)


end



function Ku = getMatrix(n,h)

if isscalar(h)
    h = h.*ones(n-1,1);
end

Dp = spdiags(h(2:end),0,n-1,n-1)\spdiags([-1,0;ones(n-3,1)*[-1,1];0,1],[0,1],n-1,n-1);
Dn = spdiags(h(1:end-1),0,n-1,n-1)\spdiags([-1,1;ones(n-3,1)*[-1,1];0,1],[-1,0],n-1,n-1);

Ku = spdiags(h(2:end)+h(1:end-1),0,n-1,n-1)\ (2*(Dp-Dn));
end

function Ku = getMatrix_stagg(n,h)

if isscalar(h)
    h = h.*ones(n-1,1);
end

%
h = [h(1);h;h(end)];

Dp = spdiags(h(2:end),0,n,n)\spdiags([-1,0;ones(n-2,1)*[-1,1];0,1],[0,1],n,n);
Dn = spdiags(h(1:end-1),0,n,n)\spdiags([-1,0;ones(n-2,1)*[-1,1];0,1],[-1,0],n,n);

Ku = spdiags(h(2:end)+h(1:end-1),0,n,n)\ (2*(Dp-Dn));
end

function A = getLaplcain(nx,ny,hx,hyc)

if isscalar(hx)
    hx = hx.*ones(nx,1);
end

if isscalar(hyc)
    hyc = hyc.*ones(ny-1,1);
end

u_Dpx = spdiags(hx(2:end),0,nx-1,nx-1)\spdiags([-1,0;ones(nx-3,1)*[-1,1];0,1],[0,1],nx-1,nx-1);
u_Dnx = spdiags(hx(1:end-1),0,nx-1,nx-1)\spdiags([-1,1;ones(nx-3,1)*[-1,1];0,1],[-1,0],nx-1,nx-1);

u_Kx = spdiags(hx(2:end)+hx(1:end-1),0,nx-1,nx-1)\ (2*(u_Dpx-u_Dnx));


hyc_append = [hyc(1);hyc;hyc(end)];

u_Dpy = spdiags(hyc_append(2:end),0,ny,ny)\spdiags([-1,0;ones(ny-2,1)*[-1,1];0,1],[0,1],ny,ny);
u_Dny = spdiags(hyc_append(1:end-1),0,ny,ny)\spdiags([-1,0;ones(ny-2,1)*[-1,1];0,1],[-1,0],ny,ny);

u_Ky = spdiags(hyc_append(2:end)+hyc_append(1:end-1),0,ny,ny)\ (2*(u_Dpy-u_Dny));
A = kron(speye(ny),u_Kx) + kron(u_Ky,speye(nx-1));
end

function varargout = getMesh

D = 0.1;
domain = [-1,1,-1, 1]; %[xmin,xmax,ymin,ymax]

nx1 = 85; nx2 = 60; nx3 = 105;  % total: 250
ny1 = 50; ny2 = 60; ny3 = 50;   % total: 160

nx = nx1+ny2+nx3;
ny = ny1+ny2+ny3;

x1 = linspace(domain(1),-D,nx1+1);
x2 = linspace(-D,D,nx2+1);
x3 = linspace(D,domain(2), nx3+1);
x = unique([x1,x2,x3]);

y1 = linspace(domain(3),-D,ny1+1);
y2 = linspace(-D,D,ny2+1);
y3 = linspace(D,domain(4), ny3+1);
y = unique([y1,y2,y3]);

hx = diff(x);
hy = diff(y);

xc = (x(2:end) + x(1:end-1))/2;
yc = (y(2:end) + y(1:end-1))/2;

hxc = diff(xc);
hyc = diff(yc);
end


Step 3.
-------

Update pressure