function mask=GenerateFZAExpMask(Nx,Ny)
dp = 0.0038;      % pixel pitch
r0=0.325*1.01;
x=-Nx*dp/2:dp:(Nx-1)*dp/2;    y=-Ny*dp/2:dp:(Ny-1)*dp/2;  
[X,Y]=meshgrid(x,y);
mask= 0.5*(1 + sign(cos(pi*(X'.^2+Y'.^2)/r0^2)));
% mask(X.^2+Y.^2>(Nx*dp/2)^2)=1;
end
