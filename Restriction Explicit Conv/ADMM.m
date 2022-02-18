function [vRecord,residualRecord,optsRecord]=ADMM(mask,b,opts,size,MaxIters,LT)
%v:output

%h:convolutional kernel
%b:captrued information
%opts:parameters
%size
%LT:transposition of linear restriction matrix
%iters:iterations

%%
% Operators and PreCalculatedMatrixs
SoftThresh=@(x,tau) sign(x).*max(0,abs(x)-tau);

ft = @(x) fftshift(fft2(ifftshift(x)));
ift = @(x) fftshift(ifft2(ifftshift(x)));


Psi=@(x) cat(3,circshift(x,1,1)-x,circshift(x,1,2)-x);
PsiT=@(x) circshift(x(:,:,1),-1,1)-x(:,:,1)+circshift(x(:,:,2),-1,2)-x(:,:,2);

PsiTPsi=ft(padarray([0,-1,0;-1,4,-1;0,-1,0],[(size.FullSizeX-3)/2 (size.FullSizeY-3)/2],0,'both'));
h=zeros(size.FullSizeX,size.FullSizeY);
h((size.FullSizeX+1)/2-size.MaskSizeX/2+1:(size.FullSizeX+1)/2+size.MaskSizeX/2,(size.FullSizeY+1)/2-size.MaskSizeY/2+1:(size.FullSizeY+1)/2+size.MaskSizeY/2)=mask;
H=ft(h);
Pre_X_Inv=1./(LT(ones(size.SensorSizeX,size.SensorSizeY)) + opts.mu1+eps);
Pre_V_Inv=1./(opts.mu1*(abs(conj(H).*H))+opts.mu2*abs(PsiTPsi)+opts.mu3+eps);

M=@(x) real(ift(H.*ft(x)));
MT=@(x) real(ift(conj(H).*ft(x)));

Full2Image=@(x) x(((size.FullSizeX+1)/2-size.ImageSizeX/2):((size.FullSizeX+1)/2+size.ImageSizeX/2-1),((size.FullSizeY+1)/2-size.ImageSizeY/2):((size.FullSizeY+1)/2+size.ImageSizeY/2)-1);
%% Initialize
x=zeros(size.FullSizeX,size.FullSizeY); xi=x;
v=zeros(size.FullSizeX,size.FullSizeY);
w=zeros(size.FullSizeX,size.FullSizeY); rho=w;
u=zeros(size.FullSizeX,size.FullSizeY,2); eta=u;

vRecord=zeros([MaxIters,size.ImageSizeX,size.ImageSizeY]);


for iter=1:MaxIters
%Update u
u0=u;
u=SoftThresh(Psi(v) + eta/opts.mu2, opts.tau/opts.mu2);
%Update x
x0=x;
x=Pre_X_Inv.*(xi+opts.mu1*M(v)+LT(b));
%Update w
w0=w;
w=rho./opts.mu3 + v;
w(w<0)=0;w(w>1)=1;
%Update v
r=(opts.mu3*w - rho)+PsiT(opts.mu2*u - eta) + MT(opts.mu1*x - xi);
v=real(ift(Pre_V_Inv.*ft(r)));


%Update xi,pho,eta
xi=xi + opts.mu1*(M(v) - x);
eta=eta+opts.mu2*(Psi(v)-u);
rho=rho+opts.mu3*(v-w);


%Update penalty parameters
res.s1=norm(-opts.mu1*MT(x-x0));
res.s2=norm(-opts.mu2*PsiT(u-u0));
res.s3=norm(-opts.mu3*(w-w0));

res.r1=norm((M(v) - x));
res.r2=(Psi(v)-u);res.r2=norm(res.r2(:));
res.r3=norm((v-w));


%storage
vRecord(iter,:,:)=Full2Image(v);
residualRecord(iter)=res;
optsRecord(iter)=opts;

%output
fprintf(['iters=',num2str(iter),'  r1=',num2str(res.r1),'  s1=',num2str(res.s1),'  r2=',num2str(res.r2),'  s2=',num2str(res.s2),'  r3=',num2str(res.r3),'  s3=',num2str(res.s3),'\n'])

end

end

function pho=penaltyParas(pho0,r,s,tol)
if r>tol*s
    pho=pho0*1.1;
elseif s>tol*r
    pho=pho0/1.1;
else
    pho=pho0;
end
end




