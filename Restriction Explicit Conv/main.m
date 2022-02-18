clear;
addpath('.\GenerateMasks')

size.ImageSizeX=256;size.ImageSizeY=256;%pixels
size.MaskSizeX=512;size.MaskSizeY=512;
size.FullSizeX=size.ImageSizeX+size.MaskSizeX-1;
size.FullSizeY=size.ImageSizeY+size.MaskSizeY-1;

%Lim Sensor Size
ratio=1;%R
size.SensorSizeX=2*round(size.MaskSizeX*ratio/2);size.SensorSizeY=2*round(size.MaskSizeY*ratio/2);%pixels

% %Lim Pixel Pitch
% rate=4;%R_SR
% size.SensorSizeX=floor(size.FullSizeX/rate)+1;size.SensorSizeY=floor(size.FullSizeY/rate)+1;%pixels


%image preprocess
image=double(imread('cameraman.tif'));
image=imresize(image,[size.ImageSizeX,size.ImageSizeY]);
image=image./max(max(image));


%mask patten
mask=GenerateFZAMask(size.MaskSizeX,size.MaskSizeY);

%captured image
% 1. nonlinear restrictions
%   1.1. noise
NoiseLevel=0.00;
bOrigin=conv2(image,mask);
bOrigin=bOrigin+NoiseLevel*mean(bOrigin,"all")*rand(size.FullSizeX,size.FullSizeY);%FullSize

%   1.2. Limited sensor bit depth
dynamicRangeExp=12;
dynamicRange=2^dynamicRangeExp;
bOriginMax=max(max(bOrigin));
bOrigin=round(bOrigin/bOriginMax*(dynamicRange))*bOriginMax/dynamicRange;

% 2. linear restrictions
%   2.1. Limited sensor Size
L=@(x) x(((size.FullSizeX+1)/2-size.SensorSizeX/2)+1:((size.FullSizeX+1)/2+size.SensorSizeX/2),((size.FullSizeY+1)/2-size.SensorSizeY/2)+1:((size.FullSizeY+1)/2+size.SensorSizeY/2));
LT=@(x) padarray(padarray(x,[size.FullSizeX-size.SensorSizeX+1,size.FullSizeY-size.SensorSizeY+1]./2-1,0,'post'),[size.FullSizeX-size.SensorSizeX+1,size.FullSizeY-size.SensorSizeY+1]./2,0,'pre');

% %   2.2. Limited pixel pitch (super-resolution)
% L=@(x) sr(x,rate);
% LT=@(x) srT(x,rate,size.FullSizeX,size.FullSizeY);

b=L(bOrigin);



%% ADMM

opts.mu1=1e-5;opts.mu2=5e1;opts.mu3=5e1;
opts.tau=5e-1;
MaxIters=300;

[vRecord,residualRecord,optsRecord]=ADMM(mask,b,opts,size,MaxIters,LT);
RecImg=squeeze(vRecord(MaxIters,:,:));

figure
subplot(1,2,1)
imagesc(RecImg)
colormap gray
box on
set(gca,'xtick',[],'ytick',[])
subplot(1,2,2)
imagesc(image)
colormap gray
box on
set(gca,'xtick',[],'ytick',[])

PSNR=zeros(MaxIters,1);
MSE=zeros(MaxIters,1);
SSIM=zeros(MaxIters,1);


for iter=1:MaxIters
RecImg=squeeze(vRecord(iter,:,:));
PSNR(iter)=psnr(RecImg,image);
MSE(iter)=sum((RecImg-image).^2,'all')/(size.ImageSizeX*size.ImageSizeY);
SSIM(iter)=ssim(RecImg,image);
end


figure
subplot(2,2,1)
plot(PSNR)
xlabel('Iterations');ylabel('PSNR')
box on;
title('psnr');
subplot(2,2,2)
plot(MSE)
xlabel('Iterations');ylabel('MSE')
box on;
subplot(2,2,3)
plot(SSIM)
xlabel('Iterations');ylabel('SSIM')
box on;
title('SSIM');
subplot(2,2,4)
imagesc(RecImg)
colorbar
colormap gray


function y=sr(x)
y=x(1:4:end,1:4:end);
end

function y=srT(x,rate,FullSizeX,FullSizeY)
y=zeros(FullSizeX,FullSizeY);
y(1:rate:end,1:rate:end)=x;
end


