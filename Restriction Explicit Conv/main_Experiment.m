clear;

addpath('.\GenerateMasks')
image = im2double(imread('THU.png'));
image=imresize(image,[526,526]);

size.ImageSizeX=526;size.ImageSizeY=526;%pixels
size.MaskSizeX=2424;size.MaskSizeY=2424;
size.FullSizeX=size.ImageSizeX+size.MaskSizeX-1;
size.FullSizeY=size.ImageSizeY+size.MaskSizeY-1;



% %Lim Sensor Size
% ratio=1;%R
% size.SensorSizeX=2*round(size.MaskSizeX*ratio/2);size.SensorSizeY=2*round(size.MaskSizeY*ratio/2);%pixels

%Lim Pixel Pitch
rate=4;%R_SR
size.SensorSizeX=floor(size.FullSizeX/rate)+1;size.SensorSizeY=floor(size.FullSizeY/rate)+1;%pixels



%mask patten
mask=GenerateFZAExpMask(size.MaskSizeX,size.MaskSizeY);

%captured image

Im = im2double(imread('raw_THU.png'));
bOrigin=Im(1700-1475:1700+1475-1,1760-1475:1760+1475-1);


%   2.1. Limited sensor Size
% L=@(x) x(((size.FullSizeX+1)/2-size.SensorSizeX/2)+1:((size.FullSizeX+1)/2+size.SensorSizeX/2),((size.FullSizeY+1)/2-size.SensorSizeY/2)+1:((size.FullSizeY+1)/2+size.SensorSizeY/2));
% LT=@(x) padarray(padarray(x,[size.FullSizeX-size.SensorSizeX+1,size.FullSizeY-size.SensorSizeY+1]./2,0,'post'),[size.FullSizeX-size.SensorSizeX+1,size.FullSizeY-size.SensorSizeY+1]./2-1,0,'pre');

%   2.2. Limited pixel pitch (super-resolution)
L=@(x) sr(x,rate);
LT=@(x) srT(x,rate,size.FullSizeX,size.FullSizeY);

b=L(bOrigin);




%% ADMM

%THU
opts.mu1=5e-6;opts.mu2=5e1;opts.mu3=5e0;
opts.tau=10;

MaxIters=50;

[vRecord,residualRecord,optsRecord]=ADMM(mask,b,opts,size,MaxIters,LT);
RecImg=squeeze(vRecord(MaxIters,:,:));

figure
subplot(1,2,1)
imagesc(RecImg)
colormap gray
box on
colorbar
set(gca,'xtick',[],'ytick',[])
subplot(1,2,2)
imagesc(image)
colormap gray
box on
colorbar
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


function y=sr(x,rate)
y=x(1:rate:end,1:rate:end);
end

function y=srT(x,rate,FullSizeX,FullSizeY)
y=zeros(FullSizeX,FullSizeY);
y(1:rate:end,1:rate:end)=x;
end



