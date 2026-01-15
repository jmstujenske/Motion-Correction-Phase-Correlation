function [fixed,dx]=correct_bidi_single(data,n_ch,whichch,lowmemory,subpixel)
%[fixed,dx]=correct_bidi_across_x(data,n_ch,whichch)
%
%Corrects for bidirectional scanning offset, differentially along the scan
%axis (assumes x-axis)
%
%data - X by Y by (C*T) stack (channels interleaved)
%n_ch - now many channels
%whichch - which channel to based the correction on
%lowmemory - reduced memory mode (default: true); processes in chunks
%
%OUTPUT:
%fixed - corrected stack of the same size as data
%dx - the calculated offsets along the x-axis
%
%requires the image processing toolbox
%
%
if nargin<5 || isempty(subpixel)
    subpixel=10;
end
if nargin<4 || isempty(lowmemory)
    lowmemory=true;
end
if nargin<2 || isempty(n_ch)
    n_ch=1;
end
if nargin<3 || isempty(whichch)
    whichch=1;
end
bg_pix=nanmean(size(data,1:2))/10; %radius of gaussian blur for background correction
upsample=subpixel;
max_shift=1.5; %assume that the pixel correction should be less than 1.5 pixels in any given location;
%this is necessary to avoid weird edge effects from limited sample points --
%maybe another solution would be better
lags=-100:100;
smooth_lag=50;
edge_remove=upsample;
[Ly,Lx_orig,nFrames]=size(data);
Lx=Lx_orig;
nFrames=nFrames/n_ch;
grn=data(:,:,whichch:n_ch:end);
xcor_dat=zeros(length(lags),Lx*upsample);

%upsample along the x-axis
img=mean(grn,3);
img=imresize(img,[Ly Lx*upsample]);

%correct for differences in fluorescence across the image, to equalize
%weighting
img2=max(img./imgaussfilt(img,bg_pix)-1,0);

%split into interleaving strips
data1=img2(1:2:end,:);
data2=img2(2:2:end,:);

%pre-initialize
min_length=min(size(data1,1),size(data2,1));
data1=double(data1(1:min_length,:,:));
data2=double(data2(1:min_length,:,:));
data1=zscore(data1,[],2);
data2=zscore(data2,[],2);
in=0;

%calculate cross-correlation across the x-axis
for lag=lags
    in=in+1;
    data_shift=data1(edge_remove+1:end-edge_remove,min(max(1+lag:end+lag,1),Lx*upsample),:);
    xcor_dat(in,:)=mean(data_shift.*data2(edge_remove+1:end-edge_remove,:,:),[1 3]);
end
%smooth across the x-axis
%editted 10/14/25 by JMS for more stable and smooth output for short
%videos:
[~,in]=max(xcor_dat);
% dx=movmean(movmedian(lags(in),smooth_lag*upsample,'omitmissing','Endpoints','shrink'),smooth_lag*upsample,'omitmissing','Endpoints','shrink');
%old solution:
% xcor_dat_smooth=movmedian(xcor_dat,smooth_lag*upsample,2,'omitmissing','Endpoints','shrink');
%[~,in]=max(xcor_dat_smooth);
% dx=lags(in);

%remove outliers -- usually it is the edges of the image due to unavoidable
%edge effects; interpolate with nearest neighbor interpolation
dx=lags(in);
dx(abs(dx-mean(dx))>=upsample*max_shift)=NaN;
if ~all(isnan(dx)) %%this shouldn't ever really happen, but just in case the dx are all NaN, then we will return the original data
dx(isnan(dx))=interp1(find(~isnan(dx)),dx(~isnan(dx)),find(isnan(dx)),'nearest','extrap');
else
    dx=zeros(size(dx));
    fixed=data;
    return;
end
dx=nanmean(dx);
% m=nanmean(dx);
% dx=conv2_symmetric(dx-m,gausskernel(smooth_lag*10,smooth_lag)','same')+m;
fixed=apply_bidi_correction_direct(data,dx,n_ch,lowmemory,subpixel);
end


function f = triangleResampling(x)
    f = (1 - abs(x)) .* (abs(x) <= 1);
end