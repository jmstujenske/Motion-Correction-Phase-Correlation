function [fixed,dx]=correct_bidi_across_x(data,n_ch,whichch)
bg_pix=50; %radius of gaussian blur for background subtraction
upsample=10;
max_shift=1;
lags=-100:100;
smooth_lag=25;
edge_remove=10;
[Ly,Lx,nFrames]=size(data);
grn=data(:,:,whichch:n_ch:end);
xcor_dat=zeros(length(lags),Lx*upsample);
upsize{whichch}=imresize(grn,[Ly Lx*upsample]);
img=mean(upsize{whichch},3);
img2=max(img./imgaussfilt(img,bg_pix)-1,0);
data1=img2(1:2:end,:);
data2=img2(2:2:end,:);
min_length=min(size(data1,1),size(data2,1));
data1=data1(1:min_length,:,:);
data2=data2(1:min_length,:,:);
data1=double(data1)-mean(data1,2);
data2=double(data2)-mean(data2,2);
in=0;
for lag=lags
    in=in+1;
    data_shift=data1(edge_remove+1:end-edge_remove,min(max(1+lag:end+lag,1),Lx*upsample),:);
    xcor_dat(in,:)=mean(data_shift.*data2(edge_remove+1:end-edge_remove,:,:),[1 3]);
end
xcor_dat_smooth=conv2(xcor_dat,ones(1,smooth_lag*upsample)/(smooth_lag*upsample),'same');
[~,in]=max(xcor_dat_smooth);
dx=lags(in);
dx(abs(dx)>=upsample*max_shift)=NaN;
dx(isnan(dx))=interp1(find(~isnan(dx)),dx(~isnan(dx)),find(isnan(dx)),'nearest','extrap');
dx=conv(dx,gausskernel(smooth_lag*20,smooth_lag*2),'same');
upsize{whichch}=upsize{whichch}(1:2:end,:,:);
[Ly Lx nFrames]=size(upsize{whichch});
% Ly=Ly/2;
idy = repmat((1:Ly)', 1, Lx);
idx = repmat((1:Lx),  Ly, 1);
dx=round(dx);
% apply offsets to indexing
DX = dx + idx;

xyInvalid = DX<0 | DX>Lx-1;

DX(xyInvalid) = 0;

%     DX = mod(DX, Lx);
%     DY = mod(DY-1, Ly) + 1;
%
ind = max(idy + DX * Ly,1);

for rep=setdiff(1:n_ch,whichch)
    upsize{rep}=imresize(data(1:2:end,:,rep:n_ch:end),[Ly Lx]);
end
fixed=data;
for rep=1:n_ch
    upsize_data1=upsize{rep};
    for i=1:nFrames
        subpix_fix=upsize_data1(ind+(Ly*Lx*(i-1)));
        fixed(1:2:end,:,n_ch*(i-1)+rep)=imresize(subpix_fix,[Ly Lx/10]);
    end
end