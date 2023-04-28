function [fixed,dx]=correct_bidi_across_x(data,n_ch,whichch,lowmemory)
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
if nargin<4 || isempty(lowmemory)
    lowmemory=true;
end
bg_pix=50; %radius of gaussian blur for background correction
upsample=10;
max_shift=1;
lags=-100:100;
smooth_lag=25;
edge_remove=10;
[Ly,Lx,nFrames]=size(data);
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
data1=data1(1:min_length,:,:);
data2=data2(1:min_length,:,:);
data1=double(data1)-mean(data1,2);
data2=double(data2)-mean(data2,2);
in=0;

%calculate cross-correlation across the x-axis
for lag=lags
    in=in+1;
    data_shift=data1(edge_remove+1:end-edge_remove,min(max(1+lag:end+lag,1),Lx*upsample),:);
    xcor_dat(in,:)=mean(data_shift.*data2(edge_remove+1:end-edge_remove,:,:),[1 3]);
end
%smooth across the x-axis
xcor_dat_smooth=conv2(xcor_dat,ones(1,smooth_lag*upsample)/(smooth_lag*upsample),'same');
[~,in]=max(xcor_dat_smooth);
dx=lags(in);

%remove outliers -- usually it is the edges of the image due to unavoidable
%edge effects; interpolate with nearest neighbor interpolation
dx(abs(dx)>=upsample*max_shift)=NaN;
dx(isnan(dx))=interp1(find(~isnan(dx)),dx(~isnan(dx)),find(isnan(dx)),'nearest','extrap');
dx=conv(dx,gausskernel(smooth_lag*20,smooth_lag*2),'same');

[Ly,Lx]=size(img);
Ly=ceil(Ly/2);

%%subpixel-registration
%
%setup the conversion
idy = repmat((1:Ly)', 1, Lx);
idx = repmat((1:Lx),  Ly, 1);
dx=round(dx);

% apply offsets to indexing
DX = dx + idx;

xyInvalid = DX<0 | DX>Lx-1;

DX(xyInvalid) = 0;

% ind = max(idy + DX * Ly,1);

%upsample the other channels
fixed=data;
if lowmemory
    batch_size=16;
else
    batch_size=Ly;
end
n_batches=ceil(Ly/batch_size);
for rep=1:n_ch
    for batch=1:n_batches
        lines=1+(batch-1)*batch_size*2:2:min(batch_size*2*batch,Ly*2);
        this_batch_size=length(lines);
        ins=(1:this_batch_size)+(batch-1)*batch_size;
        cur_ind = max(idy(ins,:)-(batch-1)*batch_size + DX(ins,:) * this_batch_size,1);
        upsized=imresize(data(lines,:,rep:n_ch:end),[this_batch_size Lx]);
%     OLD implementation: loop through every frame and apply the same correction to each
%     loop through every frame and apply the same correction to each
%         for i=1:nFrames
%             subpix_fix=upsized(cur_ind+(this_batch_size*Lx*(i-1)));
%             %downsample back to original sampling
%             fixed(lines,:,n_ch*(i-1)+rep)=imresize(subpix_fix,[Ly Lx/upsample]);
%         end
        inds=(cur_ind+this_batch_size*Lx*reshape((0:nFrames-1),1,1,[]));
        try
        subpix_fix=upsized(inds); %directly indexing the matrix is 30% faster for large matrices
        catch
            keyboard;
        end
        fixed(lines,:,rep:n_ch:end)=imresize(subpix_fix,[this_batch_size Lx/upsample]);
    end
end
end