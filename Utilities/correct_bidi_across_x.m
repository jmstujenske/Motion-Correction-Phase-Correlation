function [fixed,dx]=correct_bidi_across_x(data,n_ch,whichch,lowmemory,use_poly_fit,preserve_correlation_map,subpixel)
%[fixed,dx]=correct_bidi_across_x(data,n_ch,whichch)
%
%Corrects for bidirectional scanning offset, differentially along the scan
%axis (assumes x-axis)
%
%data - X by Y by (C*T) stack (channels interleaved)
%n_ch - now many channels
%whichch - which channel to based the correction on
%lowmemory - reduced memory mode (default: true); processes in chunks
%use_poly_fit - fit a quadratic fit rather than raw calculated values
%preserve_correlation_map - apply fixed dx value if the shift decreases
%from center out so that spatial correlation maps generated later on do not
%have abberations. Requires polyfit, so this will override use_poly_fit
%subpixl - how much to upsample in the x dimension
%
%OUTPUT:
%fixed - corrected stack of the same size as data
%dx - the calculated offsets along the x-axis
%
%requires the image processing toolbox
%
%
if nargin<7 || isempty(subpixel)
    subpixel=10;
end
if nargin<6 || isempty(preserve_correlation_map)
    preserve_correlation_map=true;
end
if nargin<4 || isempty(lowmemory)
    lowmemory=true;
end
if nargin<5 || isempty(use_poly_fit)
    use_poly_fit=true;
end
if nargin<2 || isempty(n_ch)
    n_ch=1;
end
if nargin<3 || isempty(whichch)
    whichch=1;
end
if preserve_correlation_map && ~use_poly_fit
    use_poly_fit=true;
    warning('Preserve correlation map on: Forcing poly fit on.');
end
bg_pix=50; %radius of gaussian blur for background correction
max_shift=4; %assume that the pixel correction should be less than 4 pixels in any given location;
%this is necessary to avoid weird edge effects from limited sample points --
%maybe another solution would be better
lags=-max_shift*subpixel:max_shift*subpixel;
smooth_lag=max(ceil(size(data,2)/100),1);
edge_remove=10;
[Ly,Lx_orig,nFrames]=size(data);
Lx=Lx_orig;
nFrames=nFrames/n_ch;
grn=data(:,:,whichch:n_ch:end);
xcor_dat=zeros(length(lags),Lx*subpixel);

%upsample along the x-axis
img=mean(grn,3);
img2=max(img./imgaussfilt(img,bg_pix)-1,0);
% img=imresize(img,[Ly Lx*subpixel*lag_subpix_fact]);

%correct for differences in fluorescence across the image, to equalize
%weighting
img2=max(img./imgaussfilt(img,bg_pix,'Padding','symmetric')-1,0);
img2=imresize(img2,[Ly Lx*subpixel]);

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
    data_shift=data1(edge_remove+1:end-edge_remove,min(max(1+lag:end+lag,1),Lx*subpixel),:);
    xcor_dat(in,:)=mean(data_shift.*data2(edge_remove+1:end-edge_remove,:,:),[1 3]);
end
%smooth across the x-axis
%editted 10/14/25 by JMS for more stable and smooth output for short
%videos:
[~,in]=max(xcor_dat);
dx=lags(in);
dx=dx/(subpixel);
dx(abs(dx)>=max_shift)=NaN;
% dx=movmean(movmedian(dx,smooth_lag*subpixel,'omitmissing','Endpoints','shrink'),smooth_lag*subpixel,'omitmissing','Endpoints','shrink');
dx=movmedian(dx,smooth_lag*subpixel,'omitmissing','Endpoints','shrink');

% dx(abs(dx-mean(dx))>=max_shift)=NaN;
if ~all(isnan(dx)) %%this shouldn't ever really happen, but just in case the dx are all NaN, then we will return the original data
dx(isnan(dx))=interp1(find(~isnan(dx)),dx(~isnan(dx)),find(isnan(dx)),'nearest','extrap');
else
    dx=zeros(size(dx));
    fixed=data;
    return;
end
m=nanmean(dx);
dx=conv2_symmetric(dx-m,gausskernel(smooth_lag*5*subpixel,smooth_lag/2*subpixel)','same')+m;
p=polyfit(1:length(dx),dx,2);
turningpoint=-p(2)/(2*p(1));
dx_turningpoint=polyval(p,turningpoint);

if use_poly_fit
dx=polyval(p,1:length(dx));
end

if preserve_correlation_map
    if (p(1)>0 && dx_turningpoint<0) || (p(1)<0 && dx_turningpoint>0)
    
        warning('Perserve Correlation Map On: Forcing single dx value.')
        % [fixed,dx]=correct_bidi_single(data,n_ch,whichch,lowmemory);
        % return;
        dx=m;
    end
end
if isscalar(dx) || preserve_correlation_map
dx=imresize(dx,[1 Lx_orig]);
fixed=apply_bidi_correction_direct(data,dx,n_ch,lowmemory,1);
else
fixed=apply_bidi_correction_direct(data,dx,n_ch,lowmemory,subpixel);
end
end

% function f = triangleResampling(x)
%     f = (1 - abs(x)) .* (abs(x) <= 1);
% end