function dreg=reg2P_standalone(data,mimg,kriging,numBlocks,n_ch,whichch)
%dreg=reg2P_standalone(data,mimg,kriging,numBlocks,n_ch,whichch)
%data - X by Y by (C*T) frame stack
%mimg - template image (default: 1000 frame average)
%kriging - whether to use kriging or not
%numBlocks - default is [32 1]
%n_ch - how many channels in data
%whichch - which channel to motion correct based on
%
%Based on solution from Suite2p Matlab version, now made as a standable
%implementation
%https://github.com/cortex-lab/Suite2P
%
%Please cite the original authors
%

if nargin<6 || isempty(whichch)
    whichch=1;
end
if nargin<5 || isempty(n_ch)
    n_ch=whichch;
end
if nargin<3 ||isempty(kriging)
    kriging=false;
end
pad=15;
class_data=class(data);

data=pad_expand(data,pad);
data_orig=data;

data=data(:,:,whichch:n_ch:end);

data=single(data);



[Ly,Lx,nFrames] = size(data);

if nargin<2 || isempty(mimg)
    mimg=gen_template(data,min(1000,nFrames));
%     mimg=pad_expand(mimg,pad);
    mimg=single(mimg);
else
    mimg=pad_expand(mimg,pad);
    mimg=single(mimg);
end
mimg=imgaussfilt(mimg,.5);%-1.1*imgaussfilt(mimg,2);
mimg(mimg<0)=0;
if nargin<4 || isempty(numBlocks)
    numBlocks = [32 1];
elseif length(numBlocks)==1
    numBlocks=repmat(numBlocks,1,2);
end

% fraction of block of total image
bfrac     = 1./max(2,ceil((numBlocks-3)/3));
bfrac(numBlocks==1) = 1;
% blockFrac = getOr(ops, {'blockFrac'}, bfrac);
bpix      = ceil(bfrac .* [Ly Lx]);

% pixel overlap of blocks
% pixoverlap    = [];
% pixoverlap = round((bpix.*numBlocks-[Ly Lx])./(numBlocks-2));

% edges of blocks
yB        = linspace(0, Ly, numBlocks(1)+1);
yB        = round((yB(1:end-1) + yB(2:end)) / 2);

xB        = linspace(0, Lx, numBlocks(2)+1);
xB        = round((xB(1:end-1) + xB(2:end)) / 2);

% compute indices of blocks
ib        = 0;
yBL=cell(prod(numBlocks),1);
xBL=yBL;
for iy = 1:numBlocks(1)
%     if iy == numBlocks(1)
%             yB(iy)  = Ly - floor(bpix(1)/2);
%     elseif iy == 1
%         yB(iy)  = floor(bpix(1)/2) + 1;
%     end
    if numBlocks(2) > 1
        for ix = 1:numBlocks(2)
            ib = ib+1;
%             if ix == numBlocks(2)
%                 xB(ix)  = Lx - floor(bpix(2)/2);
%             elseif ix == 1
%                 xB(ix)  = floor(bpix(2)/2) + 1;
%             end

            yBL{ib} = max(1,yB(iy)-floor(bpix(1)/2)) : ...
                min(Ly,yB(iy)+floor(bpix(1)/2));
            xBL{ib} = max(1,xB(ix)-floor(bpix(2)/2)) : ...
                min(Lx,xB(ix)+floor(bpix(2)/2));
        end
    else
        ib = ib+1;
        yBL{ib} = max(1,yB(iy)-floor(bpix(1)/2)) : ...
            min(Ly,yB(iy)+floor(bpix(1)/2));
        xBL{ib} = 1:Lx;
    end
end
nblocks = ib;
% compute smoothing mask xyMask
% gaussian centered on block with width 2/3 the size of block
% (or specified by user as smoothBlocks)
% sT(1)        = mean(cellfun(@length,yBL)) *2/3;
% sT(2)        = mean(cellfun(@length,yBL)) *2/3;
% sT=[30 30];
% sT = max(10, sT);

sT(1)        = mean(diff(yB)) * 2/3;
sT(2)        = mean(diff(xB)) * 2/3;
sT = max(10, sT);
xyMask = zeros(Ly, Lx, nblocks, 'single');
ib=0;
for iy = 1:numBlocks(1)
    if numBlocks(2) > 1
        for ix = 1:numBlocks(2)
            ib=ib+1;
            if numBlocks(1)>1
                gausy = exp(-((1:Ly)' - yB(iy)).^2 / (2*sT(1).^2));
            else
                gausy=ones(Ly,1);
            end
            gausx = exp(-((1:Lx)' - xB(ix)).^2 / (2*sT(2).^2));
            xyMask(:, :, ib) = gausy * gausx';
        end
    else
        ib=ib+1;
        if numBlocks(1)>1
            gausy = exp(-((1:Ly)' - yB(iy)).^2 / (2*sT(1).^2));
        else
            gausy = ones(Ly,1);
        end
        xyMask(:, :, ib) = repmat(gausy, 1, Lx);
    end
end

xyMask = xyMask./repmat(sum(xyMask, 3), 1, 1, size(xyMask, 3));
xyMask = reshape(xyMask, Ly*Lx, nblocks);
xyMask(isnan(xyMask))=0;
% xyMask    = xyMask;

indframes=1:size(data,3);
ds = zeros(nblocks,numel(indframes),2,'double');
Corr = zeros(numel(indframes), nblocks,'double');
%     for ib = 1:nblocks
%         % collect ds
%         ops1{i}.mimg = ops1{i}.mimgB{ib};
% 	if kriging
% 	  [ds(:,:,ib), Corr(:,ib)]  = ...
% 	      regoffKriging(data(ops1{i}.yBL{ib},ops1{i}.xBL{ib},indframes),ops1{i}, 0);
% 	else
% 	  [ds(:,:,ib), Corr(:,ib)]  = ...
% 	      regoffLinear(data(ops1{i}.yBL{ib},ops1{i}.xBL{ib},indframes),ops1{i},0);
% 	end
%     end
%     if j==1
%         ds(1,:,:) = 0;
%     end

%         DS          = [];
%         CorrFrame   = [];
%         mimg1       = zeros(Ly, Lx);
mimgB = cell(prod(numBlocks),1);
for ib = 1:numBlocks(1)*numBlocks(2)
    mimgB{ib} = mimg(yBL{ib},xBL{ib});
end



% refImg = mimg;
subpixel = 10;
useGPU = false;
phaseCorrelation = true;
% maximum shift allowed
maxregshift = 30;
% slope on taper mask preapplied to image. was 2, then 1.2
maskSlope   = 2;
% SD pixels of gaussian smoothing applied to correlation map (MOM likes .6)
smoothSigma = 1.15;


% if subpixel is still inf, threshold it for new method
subpixel = min(10, subpixel);
% data_smooth=imgaussfilt(data,1);
% data_smooth=convn(data,ones(3,1)/3,'same');
data_smooth=imgaussfilt(data,.5);%-1.1*imgaussfilt(data,2);
data_smooth(data_smooth<0)=0;
for ib=1:nblocks
    refImg=mimgB{ib};
    subdata=data_smooth(yBL{ib},xBL{ib},:);
    %     subdata=convn(subdata,reshape(gausskernel(10,2),1,1,[]),'same');
    [ly,lx,nFrames] = size(subdata);

    % Taper mask
    [ys, xs] = ndgrid(1:ly, 1:lx);
    ys = abs(ys - mean(ys(:)));
    xs = abs(xs - mean(xs(:)));
    mY      = max(ys(:)) - 4;
    mX      = max(xs(:)) - 4;
    maskMul = single(1./(1 + exp((ys - mY)/maskSlope)) ./(1 + exp((xs - mX)/maskSlope)));
    maskOffset = mean(refImg(:))*(1 - maskMul);
    % Smoothing filter in frequency domain
    hgx = exp(-(((0:lx-1) - fix(lx/2))/smoothSigma).^2);
    hgy = exp(-(((0:ly-1) - fix(ly/2))/smoothSigma).^2);
    hg = hgy'*hgx;
    fhg = real(fftn(ifftshift(single(hg/sum(hg(:))))));

    % fft of reference image
    eps0          = single(1e-10);
    cfRefImg = conj(fftn(refImg));
    if phaseCorrelation
        absRef   = abs(cfRefImg);
        cfRefImg = cfRefImg./(eps0 + absRef) .* fhg;
    end

    if useGPU
        batchSize = getBatchSize(ly*lx);
        eps0      = gpuArray(eps0);
        cfRefImg    = gpuArray(cfRefImg);
        maskMul = gpuArray(maskMul);
        maskOffset = gpuArray(maskOffset);
    else
        batchSize = 1000;
    end

    % allow max shifts +/- lcor
    lpad   = 3;
    lcorr  = min(maxregshift, floor(min(ly,lx)/2)-lpad);

    % only need a small kernel +/- lpad for smoothing
    [x1,x2] = ndgrid([-lpad:lpad]);
    xt = [x1(:) x2(:)]';
    if useGPU
        xt = gpuArray(single(xt));
    end

    if kriging
        % compute kernels for regression
        sigL     = .85; % kernel width in pixels
        Kx = kernelD(xt,xt,sigL*[1;1]);
        linds = -lpad:1/subpixel:lpad;
        [x1,x2] = ndgrid(linds);
        xg = [x1(:) x2(:)]';
        if useGPU
            xg = gpuArray(single(xg));
        end
        Kg = kernelD(xg,xt,sigL*[1;1]);
        Kmat = Kg/Kx;
    end

    % loop over batches
    % dv = zeros(nFrames, 2);
    % corr = zeros(nFrames, 1);

    nBatches = ceil(nFrames/batchSize);
    for bi = 1:nBatches
        fi = ((bi - 1)*batchSize + 1):min(bi*batchSize, nFrames);

        if useGPU
            batchData = gpuArray(single(subdata(:,:,fi)));
        else
            batchData = single(subdata(:,:,fi));
        end

        corrMap = fft2(bsxfun(@plus, maskOffset, bsxfun(@times, maskMul, batchData)));

        %keyboard;
        if phaseCorrelation
            corrMap = bsxfun(@times, corrMap./(eps0 + abs(corrMap)), cfRefImg);
        else
            corrMap = bsxfun(@times, corrMap, cfRefImg);
        end

        % compute correlation matrix
        corrClip = real(ifft2(corrMap));
        corrClip = fftshift(fftshift(corrClip, 1), 2);
        corrClipSmooth = corrClip;
        %% subpixel registration
        if subpixel > 1
            % kriging subpixel
            % allow only +/- lcorr shifts
            cc0         = corrClipSmooth(floor(ly/2)+1+(-lcorr:lcorr),...
                floor(lx/2)+1+(-lcorr:lcorr),:);
            [cmax,ii]   = max(reshape(cc0, [], numel(fi)),[],1);

            [iy, ix] = ind2sub((2*lcorr+1) * [1 1], ii);

            dl       = single(-lpad:1:lpad);
            if useGPU
                dl   = gpuArray(dl);
                ccmat = gpuArray.zeros(numel(dl), numel(dl), numel(fi), 'single');
            else
                ccmat = zeros(numel(dl), numel(dl), numel(fi), 'single');
            end
            mxpt        = [iy(:)+floor(ly/2) ix(:)+floor(lx/2)] - lcorr;
            for j = 1:numel(fi)
                % matrix +/- lpad surrounding max point
                ccmat(:, :, j) = corrClip(min(mxpt(j,1)+dl,size(corrClip,1)), min(mxpt(j,2)+dl,size(corrClip,2)), j);
            end
            %
            ccmat = reshape(ccmat,[], numel(fi));
            if kriging
                % regress onto subsampled grid
                ccb         = Kmat * ccmat;

                % find max of grid
                [cx,ix]     = max(ccb, [], 1);
                [ix11,ix21] = ind2sub(numel(linds)*[1 1],ix);
                mdpt        = floor(numel(linds)/2)+1;
                dv0         = bsxfun(@minus, ([ix11' ix21'] - mdpt)/subpixel + mxpt, ...
                    [floor(ly/2) floor(lx/2)]) - 1;
            else
                yshift      = xt(1,:) * ccmat;
                xshift      = xt(2,:) * ccmat;
                dv0         = bsxfun(@rdivide, [yshift' xshift'], sum(ccmat, 1)') + mxpt;
                dv0         = bsxfun(@minus, dv0, [floor(ly/2) floor(lx/2)] + 1);
                if isfinite(subpixel)
                    dv0 = round(dv0 * subpixel) / subpixel;
                end
                cx     = max(ccmat, [], 1);
            end
            ds(ib,fi,:) = gather_try(dv0);
            Corr(fi,ib)  = gather_try(cx);
            % otherwise just take peak of matrix
        else
            cc0     = corrClipSmooth(floor(ly/2)+1+[-lcorr:lcorr],floor(lx/2)+1+[-lcorr:lcorr],:);
            [cmax,iy]  = max(cc0,[],1);
            [cx, ix]   = max(cmax,[],2);
            iy = reshape(iy(sub2ind([size(iy,2) size(iy,3)], ix(:), (1:size(iy,3))')),...
                1, 1, []);

            dv0 = [iy(:)-lcorr ix(:)-lcorr]-1;
            ds(ib,fi,:)  = gather_try(dv0);
            Corr(fi,ib) = gather_try(cx(:));
        end

    end
end

[Ly, Lx, NT] = size(data_orig);

%%
% ds=dv;
% xyMask=true(size(data,1),size(data,2));
% smooth offsets across blocks by xyMask

%center rigid registration
dx = (xyMask * ds(:,:,2));
dy = (xyMask * ds(:,:,1));

% max_shift=[10 30];
% dx(abs(dx)>max_shift(2))=NaN;
% dy(abs(dy)>max_shift(1))=NaN;
% filt_size=5;
% while any(isnan(dx) | isnan(dy),1:3)
%     medfilt_dx=medfilt1(dx,filt_size,'omitnan','zeropad',3);
%     medfilt_dy=medfilt1(dy,filt_size,'omitnan','zeropad',3);
%     dx(isnan(dx))=medfilt_dx(isnan(dx));
%     dy(isnan(dy))=medfilt_dy(isnan(dy));
%     filt_size=filt_size+1;
% end


dx = reshape(dx, Ly, Lx, []);
dy = reshape(dy, Ly, Lx, []);

idy = repmat([1:Ly]', 1, Lx);
idx = repmat([1:Lx],  Ly, 1);

dreg = zeros(size(data), class_data);

%Remove edge effects for dark areas without cells
% upperlim_dxy=15;
% fill_fact=20;
% mask=logical(imclose(imopen(any(abs(dx)>upperlim_dxy,3) | any(abs(dy)>upperlim_dxy,3),strel('disk',floor(fill_fact/2))),strel('disk',fill_fact)));
% [x_bad,y_bad]=find(mask);
% [x,y]=find(~mask);
% dx=medfilt1(dx,5,[],3);
% dy=medfilt1(dy,5,[],3);
for i = 1:NT
    frame_num=ceil(i/n_ch);
    Im = data_orig(:,:,i);
    dx_i=dx(:,:,frame_num);
    dy_i=dy(:,:,frame_num);

dreg(:,:,i)=imwarp(cast(Im,class_data),cat(3,dx_i,dy_i));
end
dreg=dreg(pad+1:end-pad,pad+1:end-pad,:);


% if nargin > 2 && removeMean
%     dv = bsxfun(@minus, dv, mean(dv,1));
% end
function data=pad_expand(data,pad)
data_temp=data;
data=[data(pad+1:-1:2,:,:);data;data(end:-1:end-pad+1,:,:)];
data=[data(:,pad+1:-1:2,:) data data(:,end:-1:end-pad+1,:)];

data=imgaussfilt(data,floor(pad/2));
data(pad+1:end-pad,pad+1:end-pad,:)=data_temp;
% keyboard
% data(1:pad,1:pad,:)=[data(pad+:pad*2,pad+1:-1:2,:)