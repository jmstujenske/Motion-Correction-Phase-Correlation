function [dreg_full,shifts]=reg2P_standalone(data,mimg,kriging,numBlocks,n_ch,whichch,maxregshift,fs,quick,use_subpixel_reg)
%dreg=reg2P_standalone(data,mimg,kriging,numBlocks,n_ch,whichch)
%data - X by Y by (C*T) frame stack
%mimg - template image (default: 1000 frame average)
%kriging - whether to use kriging or not
%numBlocks - default is [32 1]
%n_ch - how many channels in data
%whichch - which channel to motion correct based on
%maxregshift - maximum shift (default: 30)
%Based on solution from Suite2p Matlab version, now made as a standable
%implementation
%https://github.com/cortex-lab/Suite2P
%
%Please cite the original authors
%

subpixel = 10;
useGPU = false;
phaseCorrelation = true;
% maximum shift allowed
% slope on taper mask preapplied to image. was 2, then 1.2
maskSlope   = 2;
% SD pixels of gaussian smoothing applied to correlation map (MOM likes .6)
smoothSigma = 1.15;

% if subpixel is still inf, threshold it for new method
subpixel = min(10, subpixel);

eps0 = single(1e-10);

if nargin<10 || isempty(use_subpixel_reg)
    use_subpixel_reg=true;
end
if nargin<9 || isempty(quick)
    quick=false;
end
if nargin<6 || isempty(whichch)
    whichch=1;
end
if nargin<5 || isempty(n_ch)
    n_ch=whichch;
end
if nargin<8 || isempty(fs)
    fs=30;
end
if nargin<7 || isempty(maxregshift)
    maxregshift=30;
end

if nargin<3 ||isempty(kriging)
    kriging=false;
end
if nargin<4 || isempty(numBlocks)
    numBlocks=[1 1];
end

        % function [y,x,m]=cross_corr_memory(block_size)
        %     if nargin<1 || isempty(block_size)
        %         block_size=5000;
        %     end
        %     [Ysize,Xsize,nf]=size(data_input,1:3);
        %     n_block = ceil(nf/block_size);
        %             cfRefImg = conj(fft2(single(refImg)));
        %     x=zeros(1,nf/n_ch);
        %     y=zeros(1,nf/n_ch);
        %     m=zeros(1,nf/n_ch);
        %     for i=1:n_block
        %         frames=whichch+(i-1)*block_size*n_ch:n_ch:min(i*block_size*n_ch,nf);
        %         corrMap=fft2(single(data_input(:,:,frames)));
        % 
        %         if phaseCorrelation
        %             absRef   = abs(cfRefImg);
        %             hgx = exp(-(((0:Xsize-1) - fix(Xsize/2))/smoothSigma).^2);
        %             hgy = exp(-(((0:Ysize-1) - fix(Ysize/2))/smoothSigma).^2);
        %             hg = hgy'*hgx;
        %             fhg = real(fftn(ifftshift(single(hg/sum(hg(:))))));
        %             cfRefImg = cfRefImg./(eps0 + absRef) .* fhg;
        %         end
        %         if phaseCorrelation
        %             corrMap = bsxfun(@times, corrMap./(eps0 + abs(corrMap)).*fhg, cfRefImg);
        %         else
        %             corrMap = bsxfun(@times, corrMap, cfRefImg);
        %         end
        % 
        %         % corrMap=corrMap.*cfRefImg;
        %         corrClip = real(ifft2(corrMap));
        %         corrClip = fftshift(fftshift(corrClip, 1), 2);
        %         Y_mid=floor(Ysize/2)+1;
        %         X_mid=floor(Xsize/2)+1;
        %         corrClip=corrClip(Y_mid+(-maxregshift:maxregshift),X_mid+(-maxregshift:maxregshift),:);
        %         corrClip=imresize(corrClip,subpixel);
        %         [m((frames-whichch)/n_ch+1),vals]=max(reshape(corrClip,[],size(corrClip,3)));
        %         [y_t,x_t]=ind2sub(size(corrClip,1:2),vals);
        %         y((frames-whichch)/n_ch+1)=(y_t-maxregshift*subpixel-1)/subpixel;
        %         x((frames-whichch)/n_ch+1)=(x_t-maxregshift*subpixel-1)/subpixel;
        %     end
        % end


nf=size(data,3);
if nf/n_ch<=1+fs*4
quick=false;
end
if quick
    if nargin<2 || isempty(mimg)
        % refImg=gen_mimg(data,whichch,n_ch,n_Frames_template);
        refImg=mean(data(:,:,whichch:n_ch:end),3);
    else
        refImg=mimg;
    end
    % fft_sz=5000;
        % [y,x,m]=cross_corr_memory(fft_sz);

% [ds,m] = register_blocks_fft_subpixel( ...
%     data(:,:,whichch:n_ch:end), refImg,[], ...
%     maxregshift, subpixel, ...
%     smoothSigma, maskSlope, ...
%     phaseCorrelation, kriging, ...
%     useGPU, eps0);
[ds_rigid,m] = register_blocks_fft_subpixel( ...
    data(:,:,whichch:n_ch:end), refImg,[], ...
    maxregshift, subpixel, ...
    smoothSigma, maskSlope, ...
    phaseCorrelation, kriging, ...
    useGPU, eps0,size(data,3),true);
    % xyMask_rigid = make_xyMask(size(data,1), size(data,2), [1 1]);
    % [y,x]=ds_to_dxy(xyMask_rigid,ds,size(data,1:2));
    [y,x]=ds_to_dxy([],ds_rigid,size(data,1:2));

            % x=ds(:,:,2);
            % y=ds(:,:,1);

        % x=(x-median(x));
        % y=(y-median(y));
                    m_smooth=sgolayfilt(double(m),3,floor(fs/2)*2+1);
            m_smooth=m_smooth-movmax(sgolayfilt(m_smooth,3,floor(fs/2)*2+1),[fs*3 0]);
            th=multithresh(m_smooth)*.8;
data=apply_rigid_dx(data,x,y,n_ch,use_subpixel_reg);
% ds_rigid=cat(3,y(:)',x(:)');
% if prod(numBlocks)==1
%     dreg_full=data;
%     shifts={[],ds_rigid};
%     return;
% end
movements=m_smooth<th | sqrt(x.^2+y.^2)>2; %only need to apply non-rigid to
% cases where the maximum correlation after rigid correction dips below
% expected, or large displacements, to be conservative
movements_pad=find(conv(movements,ones(1,max(ceil(fs/4),3)),'same')>0);
clear movements
movements_pad=(movements_pad(:)-1)*n_ch+(1:n_ch)*n_ch;
movements_pad=movements_pad(:);
if isempty(movements_pad)
    shifts=[];
    return;
end
data_sub=data(:,:,movements_pad);
else
    data_sub=data;
    movements_pad=1:size(data,3);
end


pad=0;
class_data=class(data_sub);
data_movtrunc_allch=data_sub;

data_sub=data_sub(:,:,whichch:n_ch:end);
data_sub=single(data_sub);
data_sub=pad_expand(data_sub,pad);

[Ly,Lx,nFrames] = size(data_sub);

if nargin<2 || isempty(mimg)
    if quick
        channel_idx=whichch:n_ch:nf;
        mimg=gen_template(data(:,:,channel_idx(setdiff(1:nf,movements_pad))),min(1000,nFrames));
        mimg=pad_expand(mimg,pad);
    else
        mimg=gen_template(data_sub,min(1000,nFrames));
    end
    mimg=single(mimg);
else
    mimg=pad_expand(mimg,pad);
    mimg=single(mimg);
end
mimg=imgaussfilt(mimg,.5);%-1.1*imgaussfilt(mimg,2);
mimg(mimg<0)=0;

% data_smooth=imgaussfilt(data,.5);%-1.1*imgaussfilt(data,2);
% clear data
data_sub(data_sub<0)=0;

[ds_rigid2] = register_blocks_fft_subpixel( ...
    data_sub, mimg,[], ...
    maxregshift, subpixel, ...
    smoothSigma, maskSlope, ...
    phaseCorrelation, kriging, ...
    useGPU, eps0,1000,false);
    % xyMask_rigid = make_xyMask(size(data,1), size(data,2), [1 1]);
    % [y,x]=ds_to_dxy(xyMask_rigid,ds,size(data,1:2));
    [y,x]=ds_to_dxy([],ds_rigid2,size(data,1:2));

            % x=ds(:,:,2);
            % y=ds(:,:,1);

        % x=(x-median(x));
        % y=(y-median(y));
data_sub=apply_rigid_dx(data_sub,x,y,n_ch,use_subpixel_reg);
if quick
ds_rigid(:,movements_pad,:)=ds_rigid(:,movements_pad,:)+ds_rigid2;
else
    ds_rigid=ds_rigid2;
end
clear ds_rigid2
if prod(numBlocks)==1
    dreg_full=data;
    shifts={[],ds_rigid};
    return;
end

if nargin<4 || isempty(numBlocks)
    numBlocks = [32 1];
elseif length(numBlocks)==1
    numBlocks=repmat(numBlocks,1,2);
end
nblocks = prod(numBlocks);

xyMask = make_xyMask(Ly, Lx, numBlocks);   

    ds = register_blocks_fft_subpixel( ...
    data_sub, mimg,numBlocks, ...
    maxregshift, subpixel, ...
    smoothSigma, maskSlope, ...
    phaseCorrelation, kriging, ...
    useGPU, eps0,5000,false);
    %loop through every block and remove outliers
NT = size(data_movtrunc_allch,3);
[Ly,Lx]=size(data_sub,1:2);
clear data_sub
% dreg = zeros([Ly Lx], class_data);
[dy,dx]=ds_to_dxy(xyMask,ds,[Ly Lx]);

if nargout>1
    if quick
        ds_full=repmat(ds_rigid,nblocks,1,1);
        ds_full(:,movements_pad,:)=ds_full(:,movements_pad,:)+ds;
        ds=ds_full;
    end
end
shifts={xyMask,ds};

clear xyMask ds
% for idx = 1:NT
%     frame_num=ceil(idx/n_ch);
%     Im = data_movtrunc_allch(:,:,idx);
%     Im=pad_expand(Im,pad);
%     dx_i=dx(:,:,frame_num);
%     dy_i=dy(:,:,frame_num);
%     dreg(:,:,idx)=imwarp(Im,cat(3,dx_i,dy_i),'linear');
% end
        if size(dx,1)==1 || size(dx,2)==1
            dreg=apply_rigid_dx(data_movtrunc_allch,dx,dy,n_ch,use_subpixel_reg);
        else
            dreg=apply_nonrigid_dx(data_movtrunc_allch,dx,dy,n_ch,use_subpixel_reg,pad);
        end

% dreg=apply_nonrigid_dx(data_movtrunc_allch,dx,dy,n_ch,false,pad);
dreg_full=zeros(size(data),class_data);

full_mov=reshape(movements_pad+(0:n_ch-1),[],1);
dreg_full(:,:,setdiff(1:nf,full_mov))=data(:,:,setdiff(1:nf,full_mov));
dreg_full(:,:,full_mov)=dreg;
end




function out=fft2_memory(X,block_size)
if nargin<2 || isempty(block_size)
    block_size=5000;
end
nf=size(X,3);
n_block = ceil(nf/block_size);
out=cell(1,n_block);
for i=1:n_block
    fs=1+(i-1)*block_size:min(i*block_size,nf);
    newval=fft2(X(:,:,fs));
out{i}=newval; %for some reason, trying to populate a matrix directly causes memory issues, so store data into cell array and then concatenate at the end
end
out=cat(3,out{:});
end

function out=apply_rigid_dx(X,dx,dy,n_ch,subpixel)
if nargin<5
    subpixel=false;
end
orig_x=1:size(X,2);
orig_y=1:size(X,1);
nf=size(X,3)/n_ch;
changes=(abs(dx)>.5 | abs(dy)>.5);
fs=1:nf;
fs=fs(changes);
out=X;

if ~subpixel
[x_coords,y_coords]=meshgrid(orig_x,orig_y);

for i=fs
    fs_allch=(i-1)*n_ch+(1:n_ch);
    x_coords_new=x_coords+double(dx(i));
y_coords_new=y_coords+double(dy(i));
x_coords_new=min(max(round(x_coords_new),1),size(X,2));
y_coords_new=min(max(round(y_coords_new),1),size(X,1));
ins_orig=(x_coords-1)*size(X,1)+y_coords+reshape(fs_allch-1,1,1,[])*prod(size(X,1:2));
ins=(x_coords_new-1)*size(X,1)+y_coords_new+reshape(fs_allch-1,1,1,[])*prod(size(X,1:2));
out(double(ins_orig))=X(double(ins));
end

else
    for i = fs
        for j=1:n_ch
            Im = X(:,:,(i-1)*n_ch+j);
            dx_i=repmat(double(dx(i)),size(X,1:2));
            dy_i=repmat(double(dy(i)),size(X,1:2));
            out(:,:,(i-1)*n_ch+j)=imwarp(Im,cat(3,dx_i,dy_i),'linear');
        end
    end
end
end


function out=apply_nonrigid_dx(X,dx,dy,n_ch,subpixel,pad)
if nargin<5
    subpixel=false;
end
if nargin<6
    pad=15;
end
orig_x=1:size(X,2)+pad*2;
orig_y=1:size(X,1)+pad*2;
nf=size(X,3)/n_ch;
changes=find(squeeze(any(abs(dx)>.5 | abs(dy)>.5,1:2)));
fs=1:nf;
% fs=fs(changes);
out=zeros(size(dx,1),size(dx,2),size(X,3),class(X));

if ~subpixel
[x_coords,y_coords]=meshgrid(orig_x,orig_y);
for i=fs
    fs_allch=(i-1)*n_ch+(1:n_ch);
    X_sub=X(:,:,fs_allch);
    X_sub=pad_expand(X_sub,pad);
    if ismember(i,changes)
    x_coords_new=x_coords+double(dx(:,:,i));
    y_coords_new=y_coords+double(dy(:,:,i));
    x_coords_new=min(max(round(x_coords_new),1),size(X_sub,2));
    y_coords_new=min(max(round(y_coords_new),1),size(X_sub,1));
    ins_orig=(x_coords-1)*size(X_sub,1)+y_coords+reshape(fs_allch-1,1,1,[])*prod(size(X_sub,1:2));
    ins=(x_coords_new-1)*size(X_sub,1)+y_coords_new+reshape(0:n_ch-1,1,1,[])*prod(size(X_sub,1:2));
    out(double(ins_orig))=X_sub(double(ins));
    else
        out(:,:,fs_allch)=X_sub;
    end
end

else
    for i = fs
        for j=1:n_ch
            Im = X(:,:,(i-1)*n_ch+j);
            Im = pad_expand(Im,pad);
                if ismember(i,changes)

            dx_i=dx(:,:,i);
            dy_i=dy(:,:,i);
            out(:,:,(i-1)*n_ch+j)=imwarp(Im,cat(3,dx_i,dy_i),'linear');
                else
                    out(:,:,(i-1)*n_ch+j)=Im;
                end
        end
    end
end
out=out(pad+1:end-pad,pad+1:end-pad,:);
end


function [ds,m] = register_blocks_fft_subpixel( ...
    data, mimg,numBlocks, ...
    maxregshift, subpixel, ...
    smoothSigma, maskSlope, ...
    phaseCorrelation, kriging, ...
    useGPU, eps0,batchSize,useSVD)

% data        : full movie [Y X T]
% mimgB       : cell array of reference images per block
% yBL, xBL    : cell arrays of y/x indices per block
% nFrames     : number of frames
% nblocks     : number of spatial blocks
% maxregshift : maximum allowed integer shift
% subpixel    : subpixel factor (1 = integer only)
% smoothSigma : frequency-domain smoothing width
% maskSlope   : taper mask slope
% phaseCorrelation : logical
% kriging     : logical
% useGPU      : logical
% eps0        : small constant for normalization
%
% ds          : [nblocks x nFrames x 2] shifts (y,x)
if nargin<13 || isempty(useSVD)
    useSVD=false;
end
if nargin<12 || isempty(batchSize)
    batchSize=1000;
end
[Ly,Lx,nFrames] = size(data);
if isempty(numBlocks)
    numBlocks=[1 1];
end
nblocks=prod(numBlocks);
[yBL, xBL] = compute_block_indices(Ly, Lx, numBlocks);
mimgB = cell(nblocks,1);
for ib = 1:nblocks
    mimgB{ib} = mimg(yBL{ib},xBL{ib});
end 

if useGPU
    eps0 = gpuArray(eps0);
end
nblocks=numel(yBL);
ds = zeros(nblocks, nFrames, 2, 'single');
m = zeros(nblocks,nFrames,'single');
ly_old = -1;
lx_old = -1;

for ib = 1:nblocks

    refImg = mimgB{ib};
    refImg = refImg - mean(refImg,'all');
    subdata = data(yBL{ib}, xBL{ib}, :);
    subdata = imgaussfilt(subdata, 0.5, 'Padding', 'symmetric');

    ly = numel(yBL{ib});
    lx = numel(xBL{ib});

    % --- geometry-dependent precomputation ---
    if ly ~= ly_old || lx ~= lx_old
        ly_old = ly;
        lx_old = lx;

        lpad  = 3;
        lcorr = min(maxregshift, floor(min(ly,lx)/2) - lpad);

        % local coordinate grid
        [x1, x2] = ndgrid(-lpad:lpad);
        xt = [x1(:) x2(:)]';
        if useGPU
            xt = gpuArray(single(xt));
        else
            xt = single(xt);
        end

        if kriging
            sigL = 0.85;
            Kx = kernelD(xt, xt, sigL*[1;1]);

            linds = -lpad:1/subpixel:lpad;
            [x1g, x2g] = ndgrid(linds);
            xg = [x1g(:) x2g(:)]';

            if useGPU
                xg = gpuArray(single(xg));
            else
                xg = single(xg);
            end

            Kg   = kernelD(xg, xt, sigL*[1;1]);
            Kmat = Kg / Kx;
            Kmat = bsxfun(@rdivide, Kmat, sum(Kmat,2));
        end

        % taper mask
        [ys, xs] = ndgrid(1:ly, 1:lx);
        ys = abs(ys - mean(ys(:)));
        xs = abs(xs - mean(xs(:)));

        mY = max(ys(:)) - 4;
        mX = max(xs(:)) - 4;

        maskMul = single( ...
            1 ./ (1 + exp((ys - mY)/maskSlope)) .* ...
            1 ./ (1 + exp((xs - mX)/maskSlope)) );

        % frequency-domain smoothing kernel
        hgx = exp(-(((0:lx-1) - fix(lx/2)) / smoothSigma).^2);
        hgy = exp(-(((0:ly-1) - fix(ly/2)) / smoothSigma).^2);
        hg  = hgy' * hgx;
        fhg = real(fftn(ifftshift(single(hg / sum(hg(:))))));
    end

    maskOffset = mean(refImg(:)) * (1 - maskMul);

    % FFT of reference
    cfRefImg = conj(fftn(refImg));
    if phaseCorrelation
        cfRefImg = cfRefImg ./ (eps0 + abs(cfRefImg)) .* fhg;
    end

    if useGPU
        cfRefImg   = gpuArray(cfRefImg);
        maskOffset = gpuArray(maskOffset);
    end

    % --- batch loop ---
    nBatches = ceil(nFrames / batchSize);

    for bi = 1:nBatches
        fi = (bi-1)*batchSize + 1 : min(bi*batchSize, nFrames);

        if useGPU
            batchData = gpuArray(single(subdata(:,:,fi)));
        else
            batchData = single(subdata(:,:,fi));
        end
if numel(batchData)>532*532*5000
    useSVD=true; %force SVD on if data is too large
end
if useSVD
    % K_b=min(100,length(fi));
    K_b=max(10,ceil(length(fi)/50));
    if K_b==length(fi)
        K_b=0;
        useSVD=false;
    end
else
    K_b=0;
end

% if useSVD
    % cc0=calc_correlation(batchData,cfRefImg,fhg,eps0,lcorr,phaseCorrelation,K_b);
% else
    cc0=calc_correlation(bsxfun(@plus, maskOffset, bsxfun(@times, maskMul, batchData)),...
        cfRefImg,fhg,eps0,lcorr,phaseCorrelation,K_b);
% end
        % --- subpixel estimation ---
        if subpixel > 1
            [~,ii] = max(reshape(cc0,[],numel(fi)),[],1);
            [iy,ix] = ind2sub((2*lcorr+1)*[1 1], ii);

            dl = single(-lpad:lpad);
            if useGPU
                dl = gpuArray(dl);
                ccmat = gpuArray.zeros(numel(dl),numel(dl),numel(fi),'single');
            else
                ccmat = zeros(numel(dl),numel(dl),numel(fi),'single');
            end

            mxpt = [iy(:) ix(:)];

            for j = 1:numel(fi)
                ccmat(:,:,j) = cc0( ...
                    max(min(mxpt(j,1)+dl, size(cc0,1)),1), ...
                    max(min(mxpt(j,2)+dl, size(cc0,2)),1), j);
            end

            ccmat = reshape(ccmat,[],numel(fi));
            if kriging
                ccb = Kmat * ccmat;
                [m(ib,fi),ixg] = max(ccb,[],1);
                [ix11,ix21] = ind2sub(numel(linds)*[1 1],ixg);
                mdpt = floor(numel(linds)/2)+1;
                dv0 = ([ix11' ix21']-mdpt)/subpixel + mxpt - lcorr ...
                       - 1;
            else
                yshift = xt(1,:) * ccmat;
                xshift = xt(2,:) * ccmat;
                dv0 = [yshift' xshift'] ./ sum(ccmat,1)' + mxpt - lcorr ...
                       - 1;
                dv0 = round(dv0 * subpixel) / subpixel;
                m(ib,fi) = max(ccmat,[],1);
            end

            ds(ib,fi,:) = gather_try(dv0);

        else
            [my,iy] = max(cc0,[],1);
            [m(ib,fi),ix] = max(my,[],2);

            dv0 = [iy(:)-lcorr ix(:)-lcorr] - 1;
            ds(ib,fi,:) = gather_try(dv0);
        end
    end
end
    for b=1:size(ds,1)
        outlier=any(abs(ds(b,:,:))>=maxregshift,3);
        if any(outlier)
            if sum(~outlier)>=2
                ds(b,outlier,1)=interp1(find(~outlier),ds(b,~outlier,1),find(outlier),'linear','extrap');
                ds(b,outlier,2)=interp1(find(~outlier),ds(b,~outlier,2),find(outlier),'linear','extrap');
            else
                ds(b,:,:)=0;
            end
        end
    end
    ds(isnan(ds) | isinf(ds))=0;
    ds=double(ds);
end


function xyMask = make_xyMask(Ly, Lx, numBlocks)
%MAKE_XYMASK Generate normalized 2D Gaussian masks for blocks.
%
%   xyMask = make_xyMask(Ly, Lx, numBlocks, yB, xB, sT)
%
%   Inputs:
%       Ly        - number of rows in the image
%       Lx        - number of columns in the image
%       numBlocks - [nY, nX] number of blocks along y and x
%       yB        - vector of y-centers of blocks (length nY)
%       xB        - vector of x-centers of blocks (length nX)
%       sT        - [sigmaY, sigmaX] Gaussian std deviation per axis
%
%   Output:
%       xyMask    - [Ly*Lx x nBlocks] matrix of normalized masks
yB        = linspace(0, Ly, numBlocks(1)+1);
yB        = round((yB(1:end-1) + yB(2:end)) / 2);

xB        = linspace(0, Lx, numBlocks(2)+1);
xB        = round((xB(1:end-1) + xB(2:end)) / 2);
sT(1)        = mean(diff(yB)) * 2/3;
sT(2)        = mean(diff(xB)) * 2/3;
sT = max(10, sT);

nblocks = prod(numBlocks);
xyMask3D = zeros(Ly, Lx, nblocks, 'single');

ib = 0;

for iy = 1:numBlocks(1)
    if numBlocks(2) > 1
        for ix = 1:numBlocks(2)
            ib = ib + 1;
            
            if numBlocks(1) > 1
                gausy = exp(-((1:Ly)' - yB(iy)).^2 / (2*sT(1)^2));
            else
                gausy = ones(Ly,1,'single');
            end
            
            gausx = exp(-((1:Lx)' - xB(ix)).^2 / (2*sT(2)^2));
            
            xyMask3D(:,:,ib) = gausy * gausx';
        end
    else
        ib = ib + 1;
        if numBlocks(1) > 1
            gausy = exp(-((1:Ly)' - yB(iy)).^2 / (2*sT(1)^2));
        else
            gausy = ones(Ly,1,'single');
        end
        xyMask3D(:,:,ib) = repmat(gausy, 1, Lx);
    end
end

% Normalize masks so sum across blocks = 1
maskSum = sum(xyMask3D, 3);
maskSum(maskSum==0) = 1;  % avoid divide by zero
xyMask3D = xyMask3D ./ maskSum;

% Reshape into 2D [Ly*Lx x nBlocks]
xyMask = reshape(xyMask3D, Ly*Lx, nblocks);

% Replace NaNs (if any) with zeros
xyMask(isnan(xyMask)) = 0;

end

function [yBL, xBL] = compute_block_indices(Ly, Lx, numBlocks)
%COMPUTE_BLOCK_INDICES Compute blocks for image partitioning.
%
%   [bpix, yB, xB, yBL, xBL] = compute_block_indices(Ly, Lx, numBlocks)
%
%   Inputs:
%       Ly        - number of rows in image
%       Lx        - number of columns in image
%       numBlocks - [nY, nX] number of blocks along y and x
%
%   Outputs:
%       bpix - [height width] of each block (approx.)
%       yB   - vector of block center positions along Y
%       xB   - vector of block center positions along X
%       yBL  - cell array of row indices for each block
%       xBL  - cell array of column indices for each block

nblocks = prod(numBlocks);

% --- fraction of block relative to image ---
bfrac = 1 ./ max(2, ceil((numBlocks - 3)/3));
bfrac(numBlocks == 1) = 1;

% approximate block sizes in pixels
bpix = ceil(bfrac .* [Ly Lx]);

% --- compute block centers ---
yB = linspace(0, Ly, numBlocks(1)+1);
yB = round((yB(1:end-1) + yB(2:end)) / 2);

xB = linspace(0, Lx, numBlocks(2)+1);
xB = round((xB(1:end-1) + xB(2:end)) / 2);

% --- compute indices for each block ---
ib = 0;
yBL = cell(nblocks,1);
xBL = cell(nblocks,1);

for iy = 1:numBlocks(1)
    if numBlocks(2) > 1
        for ix = 1:numBlocks(2)
            ib = ib + 1;
            yBL{ib} = max(1, yB(iy) - floor(bpix(1)/2)) : ...
                       min(Ly, yB(iy) + floor(bpix(1)/2));
            xBL{ib} = max(1, xB(ix) - floor(bpix(2)/2)) : ...
                       min(Lx, xB(ix) + floor(bpix(2)/2));
        end
    else
        ib = ib + 1;
        yBL{ib} = max(1, yB(iy) - floor(bpix(1)/2)) : ...
                   min(Ly, yB(iy) + floor(bpix(1)/2));
        xBL{ib} = 1:Lx;
    end
end

end

function [dy,dx]=ds_to_dxy(xyMask,ds,sz)
Ly=sz(1);
Lx=sz(2);
if size(ds,1)==1 || size(ds,2)==1
    dx = ds(:,:,2);
    dy = ds(:,:,1);
    dx = dx-round(median(dx,'all'));
    dy = dy-round(median(dy,'all'));
    return;
end
try
dx = (xyMask * ds(:,:,2));
dy = (xyMask * ds(:,:,1));
dx = dx-round(median(dx,'all'));
dy = dy-round(median(dy,'all'));
catch
    batch_size=1000;
    n_rep=size(xyMask,1);
    clear dx
    dx=zeros(n_rep,size(ds,2),'single');
    clear dy
    dy=zeros(n_rep,size(ds,2),'single');

    n_batch=ceil(n_rep/batch_size);
    for i=1:n_batch
        pix=1+(i-1)*batch_size:min(i*batch_size,n_rep);
        dx(pix,:) = (xyMask(pix,:) * ds(:,:,2));
        dy(pix,:) = (xyMask(pix,:) * ds(:,:,1));
    end
end
dx = reshape(dx, Ly, Lx, []);
dy = reshape(dy, Ly, Lx, []);
end

function cc0=calc_correlation(batchData,cfRefImg,fhg,eps0,lcorr,phaseCorrelation,K_b)
if nargin<7 || isempty(K_b)
    K_b=0;
end
[ly,lx]=size(batchData,1:2);
batchData=batchData-mean(batchData,1:2);
if K_b>0
        if size(batchData,3)>1000
        [U,S,V] = bksvd(reshape(batchData,[],size(batchData,3)), K_b);
        else
            %for <1000 frames, 512x512 pixels, svdecon is faster
        [U,S,V] = svdecon(reshape(batchData,[],size(batchData,3)));
            U=U(:,1:K_b);
            S=S(1:K_b,1:K_b);
            V=V(:,1:K_b);
        end
        corrMap=fft2(reshape(U*S,ly,lx,[]));

        if phaseCorrelation
            corrMap = bsxfun(@times, ...
                corrMap ./ (eps0 + abs(corrMap)) .* fhg, cfRefImg);
        else
            corrMap = bsxfun(@times, corrMap, cfRefImg);
        end

        corrClip = real(ifft2(corrMap));
        corrClip = fftshift(fftshift(corrClip,1),2);

        % --- subpixel estimation ---
            Mode_cc0 = corrClip( ...
                floor(ly/2)+1+(-lcorr:lcorr), ...
                floor(lx/2)+1+(-lcorr:lcorr), :);

Mode_cc0 = reshape(Mode_cc0,[],K_b);
% VS = V * S;
% cc0=Mode_cc0*VS';
cc0=Mode_cc0*V';
cc0=reshape(cc0,lcorr*2+1,lcorr*2+1,[]);

% cc0=zeros(lcorr*2+1,lcorr*2+1,nf);
% corrClip=reshape(corrClip,[],size(corrClip,3));
% VS = V * S;
% for i=1:nf
% cc0_temp=reshape(corrClip*VS(i,:)',ly,lx,[]);
% cc0(:,:,i)=cc0_temp( ...
%                 floor(ly/2)+1+(-lcorr:lcorr), ...
%                 floor(lx/2)+1+(-lcorr:lcorr), :);
% end
% cc0=reshape(cc0,lcorr*2+1,lcorr*2+1,[]);

else
            corrMap = fft2(batchData);

        if phaseCorrelation
            corrMap = bsxfun(@times, ...
                corrMap ./ (eps0 + abs(corrMap)) .* fhg, cfRefImg);
        else
            corrMap = bsxfun(@times, corrMap, cfRefImg);
        end

        corrClip = real(ifft2(corrMap));
        corrClip = fftshift(fftshift(corrClip,1),2);
            cc0 = corrClip( ...
                floor(ly/2)+1+(-lcorr:lcorr), ...
                floor(lx/2)+1+(-lcorr:lcorr), :);
end


end

function refImg=gen_mimg(data,whichch,n_ch,n_Frames_template)
        div=nf/n_ch/n_Frames_template;
        first=floor(div/2)+1;
        % div=floor(div);
        frames=(round(first:div:nf/n_ch)-1)*n_ch+whichch;
        % refImg=gen_template(data(:,:,whichch:n_ch:n_ch*100),min(100,nf/n_ch));
        small_corr_stack=reg2P_standalone(data(:,:,frames),nanmean(data(:,:,frames),3));
        refImg=nanmean(small_corr_stack,3);
end