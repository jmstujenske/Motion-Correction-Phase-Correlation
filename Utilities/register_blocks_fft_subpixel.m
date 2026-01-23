function [ds,m] = register_blocks_fft_subpixel( ...
    data, mimg,whichch,n_ch,numBlocks, ...
    maxregshift, subpixel, ...
    smoothSigma, maskSlope, ...
    phaseCorrelation, kriging, ...
    useGPU, eps0,batchSize,useSVD,bidi_comp)

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
if nargin<15 || isempty(bidi_comp)
    bidi_comp=false;
end
if nargin<14 || isempty(useSVD)
    useSVD=false;
end
if nargin<13 || isempty(batchSize)
    batchSize=1000;
end
if nargin<3 || isempty(whichch)
    whichch=1;
end
if nargin<4 || isempty(n_ch)
    n_ch=1;
end
[Ly,Lx,nf] = size(data);
nFrames=nf/n_ch;
if nargin<4 || isempty(numBlocks)
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

    refImg = single(mimgB{ib});
    refImg = zscore(refImg.^2,[],1:2);

    if nblocks>1
    subdata = data(yBL{ib}, xBL{ib}, whichch:n_ch:end);
    end
    % if nblocks>1
    % subdata = max(subdata-quantile(subdata,.05,3),0);
    % else
    % data = max(data-quantile(data,.05,3),0);
    % end
    % subdata = imgaussfilt(subdata, 0.5, 'Padding', 'symmetric');

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

    if bidi_comp
    refImg=[refImg(1:2:end,:);refImg(2:2:end,:)];
    end
    % FFT of reference
    cfRefImg = conj(fftn(refImg));

    if useGPU
        cfRefImg   = gpuArray(cfRefImg);
        maskOffset = gpuArray(maskOffset);
    end

    % --- batch loop ---
    nBatches = ceil(nFrames / batchSize);

    for bi = 1:nBatches
        fi = (bi-1)*batchSize + 1 : min(bi*batchSize, nFrames);

        if nblocks>1
        if useGPU
            batchData = gpuArray(subdata(:,:,fi));
        else
            batchData = subdata(:,:,fi);
        end
        if ib==nblocks && bi == nBatches
            clear subdata data;
        end
        else
        if useGPU
            batchData = gpuArray(data(:,:,fi));
        else
            batchData = data(:,:,fi);
        end
        end
if numel(batchData)>532*532*5000
    useSVD=true; %force SVD on if data is too large
end
if useSVD
    % K_b=min(100,length(fi));
    K_b=min(max(5,ceil(length(fi)/20)),50);
    if K_b==length(fi)
        K_b=0;
        useSVD=false;
    end
else
    K_b=0;
end
    if phaseCorrelation && ~useSVD
        cfRefImg = cfRefImg ./ (eps0 + abs(cfRefImg)) .* fhg;
    end
% if useSVD
    % cc0=calc_correlation(batchData,cfRefImg,fhg,eps0,lcorr,phaseCorrelation,K_b);
% else
    [cc0,lcorr2]=calc_correlation(bsxfun(@plus, maskOffset, bsxfun(@times, maskMul, single(batchData).^2)),...
        cfRefImg,fhg,eps0,lcorr,phaseCorrelation,K_b,bidi_comp);
% end
        % --- subpixel estimation ---
        if subpixel > 1
            [~,ii] = max(reshape(cc0,[],numel(fi)),[],1);
            [iy,ix] = ind2sub(size(cc0,1:2), ii);

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
                dv0 = ([ix11' ix21']-mdpt)/subpixel + mxpt - lcorr2 ...
                       - 1;
            else
                yshift = xt(1,:) * ccmat;
                xshift = xt(2,:) * ccmat;
                dv0 = [yshift' xshift'] ./ sum(ccmat,1)' + mxpt - lcorr2 ...
                       - 1;
                dv0 = round(dv0 * subpixel) / subpixel;
                m(ib,fi) = max(ccmat,[],1);
            end

            ds(ib,fi,:) = gather_try(dv0);

        else
            [my,iy] = max(cc0,[],1);
            [m(ib,fi),ix] = max(my,[],2);

            dv0 = [iy(:)-lcorr2(1) ix(:)-lcorr2(2)] - 1;
            ds(ib,fi,:) = gather_try(dv0);
        end
    end
end
    if bidi_comp
    ds(:,:,1)=ds(:,:,1)*2;
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