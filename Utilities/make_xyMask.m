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