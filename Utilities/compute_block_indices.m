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
trim=0;
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
            yBL{ib} = max(trim+1, yB(iy) - floor(bpix(1)/2)) : ...
                       min(Ly-trim, yB(iy) + floor(bpix(1)/2));
            xBL{ib} = max(trim+1, xB(ix) - floor(bpix(2)/2)) : ...
                       min(Lx-trim, xB(ix) + floor(bpix(2)/2));
        end
    else
        ib = ib + 1;
        yBL{ib} = max(trim+1, yB(iy) - floor(bpix(1)/2)) : ...
                   min(Ly-trim, yB(iy) + floor(bpix(1)/2));
        xBL{ib} = trim+1:Lx-trim;
    end
end

end