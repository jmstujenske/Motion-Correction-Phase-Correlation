function J = imresize_rows_only(I, newNumCols, method)
%IMRESIZE_COLS_ONLY  Resize only along columns (horizontal) while keeping rows identical.
%   J = imresize_cols_only(I, newNumCols, method)
%   - I: input image/matrix (MxN or MxNxC)
%   - newNumCols: scalar number of output columns (Nout)
%   - method: 'nearest','linear','cubic','makima','pchip',...
%
%   No vertical mixing: each output row is a 1-D interpolation of the
%   corresponding input row; rows are not combined.

if nargin < 3, method = 'linear'; end
[M, N, C] = size(I);
if length(newNumCols)==2
    newNumCols=newNumCols(2);
end
Nout = newNumCols;

% coordinate vectors for columns
x_old = (1:N).';
x_new = linspace(1, N, Nout).';

% Put data into shape [N x (M*C)] so interp1(x_old, V, xq, ...) interpolates
% along the column dimension (length(x_old) == size(V,1)).
V = permute(double(I), [2,1,3]);          % [N x M x C]
V = reshape(V, N, []);            % [N x (M*C)]

% interpolate: returns [Nout x (M*C)] — vectorized for all rows & channels
Vq = interp1(x_old, V, x_new, method, 'extrap');  % [Nout x (M*C)]

% reshape back to [M x Nout x C]
Vq = reshape(Vq, Nout, M, C);     % [Nout x M x C]
J = permute(Vq, [2,1,3]);         % [M x Nout x C]

% cast back to input class if numeric type
if ~isa(I, 'double')
    J = cast(J, class(I));
end
end
