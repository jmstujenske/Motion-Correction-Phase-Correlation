function [U,S,V] = bksvd(A, k, iter, bsize, center)
%--------------------------------------------------------------------------
%Source: https://github.com/cpmusco/bksvd
%MIT License (see m file for full license)
%
% Randomized block Krylov iteration for truncated singular value decomposition
% Computes approximate top singular vectors and corresponding values
% Described in Musco, Musco, 2015 (http://arxiv.org/abs/1504.05477)
%
% usage : 
%
%  input:
%  * A : matrix to decompose
%  * k : number of singular vectors to compute, default = 6
%  * iter : number of iterations, default = 3
%  * bsize : block size, must be >= k, default = k
%  * center : set to true if A's rows should be mean centered before the
%  singular value decomposition (e.g. when performing principal component 
%  analysis), default = false
%
%
%  output:
%  k singular vector/value pairs. 
%  * U : a matrix whose columns are approximate top left singular vectors for A
%  * S : a diagonal matrix whose entries are A's approximate top singular values
%  * V : a matrix whose columns are approximate top right singular vectors for A
%
%  U*S*V' is a near optimal rank-k approximation for A
%--------------------------------------------------------------------------

% % MIT License
% % 
% % Copyright (c) 2017 Christopher Musco, Cameron Musco
% % 
% % Permission is hereby granted, free of charge, to any person obtaining a copy
% % of this software and associated documentation files (the "Software"), to deal
% % in the Software without restriction, including without limitation the rights
% % to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% % copies of the Software, and to permit persons to whom the Software is
% % furnished to do so, subject to the following conditions:
% % 
% % The above copyright notice and this permission notice shall be included in all
% % copies or substantial portions of the Software.
% % 
% % THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% % IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% % FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% % AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% % LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% % OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% % SOFTWARE.

% Check input arguments and set defaults.
if nargin > 5
    error('bksvd:TooManyInputs','requires at most 5 input arguments');
end
if nargin < 1
    error('bksvd:TooFewInputs','requires at least 1 input argument');
end
if nargin < 2
    k = 6;
end
k = min(k,min(size(A)));

if nargin < 3 || isempty(iter)
    iter = 3;
end
if nargin < 4 || isempty(bsize)
    bsize = k;
end
if nargin < 5 || isempty(center)
    center = false;
end
if(k < 1 || iter < 1 || bsize < k)
    error('bksvd:BadInput','one or more inputs outside required range');
end

% Calculate row mean if rows should be centered.
u = zeros(1,size(A,2));
if(center)
    u = mean(A);
end
l = ones(size(A,1),1);


% We want to iterate on the smaller dimension of A.
[n, ind] = min(size(A));
tpose = false;
if(ind == 1) 
    tpose = true;
    l = u'; u = ones(1,size(A,1));
    A = A';
end

% Allocate space for Krylov subspace.
K = zeros(size(A,2),bsize*iter);
% Random block initialization.
block = randn(size(A,2),bsize);
[block,R] = qr(block,0);

% Preallocate space for temporary products.
T = zeros(size(A,2),bsize);

% Construct and orthonormalize Krlov Subspace. 
% Orthogonalize at each step using economy size QR decomposition.
for i=1:iter
    T = A*block - l*(u*block);
    block = A'*T - u'*(l'*T);
    [block,R] = qr(block,0);
    K(:,(i-1)*bsize+1:i*bsize) = block;
end
[Q,R] = qr(K,0);

% Rayleigh-Ritz postprocessing with economy size dense SVD.
T = A*Q - l*(u*Q);

[Ut,St,Vt] = svdecon(T);
S = St(1:k,1:k);
if(~tpose)
    U = Ut(:,1:k);
    V = Q*Vt(:,1:k);
else
    V = Ut(:,1:k);
    U = Q*Vt(:,1:k);
end

end


