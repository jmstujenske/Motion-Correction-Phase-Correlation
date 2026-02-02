function Vx=fix_bad_edges(Vx,bad_pix,Y,X)
if nargin<3 || isempty(Y)
[X,Y]=meshgrid(1:size(Vx,2),1:size(Vx,1));
end
n=size(Vx,3);
for i=1:n
    V_slice=Vx(:,:,i);
Fx = scatteredInterpolant(X(~bad_pix),Y(~bad_pix), V_slice(~bad_pix), 'linear', 'nearest');

V_slice(bad_pix) = Fx(X(bad_pix), Y(bad_pix));
Vx(:,:,i)=V_slice;
end
end