function out=apply_nonrigid_dx(X,dx,dy,n_ch,subpixel,pad)
if nargin<5
    subpixel=false;
end
if nargin<6
    pad=0;
end
orig_x=1:size(X,2)+pad*2;
orig_y=1:size(X,1)+pad*2;
nf=size(X,3)/n_ch;
changes=find(squeeze(any(abs(dx)>.5 | abs(dy)>.5,[1 3])));
fs=1:nf;
% fs=fs(changes);
out=zeros(size(dx,1),size(dx,2),size(X,3),class(X));
[x_coords,y_coords]=meshgrid(orig_x,orig_y);

if ~subpixel
for i=fs
    fs_allch=(i-1)*n_ch+(1:n_ch);
    X_sub=X(:,:,fs_allch);
    X_sub=pad_expand(X_sub,pad);
    nPix=prod(size(X_sub,1:2));
    if ismember(i,changes)
    x_coords_new=x_coords+double(dx(:,:,i));
    y_coords_new=y_coords+double(dy(:,:,i));
    x_coords_new=min(max(round(x_coords_new),1),size(X_sub,2));
    y_coords_new=min(max(round(y_coords_new),1),size(X_sub,1));
    ins_orig=(x_coords-1)*size(X_sub,1)+y_coords+reshape(fs_allch-1,1,1,[])*nPix;
    ins=(x_coords_new-1)*size(X_sub,1)+y_coords_new+reshape(0:n_ch-1,1,1,[])*nPix;
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

            % dx_i=dx(:,:,i);
            % dy_i=dy(:,:,i);
                            x_coords_new=x_coords+double(dx(:,:,i));
y_coords_new=y_coords+double(dy(:,:,i));
x_coords_new=min(max(round(x_coords_new),1),size(X,2));
y_coords_new=min(max(round(y_coords_new),1),size(X,1));

            % out(:,:,(i-1)*n_ch+j)=imwarp(Im,cat(3,dx_i,dy_i),'linear');

            out(:,:,(i-1)*n_ch+j)= images.internal.interp2d(Im,x_coords_new,y_coords_new,'linear',0); %less overhead than imwarp
                else
                    out(:,:,(i-1)*n_ch+j)=Im;
                end
        end
    end
end
out=out(pad+1:end-pad,pad+1:end-pad,:);
end