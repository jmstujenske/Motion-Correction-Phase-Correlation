function X=apply_rigid_dx(X,dx,dy,n_ch,subpixel,bidi_dx)
if nargin<5 || isempty(subpixel)
    subpixel=false;
end
if nargin<6 || isempty(bidi_dx)
    bidi_dx=0;
elseif ~isscalar(bidi_dx)
    bidi_dx=imresize_rows_only(bidi_dx,[1 size(X,2)]);
end
orig_x=1:size(X,2);
orig_y=1:size(X,1);
nf=size(X,3)/n_ch;
dx=double(dx);
dy=double(dy);
changes=(abs(dx)>0 | abs(dy)>0);
fs=1:nf;
fs=fs(changes);
[x_coords,y_coords]=meshgrid(orig_x,orig_y);
if ~subpixel
nPixels=prod(size(X,1:2));
for i=fs
    fs_allch=(i-1)*n_ch+(1:n_ch);
    x_coords_new=x_coords+double(dx(i));
    x_coords_new(1:2:end,:)=x_coords_new(1:2:end,:)+bidi_dx;
y_coords_new=y_coords+double(dy(i));
x_coords_new=min(max(round(x_coords_new),1),size(X,2));
y_coords_new=min(max(round(y_coords_new),1),size(X,1));
ins_orig=(x_coords-1)*size(X,1)+y_coords+reshape(fs_allch-1,1,1,[])*nPixels;
ins=(x_coords_new-1)*size(X,1)+y_coords_new+reshape(fs_allch-1,1,1,[])*nPixels;
X(double(ins_orig))=X(double(ins));
end

else
    %             dx=repmat(dx,size(X,1:2));
    %         dy=repmat(dy,size(X,1:2));
    % for i = fs
    %     for j=1:n_ch
    %         Im = X(:,:,(i-1)*n_ch+j);
    % 
    %         out(:,:,(i-1)*n_ch+j)=imwarp(Im,cat(3,dx(:,:,i),dy(:,:,i)),'linear');
    %     end
    % end
    for i = fs
        for j=1:n_ch
            Im = X(:,:,(i-1)*n_ch+j);
            x_coords_new=x_coords+double(dx(i));
            y_coords_new=y_coords+double(dy(i));
            x_coords_new=min(max(round(x_coords_new),1),size(X,2));
            y_coords_new=min(max(round(y_coords_new),1),size(X,1));
            X(:,:,(i-1)*n_ch+j)= images.internal.interp2d(Im,x_coords_new,y_coords_new,'linear',0);%less overhead than imwarp  
        end
    end
end
end