function out=apply_nonrigid_dx_parallel(X,dx,dy,n_ch,subpixel,pad,use_par)
if nargin<5
    subpixel=false;
end
if nargin<6 || isempty(pad)
    pad=0;
end
if nargin<7 || isempty(use_par)
    use_par=true;
end
if use_par
    use_par=inf;
    if isempty(gcp('nocreate'))
        parpool;
    end
else
    use_par=0;
end
orig_x=1:size(X,2)+pad*2;
orig_y=1:size(X,1)+pad*2;
nf=size(X,3)/n_ch;
dx=double(dx);
dy=double(dy);
if size(dx,3)==1
    dx=reshape(dx,1,1,[]);
    dy=reshape(dy,1,1,[]);
end
changes=find(squeeze(any(abs(dx)>.5 | abs(dy)>.5,[1 2])));
fs=1:nf;
if isempty(changes)
    out=X;
    return;
end
fs=fs(changes);
[x_coords,y_coords]=meshgrid(orig_x,orig_y);
allframes=sort(reshape(reshape(fs*n_ch,[],1)+(0:n_ch-1),[],1))';
nSubFrames=length(allframes);
batchSize=50;
nBatches=ceil(nSubFrames/batchSize);
out=X;

for ib=1:nBatches
    sub_frames=allframes(1+(ib-1)*batchSize:min(ib*batchSize,nSubFrames));
        out_sub=zeros(length(orig_y),length(orig_x),length(sub_frames));

    if use_par>0
        X_broadcast=parallel.pool.Constant(X(:,:,sub_frames));
        dx_broadcast=parallel.pool.Constant(dx(:,:,ceil(sub_frames/n_ch)));
        dy_broadcast=parallel.pool.Constant(dy(:,:,ceil(sub_frames/n_ch)));
    else
        X_broadcast.Value=X(:,:,sub_frames);
        dx_broadcast.Value=dx(:,:,ceil(sub_frames/n_ch));
        dy_broadcast.Value=dy(:,:,ceil(sub_frames/n_ch));
    end
    counts=1:length(sub_frames);

    parfor (c=counts,use_par)
            Im=X_broadcast.Value(:,:,c);
        Im = pad_expand(Im,pad);
        x_coords_new=x_coords+(dx_broadcast.Value(:,:,c));
        y_coords_new=y_coords+(dy_broadcast.Value(:,:,c));
        x_coords_new=min(max(round(x_coords_new),1),size(Im,2));
        y_coords_new=min(max(round(y_coords_new),1),size(Im,1));

        % out(:,:,(i-1)*n_ch+j)=imwarp(Im,cat(3,dx_i,dy_i),'linear');
        if ~subpixel
            ins=(x_coords_new-1)*size(Im,1)+y_coords_new;
            out_sub(:,:,c)=Im(double(ins));
        else
            out_sub(:,:,c)= images.internal.interp2d(Im,x_coords_new,y_coords_new,'linear',0); %less overhead than imwarp
        end
    end
    if pad>0
        out_sub=out_sub(pad+1:end-pad,pad+1:end-pad,:);
    end
    out(:,:,sub_frames)=out_sub;
end
if use_par>0
    clear X_broadcast;
end

end