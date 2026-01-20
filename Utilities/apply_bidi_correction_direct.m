function fixed=apply_bidi_correction_direct(data,dx,n_ch,lowmemory,subpixel)
if nargin<3 || isempty(n_ch)
    n_ch=1;
end
if nargin<4 || isempty(lowmemory)
    lowmemory=false;
end
if nargin<5 || isempty(subpixel)
    % subpixel=10;
    subpixel=length(dx)/Lx;
else
    dx=imresize(dx,[1 Lx*subpixel]);
end
[Ly,Lx,nFrames]=size(data);
Ly=ceil(Ly/2);
idy = repmat((1:Ly)', 1, Lx);
idx = repmat((0:Lx-1),  Ly, 1)*subpixel+floor(subpixel/2)+1;
% dx=round(dx);
dx=dx*subpixel;
if length(dx)~=Lx
dx=imresize(dx,[1 Lx],'Antialiasing',false);
end
% apply offsets to indexing
DX = round(dx+idx);
DX(DX<0) = 0;
DX(DX>Lx*subpixel-1)=Lx*subpixel;
% ind = max(idy + DX * Ly,1);

%upsample the other channels
fixed=data;
if lowmemory
    batch_size=16;
else
    batch_size=Ly;
end
n_batches=ceil(Ly/batch_size);

idxfirst=floor(subpixel/2)+1;
frame_batch_size=5000;
n_frame_batch=ceil(nFrames/frame_batch_size);
for fr=1:n_frame_batch
for rep=1:n_ch
    for batch=1:n_batches
        lines=1+(batch-1)*batch_size*2:2:min(batch_size*2*batch,Ly*2);
        frames=rep+(fr-1)*frame_batch_size*n_ch:n_ch:min(frame_batch_size*n_ch*fr,nFrames);
        nframes_batch=length(frames);
        this_batch_size=length(lines);
        ins=(1:this_batch_size)+(batch-1)*batch_size;
        % cur_ind = max(idy(ins,:)-(batch-1)*batch_size + DX(ins,:) * this_batch_size,1);
        if subpixel~=1
            upsized=imresize_rows_only(double(data(lines,:,frames)),[this_batch_size Lx*subpixel],'linear');
        else
            upsized=double(data(lines,:,frames));
        end

                cur_ind = max(idy(ins,:)-(batch-1)*batch_size+ (DX(ins,:)-1)* this_batch_size,1);

        inds=cur_ind+reshape(((0:nframes_batch-1)),1,1,[]).*this_batch_size*(Lx*subpixel);
        fixed(lines,:,frames)=upsized(inds);
        end
    end
end
end