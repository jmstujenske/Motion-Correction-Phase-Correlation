function fixed=apply_bidi_correction(data,dx,lowmemory)
if nargin<3 || isempty(lowmemory)
    lowmemory=false;
end
[Ly,Lx,nFrames]=size(data);
Ly=ceil(Ly/2);
upsample=10;

Lx=Lx*upsample;
%%subpixel-registration
%
%setup the conversion
idy = repmat((1:Ly)', 1, Lx);
idx = repmat((1:Lx),  Ly, 1);
dx=round(dx);

% apply offsets to indexing
DX = dx + idx;

xyInvalid = DX<0 | DX>Lx-1;

DX(xyInvalid) = 0;

% ind = max(idy + DX * Ly,1);

%upsample the other channels
fixed=data;
if lowmemory
    batch_size=16;
else
    batch_size=Ly;
end
n_batches=ceil(Ly/batch_size);
    for batch=1:n_batches
        lines=1+(batch-1)*batch_size*2:2:min(batch_size*2*batch,Ly*2);
        this_batch_size=length(lines);
        ins=(1:this_batch_size)+(batch-1)*batch_size;
        cur_ind = max(idy(ins,:)-(batch-1)*batch_size + DX(ins,:) * this_batch_size,1);
        upsized=imresize(data(lines,:,:),[this_batch_size Lx]);
%     OLD implementation: loop through every frame and apply the same correction to each
%     loop through every frame and apply the same correction to each
%         for i=1:nFrames
%             subpix_fix=upsized(cur_ind+(this_batch_size*Lx*(i-1)));
%             %downsample back to original sampling
%             fixed(lines,:,n_ch*(i-1)+rep)=imresize(subpix_fix,[Ly Lx/upsample]);
%         end
        inds=(cur_ind+this_batch_size*Lx*reshape((0:nFrames-1),1,1,[]));
        subpix_fix=upsized(inds); %directly indexing the matrix is 30% faster for large matrices
        fixed(lines,:,:)=imresize(subpix_fix,[this_batch_size Lx/upsample]);
    end
end