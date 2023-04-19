function data=reg2P_standalone_fullstack(data,blocksize,n_iter,bidi,n_ch,whichch)
if nargin <4 || isempty(bidi)
    bidi=false;
end
if nargin<2 || isempty(blocksize)
    blocksize=500;
end
if nargin<3 || isempty(n_iter)
    n_iter=1;
end

[Ly,Lx,nFrames]=size(data);
nreps=ceil(nFrames/blocksize);
if nFrames==nreps*blocksize
    data_cell=mat2cell(data,Ly,Lx,blocksize*ones(1,nreps));
else
    data_cell=mat2cell(data,Ly,Lx,[blocksize*ones(1,nreps-1) mod(nFrames,blocksize)]);
end
in=squeeze(cellfun(@isempty,data_cell));
data_cell(in)=[];
mimg=gen_template(data(:,:,whichch:n_ch:end),min(1000,nFrames));
if bidi
% [col_shift] = correct_bidirectional_offset(data,100);
% mimg=apply_col_shift(mimg,col_shift);
[mimg]=correct_bidi_across_x(mimg,1,1);
    parfor rep=1:nreps
%         data_cell{rep}=apply_col_shift(data_cell{rep},col_shift);
          data_cell{rep}=correct_bidi_across_x(data_cell{rep},n_ch,whichch);
    end
end
% dreg=zeros(Ly,Lx,nFrames,'single');
% dims=size(data2);
% options_rigid = NoRMCorreSetParms('d1',size(data2,1),'d2',size(data2,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200);
% options_nonrigid = NoRMCorreSetParms('d1',dims(1),'d2',dims(2),'grid_size',[64,64],'mot_uf',4,'bin_width',800,'max_shift',[2 2],'max_dev',[50 50],'us_fac',50,'init_batch',100,'shifts_method','cubic');
dreg=cell(nreps,1);
for iter=1:n_iter
    parfor rep=1:nreps
%         frames=1+blocksize*(rep-1):min(blocksize*rep,nFrames);
    % temp=reg2P_standalone(data2(:,:,frames),mimg,false);toc;
        dreg{rep}=reg2P_standalone(data_cell{rep},mimg,false,[32 1],n_ch,whichch);
    % dreg(:,:,frames)=normcorre_batch(data2(:,:,frames),options_nonrigid);
    end
    data=cat(3,dreg{:});
end
end

function out=apply_col_shift(in,col_shift)
[Ly,Lx,nFrames]=size(in);
d_s=imresize(in,[Ly Lx*10]);
d_s(2:2:end,max(1,1-col_shift*10):min(Lx*10,Lx*10-col_shift*10),:)=d_s(2:2:end,max(1,1+col_shift*10):min(Lx*10,Lx*10+col_shift*10),:);
if col_shift<0
d_s(2:2:end,1:1-col_shift*10,:)=repmat(d_s(2:2:end,1,:),[1,1-col_shift*10 1]);
else
d_s(2:2:end,Lx*10-col_shift*10:end,:)=repmat(d_s(2:2:end,Lx*10-col_shift*10,:),[1,1+col_shift*10 1]);
end

out=imresize(d_s,[Ly Lx]);
end