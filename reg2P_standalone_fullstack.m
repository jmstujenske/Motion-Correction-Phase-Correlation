function data=reg2P_standalone_fullstack(data,batch_size,n_iter,bidi,n_ch,whichch)
%reg2P_standalone_fullstack(data,batch_size,n_iter,bidi,n_ch,whichch)
%data - X by Y by (C*T) frame stack
%batch_size - # of frames to process at one time on one computing core (default: 500)
%n_iter - # of times to repeat motion correction (default: 1)
%bidi - whether to correct for bidirectional scanning (default: false)
%n_ch - number of channels
%whichch - which channel to motion correct based on
%
%Based on solution from Suite2p Matlab version, now made as a standable
%implementation
%https://github.com/cortex-lab/Suite2P
%
%Please cite the original authors
%
%This implementation uses the parallel processing toolbox and has two
%advancements over the original suite2p scripts: 1. sub-pixel registration,
%2. correcting for bidirectional scanning, accounting for differences in
%offset along the x-axis
%
%J.M.Stujenske, April 2023

if nargin <4 || isempty(bidi)
    bidi=false;
end
if nargin<2 || isempty(batch_size)
    batch_size=500;
end
if nargin<3 || isempty(n_iter)
    n_iter=1;
end
if ischar(data) %if filename provided instead of the data
%     [data,info]=bigread4(data);
    isfile=true;
else
    isfile=false;
end

if isfile
        info=readtifftags(data);
    Ly=info(1).ImageHeight;
    Lx=info(1).ImageWidth;
    nFrames=length(info);
    fid=fopen(data,'r');
    fseek(fid,0,'eof');
    len=ftell(fid);
    fclose(fid);
            [folder,filename,ext]=fileparts(data);
            newfile=fullfile(folder,[filename,'_motcorr',ext]);
    if len/1e9<3.99
        TiffWriter=Fast_Tiff_Write(newfile,info(1).Xresolution,0,info(1).ImageDescription);
    else
        TiffWriter=Fast_BigTiff_Write(newfile,info(1).Xresolution,0,info(1).ImageDescription);
    end
%tic;
    
    
else
    [Ly,Lx,nFrames]=size(data);
end
nreps=ceil(nFrames/batch_size);

%cut data up into cells
if ~isfile
if nFrames==nreps*batch_size
    data_cell=mat2cell(data,Ly,Lx,batch_size*ones(1,nreps));
else
    data_cell=mat2cell(data,Ly,Lx,[batch_size*ones(1,nreps-1) mod(nFrames,batch_size)]);
end

%remove empty cells
in=squeeze(cellfun(@isempty,data_cell));
data_cell(in)=[];
end


% mimg=gen_template(data(:,:,whichch:n_ch:end),min(1000,nFrames));
if isfile
    mimg=bigread4(data,1,min(1000,nFrames*n_ch));
    mimg=mean(mimg(:,:,whichch:n_ch:end),3);
else
    mimg=mean(data(:,:,whichch:n_ch:min(1000,nFrames*n_ch)),3);
end
if bidi
    [mimg]=correct_bidi_across_x(mimg,1,1);
end
if bidi && ~isfile
% [col_shift] = correct_bidirectional_offset(data,100);
% mimg=apply_col_shift(mimg,col_shift);

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
%         frames=1+batch_size*(rep-1):min(batch_size*rep,nFrames);
    % temp=reg2P_standalone(data2(:,:,frames),mimg,false);toc;
        if isfile
            for rep=1:nreps
            data_cell=bigread4(data,(rep-1)*batch_size+1,min(batch_size,nFrames-batch_size*(rep-1)));
            data_cell=correct_bidi_across_x(data_cell,n_ch,whichch);
            dreg=reg2P_standalone(data_cell,mimg,false,[32 1],n_ch,whichch);
            for a=1:size(dreg,3);TiffWriter.WriteIMG(dreg(:,:,a)');end;
            end
        else
            parfor rep=1:nreps
            dreg{rep}=reg2P_standalone(data_cell{rep},mimg,false,[32 1],n_ch,whichch);
            end
        end
    % dreg(:,:,frames)=normcorre_batch(data2(:,:,frames),options_nonrigid);
    if isfile
        close(TiffWriter);
        data=newfile;
    else
        data=cat(3,dreg{:});
    end
    if n_iter>1 && isfile
        reg2P_standalone_fullstack(newfile,batch_size,n_iter-1,bidi,n_ch,whichch);
        return;
    end
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