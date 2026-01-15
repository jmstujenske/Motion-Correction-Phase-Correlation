function data=reg2P_standalone_fullstack(data,batch_size,bidi,n_ch,whichch,numBlocks,use_par,save_path)
%reg2P_standalone_fullstack(data,batch_size,n_iter,bidi,n_ch,whichch)
%data - X by Y by (C*T) frame stack OR filename (tif)
%batch_size - # of frames to process at one time on one computing core (default: 500)
%bidi - whether to correct for bidirectional scanning (default: false)
%n_ch - number of channels
%whichch - which channel to motion correct based on
%numBlocks - blocks for motion correct (default: [32 1]);
%if you specify a second row, then it does a slow drift correction, e.g.
%[32 1; 32 32];
%use_par - which or not to use parallel computing (default: false);
%save_path - where to save new file for file input (default: location of
%orig file; if specify 0 then will output matrix)
%
%Based on solution from Suite2p Matlab version, now made as a standable
%implementation
%https://github.com/cortex-lab/Suite2P
%
%Please also cite the original authors.
%
%This implementation uses the parallel processing toolbox for matrix input
%and has three
%advancements over the original suite2p scripts: 1. sub-pixel registration,
%2. correcting for bidirectional scanning, accounting for differences in
%offset along the x-axis, 3. Optional parallel computing (memory demanding)
%
%J.M.Stujenske, April 2023
%
if ischar(data) %if filename provided instead of the data
    %     [data,info]=bigread4(data);
    isfile=true;
else
    isfile=false;
end
if nargin < 8 || isempty(save_path)
    if isfile
    folder=fileparts(data);
    else
        folder=dir;
    end
    save_path=folder;
end
if nargin < 7 || isempty(use_par) || ~use_par
    use_par=0;
end
if use_par
    use_par=inf;
end
if nargin < 6 || isempty(numBlocks)
    numBlocks=[32 1];
end
if nargin < 4 || isempty(n_ch)
    n_ch=1;
end
if nargin <5 || isempty(whichch)
    whichch=1;
end
if nargin <3 || isempty(bidi)
    bidi=false;
end
if nargin<2 || isempty(batch_size)
    batch_size=500;
end
if numel(numBlocks)>2
    [data]=reg2P_standalone_fullstack_slowdrift(data,batch_size,bidi,n_ch,whichch,true,numBlocks,use_par);
    return;
end
if isfile
    [folder,filename,ext]=fileparts(data);

    fid=fopen(data,'r');
    fseek(fid,0,'eof');
    len=ftell(fid);
    fclose(fid);
    try
        [m,~,info]=memory_map_tiff(data,'matrix',1,true);
        Ly=info(1).ImageHeight;
        Lx=info(1).ImageWidth;
        nFrames=length(info);
        memmap=true;
    catch
        memmap=false;
        info=readtifftags(data);
        Ly=info(1).ImageHeight;
        Lx=info(1).ImageWidth;
        nFrames=length(info);
    end
    newfile=fullfile(save_path,[filename,'_motcorr',ext]);
    if isfield(info,'ImageDescription')
        desc=info(1).ImageDescription;
    else
        desc=[];
    end

    if len/1e9<3.99
        TiffWriter=Fast_Tiff_Write(newfile,info(1).Xresolution,0,desc);
    else
        TiffWriter=Fast_BigTiff_Write(newfile,info(1).Xresolution,0,desc);
    end
    %tic;


else
    [Ly,Lx,nFrames]=size(data);
end
nreps=ceil(nFrames/batch_size);

%cut data up into cells
if isfile
    if ~memmap
        template_frames=bigread4(data,1,min(500*n_ch,nFrames));
        template_frames=template_frames(:,:,whichch:n_ch:end);
        [dreg,shifts]=reg2P_standalone(template_frames,median(template_frames,3),[],[1 1],1,1,200);
        mimg=median(dreg,3);
    else
        template_frames=m.Data.allchans(:,:,whichch:n_ch:min(500*n_ch,nFrames));
        [dreg,shifts]=reg2P_standalone(template_frames,median(template_frames,3),[],[1 1],1,1,200);
        mimg=median(dreg,3);
        mimg=mimg';
    end
else
    template_frames=data(:,:,whichch:n_ch:min(500*n_ch,nFrames*n_ch));
    [dreg,shifts]=reg2P_standalone(template_frames,median(template_frames,3),[],[1 1],1,1,200);
    mimg=median(dreg,3);
    if nFrames==nreps*batch_size
        data_cell=mat2cell(data,Ly,Lx,batch_size*ones(1,nreps));
    else
        data_cell=mat2cell(data,Ly,Lx,[batch_size*ones(1,nreps-1) mod(nFrames,batch_size)]);
    end
    clear data
    %remove empty cells
    in=squeeze(cellfun(@isempty,data_cell));
    data_cell(in)=[];
end
% mimg=gen_template(data(:,:,whichch:n_ch:end),min(1000,nFrames));

if bidi
    [mimg,dx_bidi]=correct_bidi_across_x(mimg,1,1,false,true,1);
end
if bidi && ~isfile
    % [col_shift] = correct_bidirectional_offset(data,100);
    % mimg=apply_col_shift(mimg,col_shift);
    for rep=1:nreps
        %         data_cell{rep}=apply_col_shift(data_cell{rep},col_shift);
        % data_cell{rep}=correct_bidi_across_x(data_cell{rep},n_ch,whichch);
        data_cell{rep}=apply_bidi_correction_direct(data_cell{rep},dx_bidi,n_ch,true,1);
    end
end
% dreg=zeros(Ly,Lx,nFrames,'single');
% dims=size(data2);
% options_rigid = NoRMCorreSetParms('d1',size(data2,1),'d2',size(data2,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200);
% options_nonrigid = NoRMCorreSetParms('d1',dims(1),'d2',dims(2),'grid_size',[64,64],'mot_uf',4,'bin_width',800,'max_shift',[2 2],'max_dev',[50 50],'us_fac',50,'init_batch',100,'shifts_method','cubic');
dreg=cell(nreps,1);

%         frames=1+batch_size*(rep-1):min(batch_size*rep,nFrames);
% temp=reg2P_standalone(data2(:,:,frames),mimg,false);toc;
if isfile
    if use_par
        myPool=gcp('nocreate');
        if isempty(myPool)
            parpool();
            myPool=gcp('nocreate');
        end
        nw=myPool.NumWorkers;
        for outrep=1:nw:nreps
            dreg=cell(nw,1);
            data_cell=cell(nw,1);
            for rep=1:(min(nw,nreps-outrep+1))
                if memmap
                    data_cell{rep}=m.Data.allchans(:,:,(rep+outrep-2)*batch_size+1:(rep+outrep-2)*batch_size+min(batch_size,nFrames-batch_size*(rep+outrep-2)));
%                     data_cell{rep}=reshape(data_cell{rep},Lx,Ly,[]);
                    data_cell{rep}=permute(data_cell{rep},[2 1 3]);
                else
                    data_cell{rep}=bigread4(data,(rep+outrep-2)*batch_size+1,min(batch_size,nFrames-batch_size*(rep+outrep-2)));
                end
                if bidi
                    data_cell{rep}=apply_bidi_correction_direct(data_cell{rep},dx_bidi,n_ch,true,1);
                end
            end
            try
            parfor (rep=1:(min(nw,nreps-outrep+1)),use_par)
                % data_cell{rep}=correct_bidi_across_x(data_cell{rep},n_ch,whichch,true); %low memory mode

                dreg{rep}=reg2P_standalone_twostep(data_cell{rep},mimg,false,numBlocks,n_ch,whichch);
                % if size(numBlocks,1)>1
                % [~,shifts]=reg2P_standalone(mean(dreg{rep},3),mimg,false,numBlocks(2,:),n_ch,whichch,10);
                % dreg{rep}=apply_reg2P_shifts(dreg{rep},shifts);
                % end
            end
            catch
                use_par=false;
             for rep=1:(min(nw,nreps-outrep+1))
                dreg{rep}=reg2P_standalone_twostep(data_cell{rep},mimg,false,numBlocks,n_ch,whichch);
            end
            end
            for a=1:nw
                if ~isempty(dreg{a})
                    for f=1:size(dreg{a},3)
                        TiffWriter.WriteIMG(dreg{a}(:,:,f)');
                    end
                end
            end
        end
    else
        for rep=1:nreps
            if memmap
                data_cell=m.Data.allchans(:,:,(rep-1)*batch_size+1:(rep-1)*batch_size+min(batch_size,nFrames-batch_size*(rep-1)));
                data_cell=permute(data_cell,[2 1 3]);
            else
                data_cell=bigread4(data,(rep-1)*batch_size+1,min(batch_size,nFrames-batch_size*(rep-1)));
            end
            % data_cell=correct_bidi_across_x(data_cell,n_ch,whichch,true); %low memory mode
            data_cell=apply_bidi_correction_direct(data_cell,dx_bidi,n_ch,true,1);
            dreg=reg2P_standalone_twostep(data_cell,mimg,false,numBlocks,n_ch,whichch);
            for a=1:size(dreg,3);TiffWriter.WriteIMG(dreg(:,:,a)');end;
        end
    end
else
    try
    parfor (rep=1:nreps, use_par)
        dreg{rep}=reg2P_standalone_twostep(data_cell{rep},mimg,false,numBlocks,n_ch,whichch);
    end
    catch
        use_par=false;
        for rep=1:nreps
            dreg{rep}=reg2P_standalone_twostep(data_cell{rep},mimg,false,numBlocks,n_ch,whichch);
        end
    end
end
% dreg(:,:,frames)=normcorre_batch(data2(:,:,frames),options_nonrigid);
if isfile
    close(TiffWriter);
    if memmap
        clear m
    end
    if save_path~=0
        data=newfile;
    else
        data=bigread4(newfile);
    end

elseif ~isfile
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
