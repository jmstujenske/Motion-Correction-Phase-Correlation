function [data]=reg2P_standalone_fullstack_slowdrift(data,batch_size,bidi,n_ch,whichch,memmap,numBlocks,use_par,save_path,fs)
%data_out=reg2P_standalone_fullstack_slowdrift(data,batch_size,n_iter,bidi,n_ch,whichch)
%Motion corrects for line scanning and then drift corrects with a more fine
%nonrigid correction every batch_size/2
%The difference of this from reg2P_standalone_fullstack with two numBlocks
%specified is that the slowdrift is smoothed to removed
%
%data - X by Y by (C*T) frame stack OR filename (tif)
%batch_size - # of frames to process at one time on one computing core (default: 500)
%bidi - whether to correct for bidirectional scanning (default: false)
%n_ch - number of channels
%whichch - which channel to motion correct based on
%memmap - try to memory map tiff and edit directly (default: false)
%numBlocks - two steps; default is [32 1;32 32];
%
%OUTPUT:
%data_out - either the corrected frame stack OR new filename (tif)
%
%Based on solution from Suite2p Matlab version, now made as a standable
%implementation
%https://github.com/cortex-lab/Suite2P
%
%Please also cite the original authors.
%
%This implementation uses the parallel processing toolbox and has two
%advancements over the original suite2p scripts: 1. sub-pixel registration,
%2. correcting for bidirectional scanning, accounting for differences in
%offset along the x-axis
%
%J.M.Stujenske, April 2023
%
if nargin < 4 || isempty(n_ch)
    n_ch=1;
end
if nargin <5 || isempty(whichch)
    whichch=1;
end
if nargin <6 || isempty(memmap)
    memmap=false;
end
if nargin<7 || isempty(numBlocks)
    numBlocks=[32 1;20 20];
end
if nargin<8 || isempty(use_par)
    use_par=0;
end
if use_par
    use_par=inf;
end

if numel(numBlocks)>=4
    slowdriftmode=true;
else
    slowdriftmode=false;
end
if nargin <3 || isempty(bidi)
    bidi=false;
end
if nargin<2 || isempty(batch_size)
    batch_size=500;
end
% if numel(numBlocks)<4
%     % error('For one step correction, use reg2P_standalone_fullstack');
%     data=reg2P_standalone_fullstack(data,batch_size,bidi,n_ch,whichch,numBlocks,true);
%     return;
% end
if ischar(data) %if filename provided instead of the data
    %[data,info]=bigread4(data);
    isfile=true;
else
    isfile=false;
end
if nargin < 9 || isempty(save_path)
    if isfile
    [folder,filename,ext]=fileparts(data);
    else
        folder=dir;
    end
    save_path=folder;
end
if nargin<10 || isempty(fs)
    fs=30;
end
if isfile
    if memmap
        try
            % [m,~,info]=memory_map_tiff(data,'matrix',n_ch,false);
            info=readtifftags(data);
            Ly=info(1).ImageHeight;
            Lx=info(1).ImageWidth;
            nFrames=length(info);
            [m]=memory_map_tiff(data,'matrix',1,false,nFrames,nFrames);
            memmap=true;
        catch
            warning('Memory mapping failed.')
            memmap=false;
        end
    end
    newfile=fullfile(save_path,[filename,'_sdmotcorr',ext]);

    % if ~memmap
        info=readtifftags(data);
        Ly=info(1).ImageHeight;
        Lx=info(1).ImageWidth;
        nFrames=length(info);
        fid=fopen(data,'r');
        fseek(fid,0,'eof');
        len=ftell(fid);
        fclose(fid);
        % [folder,filename,ext]=fileparts(data);
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
    % end
    
else
    memmap=false;
    [Ly,Lx,nFrames]=size(data);
end
batch_size=min(batch_size,nFrames);
nreps=ceil(nFrames/batch_size);

%cut data up into cells

if isfile
    if ~memmap
        template_frames=bigread4(data,1,min(500*n_ch,nFrames));
        template_frames=template_frames(:,:,whichch:n_ch:end);
    else
        template_frames=m.Data.allchans(:,:,whichch:n_ch:min(500*n_ch,nFrames));
        template_frames=permute(template_frames,[2 1 3]);
    end
else
    template_frames=data(:,:,whichch:n_ch:min(500*n_ch,nFrames));
end
        if bidi
            [template_frames,bidi_dx]=correct_bidi_across_x(template_frames,1,1);
        end
    [dreg]=reg2P_standalone(template_frames,median(template_frames,3),[],[1 1],1,1,200);
    mimg=median(dreg,3);

if ~isfile
    if nFrames==nreps*batch_size
        data_cell=mat2cell(data,Ly,Lx,batch_size*ones(1,nreps));
    else
        data_cell=mat2cell(data,Ly,Lx,[batch_size*ones(1,nreps-1) mod(nFrames,batch_size)]);
    end
    %remove empty cells
    in=squeeze(cellfun(@isempty,data_cell));
    data_cell(in)=[];
    clear data %data is now redundant
end

if bidi && ~isfile
    parfor (rep=1:nreps,use_par)
        data_cell{rep}=apply_bidi_correction_direct(data_cell{rep},bidi_dx,n_ch,true,length(bidi_dx)/Lx);        
    end
end

if ~isfile
    dreg_full=cell(nreps,1);
end

shifts={[],[];[],[]};
for rep=1:nreps
    disp(['Processing Chunk ',num2str(rep),' of ',num2str(nreps)])
    %load a block of data
    if isfile
        if ~memmap
            data_cell=bigread4(data,(rep-1)*batch_size+1,min(batch_size,nFrames-batch_size*(rep-1)));
            data_cell=apply_bidi_correction_direct(data_cell{rep},bidi_dx,n_ch,true,length(bidi_dx)/Lx);
        else
            frames=(rep-1)*batch_size+1:(rep-1)*batch_size+min(batch_size,nFrames-batch_size*(rep-1));
            data_cell=apply_bidi_correction_direct(permute(m.Data.allchans(:,:,frames),[2 1 3]),bidi_dx,n_ch,true,length(bidi_dx)/Lx);
        end
    end
    
    %correct for scanning line shifts
    if isfile
        % [dreg]=reg2P_standalone_twostep(data_cell,mimg,false,numBlocks(1,:),n_ch,whichch);
        dreg=reg2P_standalone(data_cell,mimg,false,numBlocks(1,:),n_ch,whichch,[],fs,true);
        clear data_cell
    else
        % [dreg]=reg2P_standalone_twostep(data_cell{rep},mimg,false,numBlocks(1,:),n_ch,whichch);
        dreg=reg2P_standalone(data_cell{rep},mimg,false,numBlocks(1,:),n_ch,whichch,[],fs,true);
        data_cell{rep}=[];%delete data that we are done with
    end

    %     shifts{1}=cat(1,shifts{1},shift_temp_first);
    nF_dreg=size(dreg,3);
    %correct for bigger distortions, which only really occur gradually over
    %time, by breaking into halves, and then we will smooth over three half-batch weighted running
    %average
    
    %do not recalculate shifts if there is too little data at the end of
    %the file; otherwise calculate shifts
    % if rep~=nreps || size(dreg,3)>=batch_size/2
    if size(dreg,3)>=batch_size/2
        if slowdriftmode
          [~,shift_temp]=reg2P_standalone(mean(dreg(:,:,whichch:n_ch:floor(nF_dreg/2)),3),mimg,false,numBlocks(2,:),1,1,10,[],false);
          [~,shift_temp2]=reg2P_standalone(mean(dreg(:,:,floor(nF_dreg/2)+whichch:n_ch:end),3),mimg,false,numBlocks(2,:),1,1,10,[],false);
            shifts=cat(1,shifts(end-1:end,:),shift_temp,shift_temp2);
        end
    end
    
    %apply shifts
    if rep~=1
        %if not the first repetition, we have half of the data from the
        %prior batch that still needs to be shifted
        
        %smooth over three batches
        if slowdriftmode
                shifts_temp=cat(4,shifts{1:3,2});
                shifts_temp=sum(shifts_temp.*reshape([1 2 1]/4,[1 1 1 3]),4);
%         shifts_temp=shifts{2,2};
        %apply shifts
        dreg_prior=apply_reg2P_shifts(dreg_prior,{shifts{2,1},shifts_temp});
        end
        
        if isfile
            %write data to tiff file
            % if ~memmap
                for a=1:size(dreg_prior,3)
                    TiffWriter.WriteIMG(dreg_prior(:,:,a)');
                end
            % else
            %     m.Data.allchans(:,:,frames(floor(batch_size/2)+1:end)-batch_size)=permute(dreg_prior,[2 1 3]);
            % end
        else
            %prior batch data needs shifting, and we only saved the first
            %half
            dreg_full{rep-1}=cat(3,dreg_full{rep-1},dreg_prior);
        end
        
        %smooth the first half of the current batch now
        if slowdriftmode  
                shifts_temp=cat(4,shifts{2:4,2});
%         shifts_temp=shifts{3,2};
                shifts_temp=sum(shifts_temp.*reshape([1 2 1]/4,[1 1 1 3]),4);
        dreg(:,:,1:floor(nF_dreg/2))=apply_reg2P_shifts(dreg(:,:,1:floor(nF_dreg/2)),{shifts{3,1},shifts_temp});
        end
    else
        %if the first batch, then there is no prior to deal with
        if slowdriftmode
                shifts_temp=cat(4,shifts{:,2});
%         shifts_temp=shifts{3,2};
        
                shifts_temp=sum(shifts_temp.*reshape([2 1]/3,[1 1 1 2]),4);
        dreg(:,:,1:floor(nF_dreg/2))=apply_reg2P_shifts(dreg(:,:,1:floor(nF_dreg/2)),{shifts{3,1},shifts_temp});
        end
    end
    
    %update template in an ongoing manner
    if rep==1
        mimg=mean(dreg(:,:,whichch:n_ch:end),3);
    else
    mimg=nanmean(dreg_prior,3)/2+nanmean(dreg(:,:,1:floor(nF_dreg/2)),3)/2;
    end
    if rep~=nreps
        %if not the last batch, then cut data in half
        dreg_prior=dreg(:,:,floor(nF_dreg/2)+1:end);
    else
        %if it is the last batch, then we correct the last batch just based on
        %the last two corrections
        if slowdriftmode
                shifts_temp=cat(4,shifts{end-1:end,2});
                shifts_temp=mean(shifts_temp.*reshape([1 2]/3*2,[1 1 1 2]),4);
%         shifts_temp=shifts{end,2};
        dreg(:,:,floor(nF_dreg/2)+1:end)=apply_reg2P_shifts(dreg(:,:,floor(nF_dreg/2)+1:end),{shifts{end,1},shifts_temp});
        end
    end
    if isfile
        %write the first half of the data to tif
        % if ~memmap
            for a=1:floor(nF_dreg/2)
                TiffWriter.WriteIMG(dreg(:,:,a)');
            end
            %if it is the last batch, then write the rest of the data
            if rep==nreps
                for a=floor(nF_dreg/2)+1:nF_dreg
                    TiffWriter.WriteIMG(dreg(:,:,a)');
                end
            end
        % else
        %     m.Data.allchans(:,:,frames(1:floor(nF_dreg/2)))=permute(dreg(:,:,1:floor(nF_dreg/2)),[2 1 3]);
        %     if rep==nreps
        %         m.Data.allchans(:,:,frames(floor(nF_dreg/2)+1:end))=permute(dreg(:,:,floor(nF_dreg/2)+1:end),[2 1 3]);
        %     end
        % end
    else
        if rep~=nreps
            %if not the last batch, just save half the data
            dreg_full{rep}=dreg(:,:,1:floor(nF_dreg/2));
        else
            %if last batch, save all of the data
            dreg_full{rep}=dreg;
        end
    end
end
% dreg(:,:,frames)=normcorre_batch(data2(:,:,frames),options_nonrigid);

if isfile
    close(TiffWriter);
    data=newfile;
else
    %combine all of the data together
    data=cat(3,dreg_full{:});
end
end
% 
% function out=apply_col_shift(in,col_shift)
% [Ly,Lx,nFrames]=size(in);
% d_s=imresize(in,[Ly Lx*10]);
% d_s(2:2:end,max(1,1-col_shift*10):min(Lx*10,Lx*10-col_shift*10),:)=d_s(2:2:end,max(1,1+col_shift*10):min(Lx*10,Lx*10+col_shift*10),:);
% if col_shift<0
%     d_s(2:2:end,1:1-col_shift*10,:)=repmat(d_s(2:2:end,1,:),[1,1-col_shift*10 1]);
% else
%     d_s(2:2:end,Lx*10-col_shift*10:end,:)=repmat(d_s(2:2:end,Lx*10-col_shift*10,:),[1,1+col_shift*10 1]);
% end
% 
% out=imresize(d_s,[Ly Lx]);
% end
