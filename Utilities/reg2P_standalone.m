function [dreg,shifts,m]=reg2P_standalone(data,mimg,varargin)
%kriging,numBlocks,n_ch,whichch,maxregshift,fs,quick,use_subpixel_reg,bidi_comp,bidi_correct)
%
%Two input options:
%dreg=reg2P_standalone(data,mimg,kriging,numBlocks,n_ch,whichch,maxregshift,fs,quick,use_subpixel_reg,bidi_comp,bidi_correct)
%
%OR
%
%dreg=reg2P_standalone(data,mimg,options)
%where options is a structure with fields corresponding to above variables
%
%data - X by Y by (C*T) frame stack
%mimg - template image (default: 1000 frame average)
%kriging - whether to use kriging or not (default: true)
%numBlocks - how many blocks to break image into (default: [1 1])
%n_ch - how many channels in data (default: 1)
%whichch - which channel to motion correct based on (default: 1)
%maxregshift - maximum shift (default: 30)
%fs - frame rate (default: 30)
%quick - run a quick rigid correction and then only apply nonrigid
%correction to frames with dips in cross-correlation (default: true)
%use_subpixel_reg - whether to use subpixel registration when applying
%motion correction (default: true)
%bidi_comp - whether to compensate for differences introduced by
%bidirectional scanning when calculating cross-correlation (only necessary
%for more hardware configurations; default: true)
%
%Expansion of solution from Suite2p Matlab version, now made as a
%standalone
%implementation
%https://github.com/cortex-lab/Suite2P
%
%Stujenske, JM Jan 2026
%
%
if isstruct(varargin{1})
    options_in=varargin{1};
    options=parse_options(options_in);
else
    options=parse_inputs(varargin);
end
[kriging,numBlocks,n_ch,whichch,maxregshift,fs,quick,use_subpixel_reg,...
    bidi_comp,bidi_correct,trim,subpixel,useGPU,phaseCorrelation,maskSlope,...
    smoothSigma,eps0]=eval_options(options);

% if subpixel is still inf, threshold it for new method
subpixel = min(10, subpixel);

nf=size(data,3);
if nf/n_ch<=1+fs*4
quick=false;
end
if ~isempty(mimg) && bidi_correct
    [~,dx_bidi]=correct_bidi_across_x(mimg,1,1,false,false,false);
end

if quick
    if nargin<2 || isempty(mimg)
        % refImg=gen_mimg(data,whichch,n_ch,n_Frames_template);
        % refImg=mean(data(:,:,whichch:n_ch:end),3);
        refImg=gen_template(data,min(1000,nf/n_ch),n_ch,whichch);
        if bidi_correct
            if ~bidi_comp
            [refImg,dx_bidi]=correct_bidi_across_x(refImg,1,1,false,false,false);
            refImg=refImg-movmin(refImg,trim,2)+mean(movmin(refImg,trim,2),'all');
            % dx_bidi=nanmedian(dx_bidi);
            data=apply_bidi_correction_direct(data,dx_bidi,n_ch,true,1);
            else
            [~,dx_bidi]=correct_bidi_across_x(refImg,1,1,false,false,false);
            end
        else
            dx_bidi=0;
        end
    else
        refImg=mimg;
    end
    % fft_sz=5000;
        % [y,x,m]=cross_corr_memory(fft_sz);

% [ds,m] = register_blocks_fft_subpixel( ...
%     data(:,:,whichch:n_ch:end), refImg,[], ...
%     maxregshift, subpixel, ...
%     smoothSigma, maskSlope, ...
%     phaseCorrelation, kriging, ...
%     useGPU, eps0);
[ds_rigid,m] = register_blocks_fft_subpixel( ...
    data(trim+1:end-trim,trim+1:end-trim,:), refImg(trim+1:end-trim,trim+1:end-trim),whichch,n_ch,[], ...
    maxregshift, subpixel, ...
    smoothSigma, maskSlope*10, ...
    phaseCorrelation, kriging, ...
    useGPU, eps0,50000,true,bidi_comp);
% ds_rigid=movmean(movmedian(ds_rigid,3,2),3,2);
    % xyMask_rigid = make_xyMask(size(data,1), size(data,2), [1 1]);
    % [y,x]=ds_to_dxy(xyMask_rigid,ds,size(data,1:2));
    [y,x]=ds_to_dxy([],ds_rigid,size(data,1:2));
        if bidi_comp && bidi_correct
            data=apply_rigid_dx(data,x,y,n_ch,use_subpixel_reg,dx_bidi);
        else
            data=apply_rigid_dx(data,x,y,n_ch,use_subpixel_reg);
        end
            % x=ds(:,:,2);
            % y=ds(:,:,1);

        % x=(x-median(x));
        % y=(y-median(y));
                    m_smooth=sgolayfilt(double(m),3,floor(fs/2)*2+1);
            m_smooth=m_smooth-movmax(sgolayfilt(m_smooth,3,floor(fs/2)*2+1),[fs*3 0]);
            th=multithresh(m_smooth)*.8;
% ds_rigid=cat(3,y(:)',x(:)');
% if prod(numBlocks)==1
%     dreg_full=data;
%     shifts={[],ds_rigid};
%     return;
% end
movements=m_smooth<th;% | sqrt(x.^2+y.^2)>2; %only need to apply non-rigid to
% cases where the maximum correlation after rigid correction dips below
% expected, or large displacements, to be conservative
movements_pad=find(conv(movements,ones(1,max(ceil(fs/10),3)),'same')>0);
clear movements
movements_pad=(movements_pad(:)-1)*n_ch+(1:n_ch)*n_ch;
movements_pad=movements_pad(:);
if isempty(movements_pad)
    shifts=[];
    return;
end
data_sub=data(:,:,movements_pad);
NT=length(movements_pad);
else
    NT=size(data,3);
    movements_pad=1:size(data,3);
    m=[];
end
    [Ly,Lx] = size(data,1:2);


class_data=class(data);

nFrames=NT/n_ch;
if nargin<2 || isempty(mimg)
    if quick
        channel_idx=whichch:n_ch:nf;
        mimg=gen_template(data(:,:,channel_idx(setdiff(1:nf,movements_pad))),min(1000,nFrames));
    else
        mimg=gen_template(data,min(1000,nFrames),whichch,n_ch);
        if bidi_correct
            [~,dx_bidi]=correct_bidi_across_x(mimg,1,1,false,false,false);
        else
            dx_bidi=0;
        end
        data=apply_bidi_correction_direct(data,dx_bidi,n_ch,true);
    end
    mimg=single(mimg);

else
    mimg=single(mimg);
end
mimg(mimg<0)=0;

nblocks = prod(numBlocks);
xyMask = make_xyMask(Ly, Lx, numBlocks);   
if nblocks>1 || ~quick
if quick
    ds = register_blocks_fft_subpixel( ...
    data_sub, mimg,whichch,n_ch,numBlocks, ...
    maxregshift, subpixel, ...
    smoothSigma, maskSlope, ...
    phaseCorrelation, kriging, ...
    useGPU, eps0,5000,false,bidi_comp);
else
    ds = register_blocks_fft_subpixel( ...
    data, mimg,whichch,n_ch,numBlocks, ...
    maxregshift, subpixel, ...
    smoothSigma, maskSlope, ...
    phaseCorrelation, kriging, ...
    useGPU, eps0,5000,false,bidi_comp);
end
        [dy,dx]=ds_to_dxy(xyMask,ds,[Ly Lx]);
else
dy=[];dx=[];ds=[];
end

        if nargout>1
            if quick
                shifts={xyMask,ds_rigid,ds};
            else
                shifts={xyMask,ds};
            end
% end
        end
clear xyMask ds
if ~isempty(dx) && ~all(dx==0 & dy==0,'all')
    if quick
        if size(dx,1)==1
            dreg=apply_rigid_dx(data_sub,dx,dy,n_ch,use_subpixel_reg);
        else
            dreg=apply_nonrigid_dx(data_sub,dx,dy,n_ch,use_subpixel_reg);
        end
    else
        if size(dx,1)==1
            dreg=apply_rigid_dx(data,dx,dy,n_ch,use_subpixel_reg);
        else
            dreg=apply_nonrigid_dx(data,dx,dy,n_ch,use_subpixel_reg);
        end
    end
else
    if quick
        dreg=data_sub;
    else
        dreg=data;
    end
end

% dreg=apply_nonrigid_dx(data_movtrunc_allch,dx,dy,n_ch,false,pad);

if quick
    dreg_full=zeros(size(data),class_data);
    full_mov=reshape(movements_pad+(0:n_ch-1),[],1);
    dreg_full(:,:,setdiff(1:nf,full_mov))=data(:,:,setdiff(1:nf,full_mov));
    dreg_full(:,:,full_mov)=dreg;
    dreg=dreg_full;clear dreg_full;
else
    dreg=cast(dreg,class_data);
end
end

function defaults=get_defaults()
defaults={'kriging','numBlocks','n_ch','whichch','maxregshift','fs','quick','use_subpixel_reg','bidi_comp','bidi_correct','trim','subpixel','useGPU','phaseCorrelation','maskSlope','smoothSigma','eps0';...
            true,   [1 1],       1,    1,       30,           27,    true,   true,              false,       false,        30,    10,        false,   true,              5,          1.15,         1e-10};
end

function options_out=parse_options(options)
defaults=get_defaults();
options_out=struct();
nf=size(defaults,2);
fnames=fieldnames(options);
for i=1:nf
    field_i=defaults{1,i};
    detect_match=find(strcmp(field_i,fnames));
    if ~isempty(detect_match)
        options_out.(field_i)=options.(fnames{detect_match});
    end
    if isempty(detect_match) | isempty(options_out.(field_i))
        options_out.(field_i)=defaults{2,i};
    end
end

end

function options=parse_inputs(inputs)
defaults=get_defaults();
n_inputs=length(inputs);
nf=size(defaults,2);
for i=1:nf
    field_i=defaults{1,i};
    if i<=n_inputs
        options.(field_i)=inputs{i};
        if isempty(options.(field_i))
            options.(field_i)=defaults{2,1};
        end
    else
        options.(field_i)=defaults{2,i};
    end
end
end

function [kriging,numBlocks,n_ch,whichch,maxregshift,fs,quick,...
    use_subpixel_reg,bidi_comp,bidi_correct,trim,subpixel,useGPU,...
    phaseCorrelation,maskSlope,smoothSigma,eps0]=eval_options(options)
kriging=options.kriging;
numBlocks=options.numBlocks;
n_ch=options.n_ch;
whichch=options.whichch;
maxregshift=options.maxregshift;
fs=options.fs;
quick=options.quick;
use_subpixel_reg=options.use_subpixel_reg;
bidi_comp=options.bidi_comp;
bidi_correct=options.bidi_correct;
trim=options.trim;
subpixel=options.subpixel;
useGPU=options.useGPU;
phaseCorrelation=options.phaseCorrelation;
maskSlope=options.maskSlope;
smoothSigma=options.smoothSigma;
eps0=options.eps0;
end
% function out=fft2_memory(X,block_size)
% if nargin<2 || isempty(block_size)
%     block_size=5000;
% end
% nf=size(X,3);
% n_block = ceil(nf/block_size);
% out=cell(1,n_block);
% for i=1:n_block
%     fs=1+(i-1)*block_size:min(i*block_size,nf);
%     newval=fft2(X(:,:,fs));
% out{i}=newval; %for some reason, trying to populate a matrix directly causes memory issues, so store data into cell array and then concatenate at the end
% end
% out=cat(3,out{:});
% end

        % function [y,x,m]=cross_corr_memory(block_size)
        %     if nargin<1 || isempty(block_size)
        %         block_size=5000;
        %     end
        %     [Ysize,Xsize,nf]=size(data_input,1:3);
        %     n_block = ceil(nf/block_size);
        %             cfRefImg = conj(fft2(single(refImg)));
        %     x=zeros(1,nf/n_ch);
        %     y=zeros(1,nf/n_ch);
        %     m=zeros(1,nf/n_ch);
        %     for i=1:n_block
        %         frames=whichch+(i-1)*block_size*n_ch:n_ch:min(i*block_size*n_ch,nf);
        %         corrMap=fft2(single(data_input(:,:,frames)));
        % 
        %         if phaseCorrelation
        %             absRef   = abs(cfRefImg);
        %             hgx = exp(-(((0:Xsize-1) - fix(Xsize/2))/smoothSigma).^2);
        %             hgy = exp(-(((0:Ysize-1) - fix(Ysize/2))/smoothSigma).^2);
        %             hg = hgy'*hgx;
        %             fhg = real(fftn(ifftshift(single(hg/sum(hg(:))))));
        %             cfRefImg = cfRefImg./(eps0 + absRef) .* fhg;
        %         end
        %         if phaseCorrelation
        %             corrMap = bsxfun(@times, corrMap./(eps0 + abs(corrMap)).*fhg, cfRefImg);
        %         else
        %             corrMap = bsxfun(@times, corrMap, cfRefImg);
        %         end
        % 
        %         % corrMap=corrMap.*cfRefImg;
        %         corrClip = real(ifft2(corrMap));
        %         corrClip = fftshift(fftshift(corrClip, 1), 2);
        %         Y_mid=floor(Ysize/2)+1;
        %         X_mid=floor(Xsize/2)+1;
        %         corrClip=corrClip(Y_mid+(-maxregshift:maxregshift),X_mid+(-maxregshift:maxregshift),:);
        %         corrClip=imresize(corrClip,subpixel);
        %         [m((frames-whichch)/n_ch+1),vals]=max(reshape(corrClip,[],size(corrClip,3)));
        %         [y_t,x_t]=ind2sub(size(corrClip,1:2),vals);
        %         y((frames-whichch)/n_ch+1)=(y_t-maxregshift*subpixel-1)/subpixel;
        %         x((frames-whichch)/n_ch+1)=(x_t-maxregshift*subpixel-1)/subpixel;
        %     end
        % end

