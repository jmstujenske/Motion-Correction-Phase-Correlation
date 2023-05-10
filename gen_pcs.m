function [U,lowdim]=gen_pcs(data,n_Frames_template,n_ch,whichch)
if nargin<4 || isempty(n_ch)
    n_ch=1;
end
if nargin<5 || isempty(whichch)
    whichch=1;
end
if ischar(data)
    isfile=true;
else
    isfile=false;
end
if isfile
        info=readtifftags(data);
            Ly=info(1).ImageHeight;
    Lx=info(1).ImageWidth;
    nFrames=length(info);
else
    [~,~,nFrames]=size(data);
end
nFrames=floor(nFrames/n_ch);
if nargin<2 || isempty(n_Frames_template)
    n_Frames_template=1000;
end
div=nFrames/n_Frames_template;
first=floor(div/2)+1;
% div=floor(div);
frames=(round(first:div:nFrames)-1)*n_ch+whichch;
if ~isfile
    mimg=data(:,:,frames);
else
    frame=bigread4(data,1,1,info);
    classdata=class(frame);

    mimg=zeros(Ly,Lx,n_Frames_template,classdata);
    count=0;
    for a=frames
        count=count+1;
%         keyboard
        frame=bigread4(data,a,1,info);
        mimg(:,:,count)=frame;
    end
end
[Ly,Lx]=size(mimg,1:2);
[U,S,V]=svdecon(single(reshape(mimg,[],n_Frames_template)));
if nargout>1
    nPC=10;
bord_pix=30;
PCs=U(:,1:nPC)';
U=reshape(U,Ly,Lx,[]);
U=U(1+bord_pix:end-bord_pix,1+bord_pix:end-bord_pix,:);
PCs=reshape(U(:,:,1:nPC),(Ly-bord_pix*2)*(Lx-bord_pix*2),[])';
n_rep=ceil(div);
lowdim=zeros(nPC,nFrames);
tic;
for batch=1:n_rep
    frames=1+(batch-1)*n_Frames_template:min(batch*n_Frames_template,nFrames);
    if ~isfile
        batch_frame=data(:,:,frames);
    else
        batch_frame=bigread4(data,frames(1),length(frames),info);
    end
    batch_frame=batch_frame(1+bord_pix:end-bord_pix,1+bord_pix:end-bord_pix,:);
    batch_frame=single(batch_frame);
    lowdim(:,frames)=PCs*reshape(batch_frame,[],length(frames));
end
toc;
end