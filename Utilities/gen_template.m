function [mimg,frames]=gen_template(data,n_Frames_template,n_ch,whichch)
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
    mimg=nanmean(single(data(:,:,frames)),3);
else
    mimg=zeros(Ly,Lx,'single');
    for a=frames
%         keyboard
        frame=bigread4(data,a,1,info);
        if a==frames(1)
            classdata=class(frame);
        end
        mimg=imadd(mimg,single(frame)/n_Frames_template);
    end
    mimg=cast(mimg,classdata);
end