function mimg=gen_template(data,n_Frames_template)
[Ly,Lx,nFrames]=size(data);
if nargin<2 || isempty(n_Frames_template)
n_Frames_template=1000;
end
div=nFrames/n_Frames_template;
first=floor(div/2)+1;
% div=floor(div);
mimg=nanmean(data(:,:,round(first:div:end)),3);