function dreg=apply_reg2P_shifts(data,varargin)
%dreg=apply_reg2P_shifts(data,xyMask,ds)
%dreg=apply_reg2p_shifts(data,shifts)
%
if iscell(varargin{1})
    xyMask=varargin{1}{1};
    if length(varargin{1})>2
        ds_rigid=varargin{1}{2};
        ds=varargin{1}{3};
    else
        ds=varargin{1}{2};
%         if size(ds,1)==1
%             ds_rigid=ds;
%             ds=[];
%         end
    end
    if nargin>2
        pad=varargin{2};
    else
        pad=[];
    end
else
    xyMask=varargin{1};
    ds=varargin{2};
    if nargin>3
        pad=varargin{3};
    else
        pad=[];
    end
end
[Ly,Lx,NT]=size(data);
if isempty(pad)
    %numel_diff = (x+pad)*(y+pad)-x*y = x*y +x*pad+y*pad+pad^2 - x *y = x*pad+y*pad+pad^2 = pad*(x+y)+pad^2
    %apply quadratic formula
    %0 = pad^2+(x+y)*pad-numel_diff = (-(x+y)+sqrt((x+y)^2+4*numel_diff))/2;
    numel_diff=size(xyMask,1)-numel(data(:,:,1));
    pad=(-Lx-Ly+sqrt((Lx+Ly)^2+4*numel_diff))/2;
    pad=pad/2;
end
Ly=Ly+pad*2;
Lx=Lx+pad*2;
data=pad_expand(data,pad);

if exist('ds_rigid','var')
    ds_rigid=reshape(ds_rigid,1,[],2);
    nshifts=size(ds_rigid,2);
    n_ch=NT/nshifts;
for i = 1:NT
    frame_num=ceil(i/n_ch);
    Im = data(:,:,i);
    data(:,:,i)=imwarp(Im,repmat(ds_rigid(1,frame_num,:),size(data,1:2)));
end
end
if ~isempty(ds)
dx = (xyMask * ds(:,:,2));
dy = (xyMask * ds(:,:,1));
class_data=class(data);
dx = reshape(dx, Ly, Lx, []);
dy = reshape(dy, Ly, Lx, []);
nshifts=size(dx,3);
n_ch=NT/nshifts;
idy = repmat([1:Ly]', 1, Lx);
idx = repmat([1:Lx],  Ly, 1);

dreg = zeros(size(data), class_data);
        dx_i=dx(:,:,1);
        dy_i=dy(:,:,1);
for i = 1:NT
    frame_num=ceil(i/n_ch);
    Im = data(:,:,i);
    if nshifts>1
        dx_i=dx(:,:,frame_num);
        dy_i=dy(:,:,frame_num);
    end
    dreg(:,:,i)=imwarp(cast(Im,class_data),cat(3,dx_i,dy_i));
end
else
    dreg=data;
end
dreg=dreg(pad+1:end-pad,pad+1:end-pad,:);
