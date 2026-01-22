function [dy,dx]=ds_to_dxy(xyMask,ds,sz)
Ly=sz(1);
Lx=sz(2);
if size(ds,1)==1
    dx = ds(:,:,2);
    dy = ds(:,:,1);
    dx = dx-round(median(dx,'all'));
    dy = dy-round(median(dy,'all'));
    dx = reshape(dx,1,1,[]);
    dy = reshape(dy,1,1,[]);
    return;
end
try
dx = (xyMask * ds(:,:,2));
dy = (xyMask * ds(:,:,1));
dx = dx-round(median(dx,'all'));
dy = dy-round(median(dy,'all'));
catch
    batch_size=1000;
    n_rep=size(xyMask,1);
    clear dx
    dx=zeros(n_rep,size(ds,2),'single');
    clear dy
    dy=zeros(n_rep,size(ds,2),'single');

    n_batch=ceil(n_rep/batch_size);
    for i=1:n_batch
        pix=1+(i-1)*batch_size:min(i*batch_size,n_rep);
        dx(pix,:) = (xyMask(pix,:) * ds(:,:,2));
        dy(pix,:) = (xyMask(pix,:) * ds(:,:,1));
    end
end
dx = reshape(dx, Ly, Lx, []);
dy = reshape(dy, Ly, Lx, []);
end