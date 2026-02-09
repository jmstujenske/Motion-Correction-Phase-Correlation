
function [out,D_field]=align_pyramids(to_align,new_template,blocks,type,orig_type,smooth)
max_shift_init=50;
max_shift=20;
bg_filt_sz=1;
to_align_orig=to_align;
% new_template=new_template.^2;
% to_align=to_align.^2;
if nargin<6 || isempty(smooth)
    smooth=true;
end
if nargin<4 || isempty(type)
    type='affine';
end
if nargin<5 || isempty(orig_type)
    orig_type='affine';
end
if nargin<3 || isempty(blocks)
    blocks=[25 25];
end

if smooth
    sz=ceil(size(new_template)./blocks);
    pad=max(sz);
    pad=ceil(pad/2)*2;
    to_align=pad_expand_0(to_align,pad);
    new_template=pad_expand_0(new_template,pad);
    blocks=blocks+ceil(pad./sz);
end
to_align_nan=to_align;
to_align_nan(to_align==0)=NaN;
to_align_nan(bad_edge_mask(to_align))=NaN;
new_template_nan=new_template;
new_template_nan(new_template==0)=NaN;
new_template_nan(bad_edge_mask(new_template))=NaN;

new_template=fix_bad_edges(new_template,isnan(new_template_nan));
to_align=fix_bad_edges(to_align,isnan(to_align_nan));

to_align=to_align-imgaussfilt_nanignore(to_align_nan,bg_filt_sz,'Padding','symmetric');
new_template=new_template-imgaussfilt_nanignore(new_template_nan,bg_filt_sz,'Padding','symmetric');
new_template(isnan(new_template_nan))=0;
to_align(isnan(to_align_nan))=0;
to_align=max(to_align,0);
new_template=max(new_template,0);
to_align=to_align.^2;
new_template=new_template.^2;
% [ax,ay]=gradient(to_align);
% [tx,ty]=gradient(new_template);
% new_template=sqrt(tx.^2+ty.^2);
% to_align=sqrt(ax.^2+ay.^2);
to_align(isnan(to_align))=1;
new_template(isnan(new_template))=1;


    % tform = imregtform(to_align, new_template, 'affine', optimizer, metric);
    % tform = imregcorr(to_align,imref2d(size(to_align)),new_template,imref2d(size(new_template)),'transformtype','similarity','Window',true);
tform=align_pyramids_initial(to_align,new_template,orig_type);

    if strcmp(type,'demons')
        to_align=gpuArray(to_align);
        new_template=gpuArray(new_template);
    end
    [optimizer, metric] = imregconfig('multimodal');

    optimizer.InitialRadius = 0.0001;
    optimizer.Epsilon = 1.5e-4;
    optimizer.GrowthFactor = 1.01;
    optimizer.MaximumIterations = 300;
[ydim,xdim]=size(new_template,1:2);

[X,Y]=meshgrid(1:xdim,1:ydim);
if ~isempty(tform)
    D_init=get_D_field(tform,X,Y,smooth,max_shift_init);
    to_align = apply_D_field_nan(to_align,D_init);
else
    D_init=zeros([size(X),2]);
end
% spacing=20;

% smooth_func=@(x) imgaussfilt_nanignore(double(x),spacing*2/3,'Padding',0);

spacing=[(ydim+1)/(blocks(1)) (xdim+1)/(blocks(2))];
xcoords=floor(floor(spacing(2)/2)+1:spacing(2):xdim);
ycoords=floor(floor(spacing(1)/2)+1:spacing(1):ydim);
ncuts=length(xcoords)*length(ycoords);
average_d=cell(1,ncuts);
count=cell(1,ncuts);
max_spacing=max(spacing);

xs_func=@(x) max(xcoords(x)-ceil(spacing(2)*10/3)-1,1):min(xcoords(x)+ceil(spacing(2)*10/3),xdim);
ys_func=@(x) max(ycoords(x)-ceil(spacing(1)*10/3)-1,1):min(ycoords(x)+ceil(spacing(1)*10/3),ydim);
% xs_func_template=@(x) max(xcoords(x)-ceil(spacing(2)*7/3)-1,1):min(xcoords(x)+ceil(spacing(2)*7/3),xdim);
% ys_func_template=@(x) max(ycoords(x)-ceil(spacing(1)*7/3)-1,1):min(ycoords(x)+ceil(spacing(1)*7/3),ydim);

parfor cut=1:ncuts
    [y_cut,x_cut]=ind2sub([length(ycoords) length(xcoords)],cut);
xs=xs_func(x_cut);
ys=ys_func(y_cut);
% min_dim=min(length(ys),length(xs));
l=max(length(ys),length(xs));
kernel=fspecial('gaussian',[l l],ceil(max_spacing*2/3));
kernel=imresize(kernel,[length(ys) length(xs)]);
kernel=kernel./sum(kernel,'all');
% kernel=1;
smooth_func=@(x) kernel.*double(x);
sub_temp=new_template(ys,xs);
moving=to_align(ys,xs);
[average_d{cut},count{cut}]=get_D_field_patch(sub_temp,moving,smooth_func,type,smooth,max_shift);
%only want to calculate transform based on valid, not boundary conditioned
%data; data over boundary should be NaN or 0.

end

count_ag=nan(ydim,xdim,2);
average_d_ag=nan(ydim,xdim,2);
for cut=1:ncuts
[y_cut,x_cut]=ind2sub([length(ycoords) length(xcoords)],cut);
xs=xs_func(x_cut);
ys=ys_func(y_cut);
ds=average_d{cut};
cs=count{cut};
% cs(abs(ds)>10)=NaN;
% ds(abs(ds)>10)=NaN;
average_d_ag(ys,xs,:)=nansum(cat(4,average_d_ag(ys,xs,:),ds),4);
% average_d_ag(ys,xs,:,cut)=ds;
count_ag(ys,xs,:)=nansum(cat(4,count_ag(ys,xs,:),cs),4);
% count_ag(ys,xs,:,cut)=cs;
end
y_D=average_d_ag(:,:,1)./count_ag(:,:,1);
% y_D(y_D==0)=NaN;
x_D=average_d_ag(:,:,2)./count_ag(:,:,2);
% x_D(x_D==0)=NaN;

% if smooth
%     mask=bad_edge_mask_xy(x_D,y_D);
%     x_D(mask)=NaN;
%     y_D(mask)=NaN;
% end
% y_D=imgaussfilt_nanignore(y_D,.1,'Padding','symmetric');
% x_D=imgaussfilt_nanignore(x_D,.1,'Padding','symmetric');
D_field=cat(3,x_D,y_D);
if smooth
D_field=fix_bad_edges(D_field,isnan(x_D) | isnan(y_D),Y,X);
end
if smooth
bad_pix=isnan(y_D) | isnan(x_D) | sqrt(sum((D_field).^2,3))>max_shift;
D_field=fix_bad_edges(D_field,bad_pix,Y,X);
for i=1:2
    D_field(:,:,i)=medfilt2(D_field(:,:,i),[5 5]);
end
end
% if smooth
%     for i=1:2
%     D_field(:,:,i)=medfilt2(D_field(:,:,i),[5 5]);
%     end
%     D_field=imgaussfilt_nanignore(D_field,.2,'Padding','symmetric');
% end
D_field_final=zeros(size(D_field));
% D_field=D_field_final;

X_new_first = apply_D_field_nan(X,D_init);
X_new = apply_D_field_nan(X_new_first,D_field);
% X_new=X+D_init(:,:,1)+D_field(:,:,1);
% Y_new=Y+D_init(:,:,2)+D_field(:,:,2);
% X_new = imwarp(X_new,D_field);
Y_new_first = apply_D_field_nan(Y,D_init);
Y_new = apply_D_field_nan(Y_new_first,D_field);
% X_new=min(nanmax(X_new,1),max(X));
% Y_new=min(nanmax(Y_new,1),max(Y));

D_field_final(:,:,1)=X_new-X;
D_field_final(:,:,2)=Y_new-Y;

% if smooth
%     mask=bad_edge_mask_xy(D_field_final(:,:,1),D_field_final(:,:,2)) | X_new==0 | Y_new==0;
%     D_field_final(reshape(find(mask),[],1)+(0:1)*numel(mask))=NaN;
% end
% y_D=imgaussfilt_nanignore(y_D,.1,'Padding','symmetric');
% x_D=imgaussfilt_nanignore(x_D,.1,'Padding','symmetric');
if smooth
bad_pix=X_new_first==0 | Y_new_first==0 | X_new==0 | Y_new==0 | (sqrt(sum((D_field_final).^2,3))>max_shift);
D_field_final=fix_bad_edges(D_field_final,bad_pix,Y,X);
end
if smooth
    for i=1:2
    D_field_final(:,:,i)=medfilt2(D_field_final(:,:,i),[5 5]);
    end
    D_field_final=imgaussfilt_nanignore(D_field_final,.2,'Padding','symmetric');
end
D_field=D_field_final(pad+1:end-pad,pad+1:end-pad,:);
D_field(isnan(D_field))=0;
out=apply_D_field_nan(to_align_orig,D_field);
out=gather_try(out);
D_field=gather_try(D_field);
end

function D_init=get_D_field(tform,X,Y,smooth,max_shift_init)
[X_new_init,Y_new_init] = transformPointsInverse(tform,X,Y);
X_new_init=double(X_new_init);
Y_new_init=double(Y_new_init);
X_init=X_new_init-X;
X_new_init(X_new_init==0)=NaN;
Y_init=Y_new_init-Y;
X_new_init(X_new_init==0)=NaN;

D_init=cat(3,X_init,Y_init);
if smooth
bad_pix=X_new_init==0 | Y_new_init==0 | sqrt(sum(D_init.^2,3))>max_shift_init;
D_init=fix_bad_edges(D_init,bad_pix,Y,X);
end
end

function out=apply_D_field_nan(to_align,D_field)
        bad=isnan(D_field);
        D_field(bad)=0;
        out = imwarp(to_align,D_field);
        out(bad)=0;
end

function tform=align_pyramids_initial(to_align,new_template,orig_type)
    [optimizer, metric] = imregconfig('multimodal');
    switch lower(orig_type)
        case 'kaze'
            tform = align_with_KAZE(to_align,new_template);
            if isempty(tform)
                warning('KAZE features alignment failed. Using affine as backup.')
                tform = imregtform(to_align, new_template, 'affine', optimizer, metric);
            end
        case 'phase'
            tform = imregcorr(to_align,imref2d(size(to_align)),new_template,imref2d(size(new_template)),'transformtype','similarity','Window',true);
        otherwise
            try
            tform = imregtform(to_align, new_template, orig_type, optimizer, metric);
            catch
                warning('orig_type not recognized.');
                tform=[];
            end
    end
end

function [average_d,count]=get_D_field_patch(sub_temp,moving,smooth_func,type,smooth,max_shift)
    [optimizer, metric] = imregconfig('multimodal');
sub_temp(isnan(sub_temp))=0;
moving(isnan(moving))=0;
bad=sub_temp==0 | moving==0;
good_x=any(~bad,1);
bad_y=any(sub_temp(:,good_x)==0 | moving(:,good_x)==0,2);
%%
xs_sub=good_x;
ys_sub=~bad_y;
if sum(ys_sub)>16 && sum(xs_sub)>16 %%for pyramid level 3, need at least 16 in each dim
    switch type
        case 'demons'
            D_sub_demons=imregdemons(moving(ys_sub,xs_sub), sub_temp(ys_sub,xs_sub),'DisplayWaitbar',false,'AccumulatedFieldSmoothing', 5);
            D_sub=zeros([size(moving), 2]);
            D_sub(ys_sub,xs_sub,:)=D_sub_demons;
        otherwise
            tform_sub = imregtform(moving(ys_sub,xs_sub), sub_temp(ys_sub,xs_sub),type, optimizer, metric);
                % tform_sub = imregcorr(moving(ys_sub,xs_sub),imref2d(size(moving(ys_sub,xs_sub))),sub_temp(ys_sub,xs_sub),imref2d(size(sub_temp(ys_sub,xs_sub))),'transformtype','similarity','Window',true);

            % tform_sub = imregtform(moving, sub_temp,type, optimizer, metric);
            [X,Y]=meshgrid(1:size(moving,2),1:size(moving,1));
            D_sub=get_D_field(tform_sub,X,Y,smooth,max_shift);
    end
else
    switch type
        case 'demons'
            D_sub=imregdemons(moving, sub_temp,'DisplayWaitbar',false,'AccumulatedFieldSmoothing', 5);
        otherwise
            tform_sub = imregtform(moving, sub_temp,type, optimizer, metric);
            % tform_sub = imregtform(moving, sub_temp,type, optimizer, metric);
            [X,Y]=meshgrid(1:size(moving,2),1:size(moving,1));
            D_sub=get_D_field(tform_sub,X,Y,smooth,max_shift);
    end
    % D_sub=nan([size(moving,1:2) 2]); %%suppress patches for which there is insufficient valid data
end
D_sub(abs(D_sub)>max_shift)=NaN;

Xdisp=D_sub(:,:,1);
Ydisp=D_sub(:,:,2);
average_d=cat(3,smooth_func(Ydisp),smooth_func(Xdisp));
count=cat(3,smooth_func(Ydisp./Ydisp),smooth_func(Xdisp./Xdisp));
end
