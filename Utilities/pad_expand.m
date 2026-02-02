function data=pad_expand(data,pad)
sz=size(data,1:3);
if pad>0
% data_temp=data;
data=[data(pad+1:-1:2,:,:);data;data(end:-1:end-pad+1,:,:)];
data=[data(:,pad+1:-1:2,:) data data(:,end:-1:end-pad+1,:)];
end
% data=imgaussfilt(data,pad,'Padding','symmetric');
% data(pad+1:end-pad,pad+1:end-pad,:)=data_temp;