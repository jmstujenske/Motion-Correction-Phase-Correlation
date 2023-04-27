function data=pad_expand(data,pad)
data_temp=data;
data=[data(pad+1:-1:2,:,:);data;data(end:-1:end-pad+1,:,:)];
data=[data(:,pad+1:-1:2,:) data data(:,end:-1:end-pad+1,:)];
data=imgaussfilt(data,floor(pad/2));
data(pad+1:end-pad,pad+1:end-pad,:)=data_temp;