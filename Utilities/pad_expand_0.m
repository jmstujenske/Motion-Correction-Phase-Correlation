function data=pad_expand_0(data,pad)
sz=size(data,1:3);
if pad>0
% data_temp=data;
data=[zeros(pad,sz(2),sz(3));data;zeros(pad,sz(2),sz(3))];
data=[zeros(sz(1)+pad*2,pad,sz(3)) data zeros(sz(1)+pad*2,pad,sz(3))];
end
% data=imgaussfilt(data,pad,'Padding','symmetric');
% data(pad+1:end-pad,pad+1:end-pad,:)=data_temp;