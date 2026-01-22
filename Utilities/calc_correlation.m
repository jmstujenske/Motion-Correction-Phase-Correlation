function cc0=calc_correlation(batchData,cfRefImg,fhg,eps0,lcorr,phaseCorrelation,K_b,bidi_comp)
if nargin<7 || isempty(K_b)
    K_b=0;
end
if nargin<8 || isempty(bidi_comp)
    bidi_comp=false;
end
if bidi_comp
batchData=batchData([1:2:end 2:2:end],:,:);
end
[ly,lx]=size(batchData,1:2);
m=mean(batchData,1:2);
batchData=batchData-m;
if K_b>0
        if size(batchData,3)>1000
        [U,S,V] = bksvd(reshape(batchData,[],size(batchData,3)), K_b);
        else
            %for <1000 frames, 512x512 pixels, svdecon is faster
        [U,S,V] = svdecon(reshape(batchData,[],size(batchData,3)));
            U=U(:,1:K_b);
            S=S(1:K_b,1:K_b);
            V=V(:,1:K_b);
        end
        corrMap=fft2(reshape(U,ly,lx,[]));

       % if phaseCorrelation
       %      corrMap = bsxfun(@times, ...
       %          corrMap ./ (eps0 + abs(corrMap)) .* fhg, cfRefImg);
       %  else
            corrMap = bsxfun(@times, corrMap, cfRefImg);
        % end

        corrClip = real(ifft2(corrMap));
        corrClip = fftshift(fftshift(corrClip,1),2);

        % --- subpixel estimation ---
            Mode_cc0 = corrClip( ...
                floor(ly/2)+1+(-lcorr:lcorr), ...
                floor(lx/2)+1+(-lcorr:lcorr), :);

Mode_cc0 = reshape(Mode_cc0,[],K_b);
% VS = V * S;
% cc0=Mode_cc0*VS';
cc0=Mode_cc0*S*V';
cc0=reshape(cc0,lcorr*2+1,lcorr*2+1,[]);

% cc0=zeros(lcorr*2+1,lcorr*2+1,nf);
% corrClip=reshape(corrClip,[],size(corrClip,3));
% VS = V * S;
% for i=1:nf
% cc0_temp=reshape(corrClip*VS(i,:)',ly,lx,[]);
% cc0(:,:,i)=cc0_temp( ...
%                 floor(ly/2)+1+(-lcorr:lcorr), ...
%                 floor(lx/2)+1+(-lcorr:lcorr), :);
% end
% cc0=reshape(cc0,lcorr*2+1,lcorr*2+1,[]);

else
            corrMap = fft2(batchData);

        if phaseCorrelation
            corrMap = bsxfun(@times, ...
                corrMap ./ (eps0 + abs(corrMap)) .* fhg, cfRefImg);
        else
            corrMap = bsxfun(@times, corrMap, cfRefImg);
        end

        corrClip = real(ifft2(corrMap));
        corrClip = fftshift(fftshift(corrClip,1),2);
            cc0 = corrClip( ...
                floor(ly/2)+1+(-lcorr:lcorr), ...
                floor(lx/2)+1+(-lcorr:lcorr), :);
end


end