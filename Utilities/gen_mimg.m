function refImg=gen_mimg(data,whichch,n_ch,n_Frames_template)
        div=size(data,3)/n_ch/n_Frames_template;
        first=floor(div/2)+1;
        % div=floor(div);
        frames=(round(first:div:nf/n_ch)-1)*n_ch+whichch;
        % refImg=gen_template(data(:,:,whichch:n_ch:n_ch*100),min(100,nf/n_ch));
        small_corr_stack=reg2P_standalone(data(:,:,frames),nanmean(data(:,:,frames),3));
        refImg=nanmean(small_corr_stack,3);
end