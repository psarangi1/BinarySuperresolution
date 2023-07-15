function spk_est=decode_spike_fir(y_down,h,flip,m,N,p)
    %% Function for Decoding uniformly downsampled measurements from FIR filter

    spk_est=binornd(1,p,[N,1]);
    r=length(h);
    for i=1:length(y_down)
       
            if i==1
                spk_est(1:r)=beta_exp_decon(y_down(i),h,flip); %Decode first block differently
            else
                %%Decode other blocks of length m
                temp=beta_exp_decon(y_down(i),h,flip); 
                LL=length((i-1)*m+1:min((i-1)*m+r,N));
                spk_est((i-1)*m+1:min((i-1)*m+r,N))=temp(1:LL);
            end
       
    end
end