function [spike_est,y_est_ca,t_est_,y_up,amp_algo]=deconv_full(y_corrected,y_est_den,alpha,m,t_frame_all)
%%Full code for the deconvolution pipeline
% Inputs: 
% y_corrected,y_est_MR2 : Low rate waveforms
%alpha          : AR model parameter
%m              : Multirate factor
%t_frame_all    : low rate frame timings
%Output:
%spike_est      : High rate spikes
%t_est_         : High resolution timing grid
%amp_algo       : estimated amplitude
%y_up           : estimate of the high rate calcium
    y_diff_den=y_est_den-alpha^m*[0;y_est_den(1:end-1)];
    y_diff_n=y_corrected-alpha^m*[0;y_corrected(1:end-1)];
    n_est=abs(y_diff_den-y_diff_n);
    [pks,locs] = findpeaks(y_diff_den);
    large_pks=find(pks>0.05*max(y_diff_den));
    large_pks=intersect(large_pks,find(pks<0.8*max(y_diff_den)));
    pks=pks(large_pks);
    locs=locs(large_pks);
    
    id_smalln=find(y_est_den>0.2*max(y_est_den));
    id_good=intersect(id_smalln-1,locs);

    y_d=y_est_den(2:end)-alpha^m*y_est_den(1:end-1);

    [vv3,id_est4]=sort(y_d(id_good-1),'descend');
             
    topn=min(10,length(vv3));
   
    
    %Amplitude estimation
    for th=1.0:0.1:10
        flag=1;
        amp12=match_amp(vv3(1:topn).',alpha,m);
        amp_est12=resolve_amp(amp12,th,0);
        [h1,h2]=hist(amp_est12);
        [m_val,idmax]=max(h1);
        amp_est(1)=h2(idmax);
        if isempty(amp_est12)
           flag =0;
        end
        amp_est12=resolve_amp(amp12,th,1);
                
        [h1,h2]=hist(amp_est12);
        [m_val,idmax]=max(h1);
        amp_est(2)=h2(idmax);
        if isempty(amp_est12)
            flag=0;
        end
        if flag==1
                break;
        end
    end
    
    
    amp_algo=median(amp_est);
    
    [c_fit,sol]=binary_prep(alpha,m,amp_algo);
    
    [spike_est,y_est_ca,t_est_]=SpikeDecodeAR_real(y_est_den,t_frame_all,c_fit,sol,m,alpha);
    
    y_up(1)=spike_est(1);
    for i=2:length(spike_est)
        y_up(i)=alpha*y_up(i-1)+spike_est(i);
    end
    
  
    
end