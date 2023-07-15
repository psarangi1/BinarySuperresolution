
clear;
close all;
rng(1);
addpath(genpath('Functions'))

load('GCaMP6f_chen.mat') %Load the desired dataset

%Demo Script will generate results for the following 3 examples
c_vec=[32,1,2,10]; 

err_rate_bin = zeros(length(c_vec),1);
hit_rate_bin = zeros(length(c_vec),1);
t_er_bin = zeros(length(c_vec),1);

err_rate_bin_d = zeros(length(c_vec),1);
hit_rate_bin_d = zeros(length(c_vec),1);
t_er_bin_d = zeros(length(c_vec),1);

%The plots will zoom into these regions to showcase improved spike recovery
%compared to OASIS
t1_vec=[9,80,170,20];
t2_vec=[9.8,82,172,22];


% Undersampling parameter: Target High-resolution grid on which spikes will be inferred
m=10;

flag=0; %Baseline correction algorithm
bp=5;

for i=1:length(c_vec)
    tau=0.5; %Decay time constant(indicator specific)
    c_i=c_vec(i);

    t_frame_all_=t_frame_all(:,c_i).';
    T_ca=t_frame_all_(2)-t_frame_all_(1);
    fmean_all_=double(fmean_all(:,c_i));
    
    % Baseline Correction
    len2=256;
    b_flag=0; %Choose 0 for Mitani and Komiyama algorithm and 1 for Percentile based correction
    if b_flag==1
        %Percentile based baseline removal
        for l_i=1:floor(size(fmean_all,1)/len2)-1
            if l_i==floor(size(fmean_all,1)/len2)
                b=prctile(fmean_all(l_i*len2+1:end,i),bp);
                y_corrected(l_i*len2+1:end)=(fmean_all(l_i*len2+1:end,i)-b);
            else
                b=prctile(fmean_all(l_i*len2+1:(l_i+1)*len2,i),bp);
                y_corrected(l_i*len2+1:(l_i+1)*len2)=(fmean_all(l_i*len2+1:(l_i+1)*len2,i)-b);
            end
        end
    elseif b_flag==0
        b = baseline_kde(fmean_all_,bp,100,20); %Baseline correction algorithm by Mitani and Komiyama
        y_corrected=fmean_all_-b;
    end

    y_corrected2=double(y_corrected); %Baseline corrected calcium waveform
    
    t_frame_all_=t_frame_all_(1:length(y_corrected2));

   

    
    T_hi=(t_frame_all_(2)-t_frame_all_(1))/(m);
    lambda=0;

    g=exp(-(t_frame_all_(2)-t_frame_all_(1))/tau);
    % OASIS Output
    [c_oasis2, s_oasis, options] = deconvolveCa(y_corrected2, 'ar1',g, 'foopsi', 'lambda', lambda,'optimize_pars'); 
    
    s_oasis_th=s_oasis;
    th=0.05*max(s_oasis);
    s_oasis_th(s_oasis_th<th)=0;
    % Denoise/Pre-process (using OASIS)
    y_est_den=denoise_AR(y_corrected2,1/T_ca,tau,0);
            
    %AR model parameter
    tau=-(t_frame_all_(2)-t_frame_all_(1))/log(options.g);
    alpha=(options.g)^(1/m);


    %Full deconvolution pipeline
    
    [spike_est22,y_est_ca,t_est_22,y_up22,amp_algo]=deconv_full(y_corrected2,y_est_den,alpha,m,t_frame_all);
   

    %Pruning 
    [spike_est_bin_supp,t_bin_est,y_est_bin]=Pruning_algo(y_up22,amp_algo,t_est_22,tau);
        
           
   
    %Performance Evaluation

    %Ground Truth
    t_ephys_all_ci=t_ephys_all(detected_spikes_all(:,c_i)>0,c_i);

    cost = struct( ...
        'miss0',        1, ...      % cost of isolated miss
        'miss1',        .5, ...     % cost of non-isolated miss
        'flsp0',        1, ...      % cost of isolated false-positive
        'flsp1',        .5, ...     % cost of non-isolated false-positive
        'nneighbors',   2, ...      % # neighbors that define a non-isolated spike
        'timeburst',    .5, ...     % time constant defining non-isolated spikes
        'timematch',    .0, ...     % maximal time distance for a perfect match
        'timedelay',    0.2, ...     % maximal time distance for matching spikes
        'maxdelay',     2, ...      % cost of a match at maximal time distance
        'sides',        [] ...      % set the start and end times of the acquisition for a 'smart' count that minimizes side effects
        );

        [~,dsec_bin]=spk_distance(t_ephys_all_ci,t_bin_est,cost);

        err_rate_bin(i)=f1score(dsec_bin.count.miss,dsec_bin.count.falsep,dsec_bin.count.true,dsec_bin.count.detections);
        hit_rate_bin(i)=dsec_bin.count.trued/dsec_bin.count.true;
        t_er_bin(i)=median(abs(dsec_bin.delays));


        fprintf('File %s: Hit rate Bin:: %f,  Error rate Bin:: %f, Timing:: %f \n',num2str(i),hit_rate_bin(i),err_rate_bin(i),t_er_bin(i));
        

         %Plot Generation (Zooming into interesting regions)

        t1 = t1_vec(i);
        t2 = t2_vec(i);
        
        
        
        figure;
        ax1=subplot(3,1,1);
        plot(t_frame_all(:,c_i),fmean_all(:,c_i),'r','LineWidth',0.4),set(gca,'FontSize',20,'XTick',[], 'YTick', [],'box','off'),xlim([t1,t2])
        title('B-OASIS (Input at 60Hz)')
        ax2=subplot(3,1,2);
        stem(t_bin_est,ones(length(t_bin_est),1),'b','Marker','none','LineWidth',2)
        xlim([t1,t2])
        xticks(t1:(t2-t1)/10:t2)
        set(gca, 'YTick', [],'box','off','fontweight','bold','fontsize',12)
        
        ax3=subplot(3,1,3);
        stem(t_ephys_all_ci,ones(length(t_ephys_all_ci),1),'k','Marker','none','LineWidth',1.2)
        
        xlim([t1,t2])
        linkaxes([ax1,ax2,ax3],'x')
        xticks(t1:(t2-t1)/10:t2)
        set(gca,'YTick', [],'box','off','fontweight','bold','fontsize',12)
        %xticks([t1:(t2-t1)/10:t2])
        
        
        
        figure;
        ax1=subplot(3,1,1);
        plot(t_frame_all(:,c_i),fmean_all(:,c_i),'r','LineWidth',0.4),set(gca,'FontSize',20,'XTick',[], 'YTick', [],'box','off'),xlim([t1,t2])
        title('OASIS (Input at 60Hz)')
        ax2=subplot(3,1,2);
        stem(t_frame_all(:,c_i),s_oasis_th,'b','Marker','none','LineWidth',2),
        xlim([t1,t2])
        xticks(t1:(t2-t1)/10:t2)
        set(gca,'YTick', [],'box','off','fontweight','bold','fontsize',12)
        
        ax3=subplot(3,1,3);
        stem(t_ephys_all_ci,ones(length(t_ephys_all_ci),1),'k','Marker','none','LineWidth',1.2),
        
        xlim([t1,t2])
        linkaxes([ax1,ax2,ax3],'x')
        xticks(t1:(t2-t1)/10:t2)
        set(gca,'YTick', [],'box','off','fontweight','bold','fontsize',12)
        %xticks([t1:(t2-t1)/10:t2])
        
        
        


        %Same Data with synthetic downsampling to 30Hz (by subselecting
        %every other sample of the fluorescence signal)

        t_frame_all_=t_frame_all(1:2:end,c_i).';
        T_ca=t_frame_all_(2)-t_frame_all_(1);
        fmean_all_=double(fmean_all(1:2:end,c_i));
    
        y_corrected2=double(y_corrected(1:2:end)); %Baseline corrected calcium waveform
    

        T=(t_frame_all_(2)-t_frame_all_(1))/(m);
        lambda=0;
        
        g=exp(-(t_frame_all_(2)-t_frame_all_(1))/tau);
        % OASIS Output
        [c_oasis2d, s_oasisd, optionsd] = deconvolveCa(y_corrected2, 'ar1',g, 'foopsi', 'lambda', lambda,'optimize_pars'); 
        
         s_oasis_thd=s_oasisd;
        thd=0.05*max(s_oasisd);
        s_oasis_thd(s_oasis_thd<thd)=0;
        % Denoise/Pre-process (using OASIS)
        y_est_dend=denoise_AR(y_corrected2,1/T_ca,tau,0);
            
        %AR model parameter
        tau=-(t_frame_all_(2)-t_frame_all_(1))/log(optionsd.g);
        alphad=(optionsd.g)^(1/m);
    
    
        %Full deconvolution pipeline
        
        [spike_est22d,y_est_cad,t_est_22d,y_up22d,amp_algod]=deconv_full(y_corrected2,y_est_dend,alphad,m,t_frame_all_);
       
    
        %Pruning 
        [spike_est_bin_suppd,t_bin_estd,y_est_bind]=Pruning_algo(y_up22d,amp_algod,t_est_22d,tau);
            
               
        %Ground Truth        
        gt_spikes_=detected_spikes_all(:,c_i);
        t_hi_sp=t_ephys_all(:,c_i);
            
        %Performance Evaluation
        t_ephys_all_ci=t_ephys_all(detected_spikes_all(:,c_i)>0,c_i);
    
        cost = struct( ...
            'miss0',        1, ...      % cost of isolated miss
            'miss1',        .5, ...     % cost of non-isolated miss
            'flsp0',        1, ...      % cost of isolated false-positive
            'flsp1',        .5, ...     % cost of non-isolated false-positive
            'nneighbors',   2, ...      % # neighbors that define a non-isolated spike
            'timeburst',    .5, ...     % time constant defining non-isolated spikes
            'timematch',    .0, ...     % maximal time distance for a perfect match
            'timedelay',    0.2, ...     % maximal time distance for matching spikes
            'maxdelay',     2, ...      % cost of a match at maximal time distance
            'sides',        [] ...      % set the start and end times of the acquisition for a 'smart' count that minimizes side effects
            );
    
            [d_bin,dsec_bin]=spk_distance(t_ephys_all_ci,t_bin_estd,cost);
    
            err_rate_bin_d(i)=f1score(dsec_bin.count.miss,dsec_bin.count.falsep,dsec_bin.count.true,dsec_bin.count.detections);
            hit_rate_bin_d(i)=dsec_bin.count.trued/dsec_bin.count.true;
            t_er_bin_d(i)=median(abs(dsec_bin.delays));


            fprintf('File %s: Hit rate Bin:: %f,  Error rate Bin:: %f, Timing:: %f \n',num2str(i),hit_rate_bin_d(i),err_rate_bin_d(i),t_er_bin_d(i));
            
            
            %Plot Generation (Zooming into interesting regions)

           

            figure;
            ax31=subplot(3,1,1);
            plot(t_frame_all_,fmean_all_,'r','LineWidth',0.4),
            set(gca,'FontSize',20,'XTick',[], 'YTick', [],'box','off'),
            xlim([t1,t2])
            title('B-OASIS (Input at 30Hz)')
            
            ax32=subplot(3,1,2);
            stem(t_bin_estd,ones(length(t_bin_estd),1),'b','Marker','none','LineWidth',2)
            
            xlim([t1,t2])
            xticks(t1:(t2-t1)/10:t2)
            set(gca, 'YTick', [],'box','off','fontweight','bold','fontsize',12)
            
            
            
            ax33=subplot(3,1,3);
            stem(t_ephys_all_ci,ones(length(t_ephys_all_ci),1),'k','Marker','none','LineWidth',1.2)
            
            xlim([t1,t2])
            linkaxes([ax31,ax32,ax33],'x')
            
            xticks(t1:(t2-t1)/10:t2)
            set(gca, 'YTick', [],'box','off','fontweight','bold','fontsize',12)
            
            
            
            


            figure;
            ax41=subplot(3,1,1);
            plot(t_frame_all_,fmean_all_,'r','LineWidth',0.4),set(gca,'FontSize',20,'XTick',[], 'YTick', [],'box','off'),xlim([t1,t2])
            title('OASIS (Input at 30Hz)')
            ax42=subplot(3,1,2);
            stem(t_frame_all_,s_oasis_thd,'b','Marker','none','LineWidth',2),
            xlim([t1,t2])
            %set(gca, 'YTick', [],'box','off','fontweight','bold','fontsize',12)
            xticks(t1:(t2-t1)/10:t2)
            set(gca, 'YTick', [],'box','off','fontweight','bold','fontsize',12)
            
            ax43=subplot(3,1,3);
            stem(t_ephys_all_ci,ones(length(t_ephys_all_ci),1),'k','Marker','none','LineWidth',1.2)
            
            xlim([t1,t2])
            linkaxes([ax41,ax42,ax43],'x')
            
            xticks(t1:(t2-t1)/10:t2)
            set(gca, 'YTick', [],'box','off','fontweight','bold','fontsize',12)
            
           
end
