%%Script for Evaluating Spike Recovery Performance for different noise
%%powers
%Evaluates performance of both Spike recovery (using F-score metric)
%and Spike count estimation (using l1 distance)

%For speeding up - Disable the Box-constraint comparison by setting the
%flag box=0


clear all;
close all;
addpath(genpath('Functions/'))

alpha_vec=[0.5,0.9]; %Two different values of alpha
N=100;
mmax=10;

%% Downsampling factor
m=5;
amp=2; %Binary Spike Amplitude (will control SNR)
p=0.35; %Spiking probability

mc_iter=20;
time_box=zeros(2,mc_iter);
time_bin=zeros(2,mc_iter);

box=1; %Box Constraint comparison is on (will be slow due to CVX)

%% Sweep over values of alpha in the vector alpha_vec ang generate results
for alpha_i=1:length(alpha_vec)
    alpha=alpha_vec(alpha_i);
    
   
    
   
    %% Initialize results variables with 0
    false_pos_l1_relax=zeros(mmax,mc_iter);
    hit_rate_l1_relax=zeros(mmax,mc_iter);
    err_rate_l1_relax=zeros(mmax,mc_iter);

    false_pos_l1_relax_box=zeros(mmax,mc_iter);
    hit_rate_l1_relax_box=zeros(mmax,mc_iter);
    err_rate_l1_relax_box=zeros(mmax,mc_iter);
    
    
    false_pos_bin_norelax=zeros(mmax,mc_iter);
    hit_rate_bin_norelax=zeros(mmax,mc_iter);
    err_rate_bin_norelax=zeros(mmax,mc_iter);
    
   
    d_l1_vec =zeros(mmax,mc_iter);
    d_bin_vec =zeros(mmax,mc_iter);
    d_l1_box_vec =zeros(mmax,mc_iter);

    count_err_=zeros(mmax,mc_iter);
    count_directl1=zeros(mmax,mc_iter);
    count_directl1_box=zeros(mmax,mc_iter);
    

 
    p1=length(alpha)+1;
    G_alpha=toeplitz([1,-alpha,zeros(1,N-p1)],[1,zeros(1,N-1)]);
    G1=inv(G_alpha);

    DD=eye(N);
    D=DD(1:m:end,:);
    
    %Prep binary search (executed only once for a given set of model
    %parameters): Construction of set Theta_alpha
    [c_fit,sol]=binary_prep(alpha,m,amp);

    %Values for sweeping Noise (will also control SNR together with p and
    %amp)
    n_vec=logspace(0,-5,mmax);
    
    %% Loop over the different noise levels
    for n_i=1:length(n_vec)

        sig=sqrt(n_vec(n_i)); %Set the noise power
       
        %%Monte Carlo for a given noise power (Different spike and noise
        %%realizations)
        parfor mc=1:mc_iter

            %% Generate Synthetic measurements

            %Ground truth spike generation 
            x=amp*binornd(1,p,[N,1]); %Bernoulli
            x(1)=0;

            %High-rate Filter output (noiseless)
            y=G1*x;

            %Low-rate measurements after m-fold undersampling
            z=y(1:m:end);
            M=length(z);
            
            %Noisy low-rate samples
            z_n = z+sig*randn(M,1); %Noise variance chosen to be sig^2
                
            tic
            %%Binary Spike Decoding Algorithm 
            x_est=SpikeDecodeAR(z_n,alpha,c_fit,sol,amp,m);
            time_bin(n_i,mc)=toc; %Computes Run-time for decoding
            

            %Extract Count Information from Decoded Spikes 
            tru_count=zeros(1,length(z_n)); %Initialize the ground truth count vector
            count_est=zeros(1,length(z_n)); %Initialize the count estimate vector
            
            %Extract the count for the first block
            tru_count(1)=x(1)/amp;
            count_est(1)=x_est(1)/amp;
            
            %% Extract the count for remaining M-1 blocks
            for jj=1:length(z_n)-1           
                tru_count(jj+1)=sum(x((jj-1)*m+2:jj*m+1)/amp);
                count_est(jj+1)=sum(x_est((jj-1)*m+2:jj*m+1)/amp);
            end

            
            
            %Extract Timing information for F-score computation
            TT=1+(M-1)*m;

            t_hi=0:1e-3:(TT-1)*1e-3;
            gt_spk=t_hi(logical(x(1:TT)));
            t_greedy=t_hi(logical(x_est));
    
            H=D*G1(:,1:TT);
    
            if box==1

                n_est=norm(z-z_n); %Use oracle noise level for realxation algorithms

                %% L1 minimization with non-negativity
                x_est_l1=l1_recovery(H,z_n,TT,n_est);
    
                %% L1 minimization with Box-relaxation for binary
                tic
                x_est_l1_box=l1_recovery_box(H,z_n,TT,n_est,amp);
                time_box(n_i,mc)= toc; %Timing computation (Significantly slower compared to binary)
    
                %Support Recovery for relaxation schemes 
                x_est_l1_supp=zeros(TT,1);
                x_est_l1_supp(x_est_l1>amp/2)=amp;
        
                x_est_l1_box_supp=zeros(TT,1);
                x_est_l1_box_supp(x_est_l1_box>0.5*amp)=amp;
                
                %Extract Count Information from Decoded Spikes obtained by
                %relaxation schemes
                count_l1=zeros(1,length(z_n));
                count_l1_box=zeros(1,length(z_n));
    
                count_l1(1)=round(x_est_l1(1));
                count_l1_box(1)=round(x_est_l1_box(1));
                for jj=1:length(z_n)-1
                    
                    count_l1(jj+1)=round(sum(x_est_l1((jj-1)*m+2:jj*m+1)/amp));
                    count_l1_box(jj+1)=round(sum(x_est_l1_box((jj-1)*m+2:jj*m+1)/amp));
                end

                t_l1=t_hi(logical(x_est_l1_supp)); %Extract support for F-score
                t_l1_box=t_hi(logical(x_est_l1_box_supp)); %Extract support for F-score
            end
        
      
            
            
         
                 
    
            %Evaluate performance by computing F-score (Borrowed from MLspike)
            j_d=2;
            t_delay=1e-3*(j_d);
            cost = struct( ...
                                    'miss0',        1, ...      % cost of isolated miss
                                    'miss1',        .5, ...     % cost of non-isolated miss
                                    'flsp0',        1, ...      % cost of isolated false-positive
                                    'flsp1',        .5, ...     % cost of non-isolated false-positive
                                    'nneighbors',   2, ...      % # neighbors that define a non-isolated spike
                                    'timeburst',    .5, ...     % time constant defining non-isolated spikes
                                    'timematch',    .0, ...     % maximal time distance for a perfect match
                                    'timedelay',    t_delay, ...     % maximal time distance for matching spikes
                                    'maxdelay',     2, ...      % cost of a match at maximal time distance
                                    'sides',        [] ...      % set the start and end times of the acquisition for a 'smart' count that minimizes side effects
                                    );
            %Spike Recovery Error computation for Proposed Algorithm
            [d_bin,dsec_bin]=spk_distance(gt_spk,t_greedy,cost);
            err_rate_bin_norelax(n_i,mc)=f1score(dsec_bin.count.miss,dsec_bin.count.falsep,dsec_bin.count.true,dsec_bin.count.detections);
            hit_rate_bin_norelax(n_i,mc)=dsec_bin.count.trued/dsec_bin.count.true;
            false_pos_bin_norelax(n_i,mc)=dsec_bin.count.falsep/dsec_bin.count.detections;
            d_bin_vec(n_i,mc)=d_bin;
                
            %Spike Recovery Error computation for Relaxation based
            %Algorithms (L1 and L1+Box)
            if box==1
                [d_l1,dsec_l1]=spk_distance(gt_spk,t_l1,cost);
                err_rate_l1_relax(n_i,mc)=f1score(dsec_l1.count.miss,dsec_l1.count.falsep,dsec_l1.count.true,dsec_l1.count.detections);
                hit_rate_l1_relax(n_i,mc)=dsec_l1.count.trued/dsec_l1.count.true;
                false_pos_l1_relax(n_i,mc)=dsec_l1.count.falsep/dsec_l1.count.detections;
                d_l1_vec(n_i,mc)=d_l1;
    
                 [d_l1_box,dsec_l1_box]=spk_distance(gt_spk,t_l1_box,cost);
                 err_rate_l1_relax_box(n_i,mc)=f1score(dsec_l1_box.count.miss,dsec_l1_box.count.falsep,dsec_l1_box.count.true,dsec_l1_box.count.detections);
                 hit_rate_l1_relax_box(n_i,mc)=dsec_l1_box.count.trued/dsec_l1_box.count.true;
                 false_pos_l1_relax_box(n_i,mc)=dsec_l1_box.count.falsep/dsec_l1_box.count.detections;
                 d_l1_box_vec(n_i,mc)=d_l1_box;
            end

            %Evaluate performance of count estimation by computing l1
            %distance 
            count_err_(n_i,mc)=norm(tru_count-count_est,1); %Proposed algorithm
            if box==1
                count_directl1(n_i,mc)=norm(tru_count-count_l1,1); %L_1 minimization with NN
                count_directl1_box(n_i,mc)=norm(tru_count-count_l1_box,1); %L_1 minimization with box
            end
        end
        
        fprintf("Completed %d\n",n_i)
    end
    [status, msg, msgID] = mkdir('Simulation Data');
    %Save data
    save(strcat('Simulation Data/Spike_and_Count_Recovery_L1_comparison_mc_noise_vary_','m_',num2str(m),'_alpha_',num2str(100*alpha),'_prob_',num2str(100*p),'_amp_',num2str(amp)))
    
    %%Plot Generation 
    figure(1) %F-score spike recovery
    semilogx(n_vec,mean(1-err_rate_bin_norelax,2)','LineWidth',2,'Marker','diamond','MarkerSize',10)
    hold on
    if box==1
        
        semilogx(n_vec,mean(1-err_rate_l1_relax_box,2)','LineWidth',2,'Marker','diamond','MarkerSize',10)
    end
    xlabel('Noise Level','FontSize',14)
    ylabel('F-score','FontSize',14)
    if alpha_i==2
        if box==1
            legend('Proposed (\alpha=0.5)','l_1 Box (\alpha=0.5)', 'Proposed (\alpha=0.9)','l_1 Box (\alpha=0.9)','FontSize',12,'NumColumns',2,Location=[0.64,0.62,0.1,0.1])
        else
            legend('Proposed (\alpha=0.5)','Proposed (\alpha=0.9)')
        end
        grid minor
    end
    

    figure(2) %L1 error for spike count estimation
    loglog(n_vec,mean(count_err_+1e-10,2),'LineWidth',2,'Marker','diamond','MarkerSize',10)
    hold on
    if box==1
        
        loglog(n_vec,mean(count_directl1_box+1e-10,2),'LineWidth',2,'Marker','diamond','MarkerSize',10)
    end
    xlabel('Noise Level','FontSize',14)
    ylabel('Count Estimation Error','FontSize',14)
    if alpha_i==2
        if box==1
            legend('Proposed (\alpha=0.5)','l_1 Box (\alpha=0.5)', 'Proposed (\alpha=0.9)','l_1 Box (\alpha=0.9)','FontSize',12,'NumColumns',2,Location=[0.64,0.62,0.1,0.1])
        else
            legend('Proposed (\alpha=0.5)','Proposed (\alpha=0.9)')
        end
        grid minor
    end
    
end
       
    