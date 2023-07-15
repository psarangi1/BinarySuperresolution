%%Script for Evaluating Spike Recovery Performance for different
%%undersampling factor
%Evaluates performance for spike recovery (using F-score metric)

clear all;
close all;
addpath(genpath('Functions/'))
mc_iter=20;
N=100;
mmax=10;

%% Model Parameters
alpha=0.5; %AR(1) filter
amp=1; %Ground truth spike parameter
p=0.35; %Spiking probability of high-rate binary spikes
sig=0.001; %Noise power

p1=length(alpha)+1; 

G_alpha=toeplitz([1,-alpha,zeros(1,N-p1)],[1,zeros(1,N-1)]); %Filtering matrix (Toepitz)
G1=inv(G_alpha);


%% Initialize result variables with 0
false_pos_l1_relax=zeros(mmax,mc_iter);
hit_rate_l1_relax=zeros(mmax,mc_iter);
err_rate_l1_relax=zeros(mmax,mc_iter);


false_pos_l1_box_relax=zeros(mmax,mc_iter);
hit_rate_l1_box_relax=zeros(mmax,mc_iter);
err_rate_l1_box_relax=zeros(mmax,mc_iter);

false_pos_bin_relax=zeros(mmax,mc_iter);
hit_rate_bin_relax=zeros(mmax,mc_iter);
err_rate_bin_relax=zeros(mmax,mc_iter);


count_err_=zeros(mmax,mc_iter);
count_err_post_l1=zeros(mmax,mc_iter);
count_directl1=zeros(mmax,mc_iter);



%% Loop for varying different downsampling factors
for m=3:mmax
    
    %Prep binary search (executed only once for a given set of model
    %parameters):
    [c_fit,sol]=binary_prep(alpha,m,amp); %Function Constructs the set Theta_alpha
   
    %% Monte carlo (for each undersampling factor)
    parfor mc=1:mc_iter
        %%Generate Synthetic measurements

        %Ground truth spike generation 
        x=amp*binornd(1,p,[N,1]); %Bernoulli with paraemter p
        x(1)=0;
        
        %High-rate Filter output (noiseless)
        yhi=G1*x;

        %Low-rate measurements after m-fold undersampling
        ylo=yhi(1:m:end);
        M=length(ylo);
          
        %Noisy low-rate samples
        z_n = ylo+sig*randn(M,1);
           

        %% Begin Decoding
        x_est=SpikeDecodeAR(z_n,alpha,c_fit,sol,amp,m);
    
        %% Prepare the measurement matrix for box-constraint
        DD=eye(N);
        D=DD(1:m:end,:);
        y_m = D*G1*x;
        TT=1+(M-1)*m;
        zz = D*G1(:,1:TT)*x(1:TT);
        %norm(y_m(1:m-1)-z(1:m-1))
        M=length(zz);

        H=D*G1(:,1:TT); %Measurement matrix for l1 minimization techniques

       
        %Use ground truth noise level to set threshold for 
        n_est=norm(zz-z_n);

        %%l1 minimization with non-negativity
        x_est_l1=l1_recovery(H,z_n,TT,n_est);

        %%l1 minimization with box-constraint
        x_est_l1_box=l1_recovery_box(H,z_n,TT,n_est,amp);
        
        %%Support extraction
        x_est_l1_supp=zeros(TT,1);
        x_est_l1_supp(x_est_l1>amp/2)=amp; %Threshold based on amplitude
        
        x_est_l1_box_supp=zeros(TT,1);
        x_est_l1_box_supp(x_est_l1_box>amp/2)=amp; %Threshold based on amplitude
        
        %Timing Extraction
        t_hi=0:1e-3:(TT-1)*1e-3;
        gt_spk=t_hi(logical(x(1:TT)));
        t_l1=t_hi(logical(x_est_l1_supp));
        t_l1_box=t_hi(logical(x_est_l1_box_supp));
        t_greedy=t_hi(logical(x_est));
             

        %% Evaluate performance by computing F-score (Borrowed from MLspike)
        j_d=2;
        t_delay=1e-3*(j_d);
        
        cost = struct( ...
                                'miss0',        1, ...      % cost of isolated miss
                                'miss1',        .5, ...     % cost of non-isolated miss
                                'flsp0',        1, ...      % cost of isolated false-positive
                                'flsp1',        .5, ...     % cost of non-isolated false-positive
                                'nneighbors',   2, ...      % # neighbors that define a non-isolated spike
                                'timeburst',    0.5, ...     % time constant defining non-isolated spikes
                                'timematch',    .0, ...     % maximal time distance for a perfect match
                                'timedelay',    t_delay, ...     % maximal time distance for matching spikes
                                'maxdelay',     2, ...      % cost of a match at maximal time distance
                                'sides',        [] ...      % set the start and end times of the acquisition for a 'smart' count that minimizes side effects
                                );
        %Error computation for Proposed Algorithm
        [d_bin,dsec_bin]=spk_distance(gt_spk,t_greedy,cost);
        err_rate_bin_relax(m,mc)=f1score(dsec_bin.count.miss,dsec_bin.count.falsep,dsec_bin.count.true,dsec_bin.count.detections);
        hit_rate_bin_relax(m,mc)=dsec_bin.count.trued/dsec_bin.count.true;
        false_pos_bin_relax(m,mc)=dsec_bin.count.falsep/dsec_bin.count.detections;
        d_bin_vec(m,mc)=d_bin;
        
        %Error computation for l1 minimization
        [d_l1,dsec_l1]=spk_distance(gt_spk,t_l1,cost);
        err_rate_l1_relax(m,mc)=f1score(dsec_l1.count.miss,dsec_l1.count.falsep,dsec_l1.count.true,dsec_l1.count.detections);
        hit_rate_l1_relax(m,mc)=dsec_l1.count.trued/dsec_l1.count.true;
        false_pos_l1_relax(m,mc)=dsec_l1.count.falsep/dsec_l1.count.detections;
        d_l1_vec(m,mc)=d_l1;

        %Error computation for l1 minimization with box-constraint
        [d_l1_box,dsec_l1_box]=spk_distance(gt_spk,t_l1_box,cost);
        err_rate_l1_box_relax(m,mc)=f1score(dsec_l1_box.count.miss,dsec_l1_box.count.falsep,dsec_l1_box.count.true,dsec_l1_box.count.detections);
        hit_rate_l1_box_relax(m,mc)=dsec_l1_box.count.trued/dsec_l1_box.count.true;
        false_pos_l1_box_relax(m,mc)=dsec_l1_box.count.falsep/dsec_l1_box.count.detections;
        d_l1_box_vec(m,mc)=d_l1_box; 
    end
end

%%Plotting the results
h1=figure;
plot(1:mmax,mean(1-err_rate_bin_relax,2)','LineWidth',2,'Marker','o','MarkerSize',10)
hold on
plot(1:mmax,mean(1-err_rate_l1_relax,2)','LineWidth',2,'Marker','*','MarkerSize',10)
plot(1:mmax,mean(1-err_rate_l1_box_relax,2)','LineWidth',2,'Marker','*','MarkerSize',10)
xlim([3,mmax])
xlabel('Undersampling Factor (D)','FontSize',15)
ylabel('F-score','FontSize',15)
legend('Proposed','l_1','Box-l_1','FontSize',15)
pbaspect([5,2,1])
grid minor
pbaspect([5,2,1])
%ylim([0.2,1.2])
xlim([3,10])
title_text = "p="+num2str(p)+",\alpha="+num2str(alpha);
title(title_text,'FontSize',12)




       
    