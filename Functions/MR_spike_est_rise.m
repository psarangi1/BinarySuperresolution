function [spike_est_MR2,spike_est_MR2_supp,y_est_MR,t_est_]=MR_spike_est_rise(F,tau1,tau2,m,m1,l_f,p_g,p,nz,len1,t_frame)
%%For Pruning spikes only decay and rise time constant
% ============================== INPUTS ============================== 
%   F           : T X 1 Calcium signal
%
%   p           : p-norm. p lies in [0,1]
%
%   p_g         : Noise Level for FOCUSS algorithm
%
%   tau1        : Rise time
%
%   tau2        : Decay Time
%
%   m           : Upsampling factor
%
%   m1          : Interpolation factor
%                   
%   len1        : Size of each window
%
%   nz          : Previous window for initilization
%
%   t_frame     : Calcium Sampling intervals
% ==============================  OUTPUTS ============================== 
% spike_est_MR2        : Estimated spikes m*Tx1 
% spike_est_MR2_supp   : Default threshold (not adaptive)
% y_est_MR             : Fit generated by algortihm Tx1
% t_est_               : Time vector for spikes mTx1
%
% ============== Examples of Commands ===============
% [Example 1]
%   Parameters shown in Demo script
%   [spike_est_MR22,~,y_est_MR2,t_est_2]=MR_spike_est_rise3(y_corrected,trise,tau,m,m1,l_f,l_g,p,pad_len,len1,t_frame_all);


    T=(t_frame(2)-t_frame(1))/(m);


    %Equivalent Digital filter representation in matrix form at a
    %given rate which is m times fater than calcium's sampling
    %frequency
    alpha1=(T-2*tau2);
    alpha2=(T+2*tau2);
    beta1=(T-2*tau1);
    beta2=(T+2*tau1);
    C_circ1=toeplitz([T^2*(tau2-tau1),2*T^2*(tau2-tau1),T^2*(tau2-tau1),zeros(1,len1*m+nz*m-3)],[T^2*(tau2-tau1),zeros(1,len1*m+nz*m-1)]);
    C_circ2=toeplitz([alpha2*beta2,alpha1*beta2+alpha2*beta1,alpha1*beta1,zeros(1,len1*m+nz*m-3)],[alpha2*beta2,zeros(1,len1*m+nz*m-1)]);

    C2_inv=inv(C_circ2);
    C_circ=C2_inv*C_circ1/T;
    %Subselecting rows to capture the downsampling operation in
    %spike and observed calcium
    C_circ_down=C_circ(1:m/m1:end,:);



    jj_0=0;
    jj_L=(length(F)/len1)-1;
    l1=length(jj_0:jj_L);

    %converting data into blocks for par-for processing
    t_est_=zeros(m*len1,l1);
    y_est_MR=zeros(len1,l1);
    spike_est_MR2=zeros(m*len1,l1);
    F1=double(reshape(F,len1,l1));


    t_frame11=reshape(t_frame,len1,l1);
    t_frame1=t_frame11-t_frame11(1,:);

    %Parallel procefssing of windowed data
    parfor jj=jj_0:jj_L
        y=F1(:,jj+1);
        %Interpolate to intermediate samplingrate
        y_up = interp1(t_frame1(:,jj+1),y,0:m*T/m1:len1*m*T-m*T/m1,'pchip','extrap')';
        %Accounting for non-zero initial condition for AR-model
        if jj==jj_0
            y_prev=y_up(1)*ones(nz*m1,1);
            y = [y_prev;y_up];
        else
            y_prev = interp1(t_frame1(1:nz+1,jj+1),[F1(end-nz+1:end,jj);y_up(1)],0:m*T/m1:nz*m*T,'pchip','extrap')';
            y_prev=y_prev(1:end-1);
            y=[y_prev;y_up];
        end


        t_est_(:,jj+1)=(0:T:len1*m*T-T);
        y_temp_prev=-C_circ2(2,1)*[y_prev(end);zeros(size(C_circ,1)-1,1)]-C_circ2(3,1)*[y_prev(end-1);y_prev(end);zeros(size(C_circ1,1)-2,1)];
        y_sub=C2_inv*y_temp_prev;
        %Non-convex FOCUSS to solve the spike deconvoltuion problem
        X_est=MFOCUSS(C_circ_down,y-y_sub(1:m/m1:end),l_f,'p',p,'prune_gamma',p_g,'MAX_ITERS', 5000);

        %Fit of dFoF at higher rate
        y1=C_circ*X_est;
        %Fit of the observed calcium signal
        y_est=y1(1:m:end);
        y_add=y_sub(1:m:end);
        y_est_MR(:,jj+1)=y_est(nz+1:end)+y_add(nz+1:end);
        spike_est_MR2(:,jj+1)=X_est(nz*m+1:end); 

    end
    t_est_=t_est_+t_frame11(1,:);
    t_est_=t_est_(:);
    spike_est_MR2=spike_est_MR2(:);
    y_est_MR=y_est_MR(:);
    spike_est_MR2_supp=zeros(size(spike_est_MR2));
    spike_est_MR2_supp(spike_est_MR2>2e-2)=1;
    %Vectorize spikes
    
end