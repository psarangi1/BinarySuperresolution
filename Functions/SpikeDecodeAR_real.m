function [spike_est,y_est_MR,t_est_]=SpikeDecodeAR_real(F,t_frame,c_fit,sol,m,alpha)
%%SUPERB multirate model based binary detection
%%One time Precomputation at the beginning 
%F       : low-rate, denoised, baseline corrected signal
%t_frame : low-resolution grid for frame times
%c_fit   : Precomputed quntities Theta_alpha
%sol     : Precomputed quntities S_all
%m       : multi-rate factor
%alpha   : AR model parameter
%%Output:
%spike_est      : High rate spikes (grid representation)
%t_est_         : High resolution timing grid
%y_est_MR       : estimate of the low rate calcium

%Processing to keep track of spike timing instances
t_est_=zeros((length(t_frame)-1)*m+1,1);
t_est_(1)=t_frame(1);
for i=2:length(t_frame)
  T=(t_frame(i)-t_frame(i-1))/(m);
  t_est_((i-2)*m+2:(i-1)*m+1)=t_frame(i-1)+(T:T:T*m);
end
index=zeros(size(F,1),1);
y_est_MR=zeros(size(F,1),1);
spike_est=zeros(m*(size(F,1)-1)+1,1);
F_diff=F-alpha^m*[0;F(1:end-1)];
sum_sol=sum(sol);

%Parallel implementation for estimating disjoing blocks (Identical to
%SpikeDecodeAR)
parfor i=1:length(F)
    %Nearest neighbor matching by exploiting the precomputed 
    %parameters from the multirate model
    [val1,id1]=sort(abs(F_diff(i)-c_fit));
    index(i)=id1(1);
end
spike_est(1)=c_fit(index(1));

%Generating estimated calcium waveform
y_est_MR(1)=c_fit(index(1));
for i=2:size(F,1)
    spike_est((i-2)*m+2:(i-1)*m+1)=flipud(sol(:,index(i)));
    y_est_MR(i)=c_fit(index(i))+alpha^m*y_est_MR(i-1);
end

end