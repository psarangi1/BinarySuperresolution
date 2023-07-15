function [spike_est_MR_supp,t_MR_est,y_est_MR]=Pruning_algo(y_up,amp_algo,t_est_,tau)
len1=64;
jj_L=floor(length(y_up)/len1);
p=0.8;
l_f=1e-2;
p_g=1e-4;
pad_len=15;
trise=15e-3;
[spike_est_MR,spike_est_MR_supp,y_est_MR2,t_est_2]=MR_spike_est_rise(y_up(1:jj_L*len1).',trise,tau,2,1,l_f,p_g,p,pad_len,len1,t_est_((1:jj_L*len1)));

y_est_MR=y_est_MR2(:);
t_MR=t_est_2(:);

spike_est_MR_supp=zeros(size(spike_est_MR));
spike_est_MR_supp(spike_est_MR>=amp_algo-0.1*amp_algo)=1;
spike_est_MR_supp(spike_est_MR<amp_algo-0.1*amp_algo)=0;

t_MR_est=t_MR (logical(spike_est_MR_supp));
end