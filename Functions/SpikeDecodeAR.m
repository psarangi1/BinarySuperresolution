function x_est=SpikeDecodeAR(z_n,alpha,c_fit,sol,amp,m)
%%Nearest Neighbor Spike Decoding algorithm
%Inputs:
%z_n: Low-rate noisy measurements
%alpha: AR(1) parameter
%amp: spike amplitude
%m: Downsampling factor
%Binary Search Offline parameters: c_fit (measurement Theta_alpha) and sol (spikes)
%Output: Decoded Binary spikes
%Decode the first sample
if abs(z_n(1)-amp)>z_n(1)
    x_est=0; 
else
    x_est=1;
end

%Construct the one-sample sub-problem 
w=z_n(2:end)-alpha^m*z_n(1:end-1);

%%Begin Decoding
for jj=1:size(w,1)
    %Search within the set Theta_alpha: For Each low-rate observation 
    %we find the nearest neighbor in the vector c_fit
    [minValue,closestIndex] = min(abs(w(jj)-c_fit));
    x_est=[x_est;sol(:,closestIndex)];
end
x_est=full(x_est);
end