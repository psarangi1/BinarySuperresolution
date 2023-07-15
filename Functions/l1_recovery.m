function x_l1=l1_recovery(H,zz,TT,noise) 
    %Spikes estimation by solving l1 minimization problem with
    %non-negativity
    %(Implemented using CVX)
    %Inputs:
    %zz: low-rate noisy observation
    %H: measurement matrix (obtained by undersampling)
    %noise: upper bound on l2 norm ball
    %Output (x_l1): Estimated spikes
    N=size(H,2);
    cvx_begin quiet
        variable x_l1(N,1)
        minimize norm(x_l1,1)
        subject to
        norm(zz-H*x_l1)<=noise
        x_l1>=0
    cvx_end
end