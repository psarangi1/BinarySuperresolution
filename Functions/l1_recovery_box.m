function x_l1_box=l1_recovery_box(H,zz,TT,noise,amp) 
    %Spikes estimation by solving l1 minimization problem with box constraint
    %(Implemented using CVX)
    %Inputs:
    %zz: low-rate noisy observation
    %H: measurement matrix (obtained by undersampling)
    %noise: upper bound on l2 norm ball
    %amp: spike amplitude for box-constraint
    %Output (x_l1): Estimated spikes
    N=size(H,2);
    cvx_begin quiet
        variable x_l1_box(N,1)
        minimize norm(x_l1_box,1)
        subject to
        norm(zz-H*x_l1_box)<=1.02*noise
        x_l1_box>=0
        x_l1_box<=amp
    cvx_end
end