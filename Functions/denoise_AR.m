function ca=denoise_AR(y_corrected,samp_rate,tau,type)

if type== 1
    ops.imageRate = samp_rate;
    ops.sensorTau = tau;
    ops.deconvType = 'OASIS';%'L0' or 'OASIS' 
    [sp, ca, coefs, B, sd, ops, baselines,F1_p] = wrapperDECONV(ops,y_corrected);
else 
    lambda=0;
     g=exp(-1/(samp_rate*tau));
    [ca, ~, ~] = deconvolveCa(y_corrected, 'ar1',g, 'foopsi', 'lambda', lambda,'optimize_pars'); 

end

end