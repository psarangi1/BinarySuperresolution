function [c_fit,sol]=binary_prep(alpha,m,amp)
%%Preparation of the set Theta_alpha and S_all
%%One time Precomputation at the beginning 
%alpha: AR model parameter (at high rate)
%m    : Undersampling factor (1<m<14)
%amp  : spike amplitude (A)
vec_alpha=alpha.^(m-1:-1:0);
sol=sparse(zeros(m,2^m));
c_fit=zeros(2^m,1);
for j = 0:2^m-1 % Tries all possibilities
       sol(:,j+1) = sparse(amp*de2bi(j,m));
       c_fit(j+1)=vec_alpha*sol(:,j+1);
end
end