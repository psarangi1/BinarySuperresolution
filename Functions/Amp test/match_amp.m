function [amp,amp_alter]=match_amp(y,alpha,m)
k=size(y,2);
vec_alpha=alpha.^(m-1:-1:0);
sol=ones(m,2^m);
amp=zeros(2^m,k);


for j = 0:2^m-1 % Tries all possibilities
       sol(:,j+1) = de2bi(j,m);
       for kk=1:k
        amp(j+1,kk)=(y(kk))/(vec_alpha*sol(:,j+1));
       end
       %amp_alter(j+1)=(y_2-alpha^m*y_1)/(vec_alpha*sol(:,j+1));
end

end