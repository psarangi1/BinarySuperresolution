clear all;
close all;
%rng(5)
N=2048;
x = zeros(N,1); %my signal
alpha=0.95;
gamma_vec=[-alpha];
G=toeplitz([1,gamma_vec,zeros(1,N-length(gamma_vec)-1)],[1,zeros(1,N-1)]);
p=0.1;
amp=1;
x=amp*binornd(1,p*ones(N,1));
x(1)=amp;
G1=inv(G);
y=G\x;
sigma=1e-2;
m=8;
off_val=G\ones(N,1)/2;
y_noise=y+sigma*randn(N,1);
y_clean=y(1:m:end);
y_down=y_noise(1:m:end);
%y_down=y_noise(1:m:end)-off_val(1:m:end);

G_d=toeplitz([1,-gamma_vec^m,zeros(1,N/m-length(gamma_vec)-1)],[1,zeros(1,N/m-1)]);
G_d1=inv(G_d);
cvx_begin
variable z(N/m,1)
minimize norm(z,1)
subject to
norm(y_down-G_d1*z)<=1.2e-2*sqrt(N/m)
z>=0
cvx_end
y_down2=G_d1*z;
c0=[0;1];
%c0=[-1/2;1/2];
x_est=[];
sol=ones(m,2^m);
for j = 0:2^m-1 % Tries all possibilities
       sol(:,j+1) = amp*de2bi(j,m);
end
%sol(sol==0)=-1/2;
%sol(sol==1)=1/2;
c_fit=zeros(size(y_down));
%alpha_t=alpha-0*rand(1);
[vv,id_est]=sort(y_down2-alpha^m*[0;y_down2(1:end-1)],'descend');
id_est=id_est(1:30);
amp1=match_amp([y_down2(id_est).';y_down2(id_est-1).'],alpha,m);
%sigma=sigma+abs(alpha_t-alpha);
%length(intersect(amp1,amp2))
% amp_est=amp1(:,1);
% [v,id_b]=min(abs(amp1(:,1)-amp));
% amp_true_match=amp1(id_b,1);
% for i=2:size(amp1,2)
%     [v,id_b]=min(abs(amp1(:,i)-amp_true_match));
%     amp_true_match=(amp1(id_b,i)+amp_true_match)/2;
%     
%     diff_mat=abs(amp1(:,i)-amp_est.');
%     temp=[];
%     for j=1:size(diff_mat,2)
%          [match_amp,id_min]=min(diff_mat(:,j));
%          %id_b,id_min
%          if match_amp<=1e-2
%              %temp=[temp;(amp_est(j)+amp1(id_min,i))/2];
%              temp=[temp;amp_est(j)];
%          end
%         %          %226,151,124,100,100,83,83,65,45
%         %          if j==45
%         %              1
%         %          end
%          
%     end
%     amp_est=temp;
%     
%     %     [id_match1,id_match2]=find(diff_mat<4e-3);
%     %     amp_est=amp_est(unique(id_match2));
% end
amp_est=resolve_amp(amp1);
[c_fit,alpha,sol]=binary_prep(1,m,0:m/N:1-m/N,amp_est(1));
[spike_est,y_est_MR,t_est_]=MR_spike_est_binary(y_down2,0:m/N:1-m/N,c_fit,sol,m,alpha);
figure,
ax1=subplot(2,1,1);
stem(t_est_,spike_est),
ax2=subplot(2,1,2);
stem(t_est_,x)
linkaxes([ax1,ax2],'x')