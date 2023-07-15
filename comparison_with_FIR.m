%%Script Compares the sample complexity for FIR vs IIR Filter
clear all;
close all;
addpath('Functions/')
%% Model parameters
N=100; %Length of binary vector
amp=1; %Binary spike amplitude   
p=0.5; %Bernoulli probability
%Filter Parameter - AR(1)
alpha=0.5; %Set with this value to ensure FIR recovery is possible using greedy algorithm (Can choose alpha<=0.5)
mc_iter=100; %Total number of monte carlo iterations
mmax=10; %Maximum undersampling factor (1-mmax)

%% Initialize the result variables with 0
avg_err_ham_m_ar=zeros(mmax,1);
avg_err_ham_m_fir=zeros(mmax,1);

err_rate_bin_iir=zeros(mmax,mc_iter,2);
err_rate_bin_fir=zeros(mmax,mc_iter,2);

%% FIR Filter Lengths
l_vec=[3,5,7];
for flen=1:length(l_vec)
    L=l_vec(flen); %Set filter length

    %%Loop for sweeping different downsampling facotrs
    for m=1:mmax

        %Construct suitable filtering matrix for generating the measurement
        p1=length(alpha)+1;
        G_alpha=toeplitz([1,-alpha,zeros(1,N-p1)],[1,zeros(1,N-1)]);
        G1=inv(G_alpha);
    
        %Initialize result variables
        err_ar=zeros(mc_iter,1);
        err_fir=zeros(mc_iter,1);

        %FIR filter (obtained by truncation of IIR AR(1) filter)
        h=alpha.^(0:L-1);

        
        %Prep binary search (executed only once for a given set of model
        %parameters):
        [c_fit,sol]=binary_prep(alpha,m,amp);

        %% Monte Carlo for each parameter
        parfor mc=1:mc_iter

            %%Generate Synthetic measurements
            
            %Ground truth spike generation 
            x=binornd(1,p,[N,1]);
            x(1)=0;

            %Output of IIR filter
            y_ar=G1*x;
            z=y_ar(1:m:end);
            
            %Decoding of Spikes (for IIR)
            x_est=SpikeDecodeAR(z,alpha,c_fit,sol,amp,m);
        
            M=length(z);
            
            TT=1+(M-1)*m;
            
            %Filtering matrix for FIR (Toeplitz)
            H=toeplitz([h(1);zeros(N-1,1)],[h';zeros(N-L,1)]);
            %Output of FIR filter
            y=H*x;
            y_down=y(1:m:end); %Uniform Downsampling

            %FIR Decoding Based on Greedy beta expansion (from ICASSP 2021) -
            spk_est=decode_spike_fir(y_down,h,0,m,N,p);
            
         
            %Compute Hamming Distance Error
            err_ar(mc)=norm(x_est-x(1:length(x_est)),1);
            err_fir(mc)=norm(spk_est(1:length(x_est))-x(1:length(x_est)),1);
    
               
        end
        avg_err_ham_m_ar(m)=mean(err_ar);
        avg_err_ham_m_fir(m)=mean(err_fir);
    
    end
    save(strcat('Simulation Data/FIR_vs_IIR_p',num2str(100*p),'_alpha_',num2str(100*alpha),'_flen_',num2str(L)))
end


%%Plot Generation

load('Simulation Data/FIR_vs_IIR_p50_alpha_50_flen_3')

M_vec=(1:10);
s=floor(N*p);
h1=figure;
plot(M_vec,avg_err_ham_m_ar/N,'LineWidth',2,'Marker','*','MarkerSize',10)
hold on,
plot(M_vec,avg_err_ham_m_fir/N,'LineWidth',2,'Marker','x','MarkerSize',10)
load('Simulation Data/FIR_vs_IIR_p50_alpha_50_flen_5')
plot(M_vec,avg_err_ham_m_fir/N,'LineWidth',2,'Marker','diamond','MarkerSize',10)
load('Simulation Data/FIR_vs_IIR_p50_alpha_50_flen_7')
plot(M_vec,avg_err_ham_m_fir/N,'LineWidth',2,'Marker','square','MarkerSize',10)

set(h1,'Units','Inches');
set(gca,'FontSize', 10)
pos = get(h1,'Position');
xlabel('Undersampling factor (D)','FontSize',14)
ylabel('Average Hamming Distance','FontSize',14)
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])


y_l=linspace(0.0,0.35,100);
plot(3*ones(1,100),y_l,'LineWidth',2,'LineStyle','--','HandleVisibility','off')
plot(5*ones(1,100),y_l,'LineWidth',2,'LineStyle','--','HandleVisibility','off')
plot(7*ones(1,100),y_l,'LineWidth',2,'LineStyle','--','HandleVisibility','off')
legend({'AR(1), \alpha=0.5', 'FIR (r=3)','FIR (r=5)','FIR (r=7)'},'FontSize',12,Location=[0.2,0.4,0.1,0.1])

text(3-0.5,0.25,'D=3','Color','red','FontSize',12)
text(5-0.5,0.25,'D=5','Color','red','FontSize',12)
text(7-0.5,0.25,'D=7','Color','red','FontSize',12)
grid minor
title('p=0.5, s=500','FontSize',15)
pbaspect([5,2,1])