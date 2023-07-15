%%Script for Demonstrating Spike Reconstruction Performance 
clear all;
close all;
addpath(genpath('Functions'))

seed=1;
rng(seed)
%Generative Model parameters

%Length of Binary Spiking Signal
N=51;
%Spiking Probability
p=0.35;
%AR(1) filter parameter
alpha=0.9;
p1=length(alpha)+1;
G_alpha=toeplitz([1,-alpha,zeros(1,N-p1)],[1,zeros(1,N-1)]);
G1=inv(G_alpha);
amp=1; %Spike amplitude

x=amp*binornd(1,p,[N-1,1]);
x(1)=0;
x=[x;0];

%Downsampling factor
m=10;
y=G1*x;
z=y(1:m:end);
M=length(z);
noise=0.0; %Noise level

%Generating noisy observations (Additive Gaussian)
z_n=z+noise*randn(length(z),1);
%Time vector 
t=linspace(0,N-1,N); %hi-res grid
t_low=t(1:m:end); %low-res grid

DD=eye(N);
D=DD(1:m:end,:);
       
        
TT=1+(M-1)*m;
   
H=D*G1(:,1:TT);

x_l1=l1_recovery(H,z_n,TT,noise*sqrt(length(z_n))); 
x_l1_box=l1_recovery_box(H,z_n,TT,noise*sqrt(length(z_n)),amp); 

amp=1;

[c_fit,sol]=binary_prep(alpha,m,amp);
x_est=SpikeDecodeAR(z_n,alpha,c_fit,sol,amp,m);
L=length(x_est);    

id11=find(x>0);
on_low=find(mod((id11-1),m)==0);
in_bet=setdiff(id11,id11(on_low));


[status, msg, msgID] = mkdir('Waveform Plots');
%Plot Paramaters
plot_d=0;
plot_up=N;
h1=figure;
stem(t(in_bet),x(in_bet),'r','LineWidth',2,'Marker', 'none')
hold on
stem(t,zeros(length(t),1),'MarkerSize',5)
if ~isempty(on_low)
    stem(t(id11(on_low)),x(id11(on_low)),'Color',[0.8500 0.3250 0.0980],'LineWidth',2,'Marker', 'none')
end
set(h1,'Units','Inches');
set(gca,'FontSize', 12,'box','off','ytick',[])
pos = get(h1,'Position');
xlabel('Time (s)','FontSize',8)
xlim([plot_d,plot_up])
%ylabel('Amplitude','FontSize',20)
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
title('Ground Truth Spike')
pbaspect([10 1 1])
%pbaspect([4 1 1])
set(gca,'XTick',t_low);
set(gca,'XTickLabel',t_low);
file1=strcat('Waveform Plots/Spikes_gt_bin_noise_',num2str(100*noise),'_m_',num2str(m),'_seed_',num2str(seed),'.eps');
print(h1,file1,'-depsc','-r0')
%subplot(2,1,2)

h1=figure;
stem(t(1:L),x_est(1:L),'b','LineWidth',2,'Marker', 'none')
hold on
stem(t,zeros(length(t),1),'MarkerSize',5)
set(h1,'Units','Inches');
set(gca,'FontSize', 10,'box','off','ytick',[])
pos = get(h1,'Position');
xlabel('Time (s)','FontSize',8)
xlim([plot_d,plot_up])
%ylabel('Amplitude','FontSize',20)
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
title('Spike Estimate (Binary')
pbaspect([10 1 1])
set(gca,'XTick',t_low);
set(gca,'XTickLabel',t_low);
%pbaspect([4 1 1])
file2=strcat('Waveform Plots/Spikes_est_bin_',num2str(100*noise),'_m_',num2str(m),'_seed_',num2str(seed),'.eps');
print(h1,file2,'-depsc','-r0')


h1=figure;
stem(t(1:TT),x_l1,'Color',[0.4940 0.1840 0.5560],'LineWidth',2,'Marker', 'none')
hold on
stem(t_low,zeros(length(t_low),1),'MarkerSize',1)
set(h1,'Units','Inches');
set(gca,'FontSize', 12,'box','off','ytick',[])
pos = get(h1,'Position');
xlabel('Time (s)','FontSize',8)
xlim([plot_d,plot_up])
%ylabel('Amplitude','FontSize',20)
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
title('Spike Estimate (l_1)')
pbaspect([10 1 1])
%pbaspect([4 1 1])
set(gca,'XTick',t_low);
set(gca,'XTickLabel',t_low);
file3=strcat('Waveform Plots/Spikes_est_l1_',num2str(100*noise),'_m_',num2str(m),'_seed_',num2str(seed),'.eps');
print(h1,file3,'-depsc','-r0')


h1=figure;
stem(t(1:TT),x_l1_box,'Color',[0.4940 0.1840 0.5560],'LineWidth',2,'Marker', 'none')
hold on
stem(t_low,zeros(length(t_low),1),'MarkerSize',1)
set(h1,'Units','Inches');
set(gca,'FontSize', 12,'box','off','ytick',[])
pos = get(h1,'Position');
xlabel('Time (s)','FontSize',8)
xlim([plot_d,plot_up])
%ylabel('Amplitude','FontSize',20)
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
title('Spike Estimate (Box-l_1)')
pbaspect([10 1 1])
%pbaspect([4 1 1])
set(gca,'XTick',t_low);
set(gca,'XTickLabel',t_low);
file3_2=strcat('Waveform Plots/Spikes_est_l1_box_',num2str(100*noise),'_m_',num2str(m),'_seed_',num2str(seed),'.eps');
print(h1,file3_2,'-depsc','-r0')



h1=figure;
plot(t,y,'b','LineWidth',2)
hold on
stem(t,zeros(length(t),1),'MarkerSize',5)
xlabel('Time (s)','FontSize',8)
%pbaspect([10 1 1])
pbaspect([10 1 1])
set(h1,'Units','Inches');
set(gca,'FontSize', 12,'box','off','ytick',[])
pos = get(h1,'Position');
%ylabel('Amplitude','FontSize',20)
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
title('Calcium (high-rate)')
xlim([plot_d,plot_up])
set(gca,'XTick',t(1:2:end));
set(gca,'XTickLabel',t(1:2:end));
file3=strcat('Waveform Plots/Cal_Waveform_hi',num2str(100*noise),'_m_',num2str(m),'_seed_',num2str(seed),'.eps');
print(h1,file3,'-depsc','-r0')

h1=figure;
stem(t_low,z_n,'k-','LineWidth',2,'MarkerSize',5)
pbaspect([10 1 1])
set(h1,'Units','Inches');
set(gca,'FontSize', 12,'box','off','ytick',[])
pos = get(h1,'Position');
xlabel('Time (s)','FontSize',8)
%ylabel('Amplitude','FontSize',20)
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
title('Calcium (low-rate)')
xlim([plot_d,plot_up])
set(gca,'XTick',t_low);
set(gca,'XTickLabel',t_low);
file4=strcat('Waveform Plots/Cal_Waveform_lo',num2str(100*noise),'_m_',num2str(m),'_seed_',num2str(seed),'.eps');
print(h1,file4,'-depsc','-r0')

