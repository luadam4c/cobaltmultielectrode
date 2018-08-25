%%%%%figure1%%%%
% typical seizure
% june_09_09_seizure_12 cut off first 12s
data=data(12*25000:length(data),:);
plot((1:length(data))/25000,data(:,1)+0.08);
set(gca,'FontSize',20)
xlabel('Time(second)','FontSize',20)
ylabel('Voltage(V)','FontSize',20)
%% june_09_09_seizure_12

%%% freqency plot%%%%%
Fs=500;
% t = 1:(1/Fs)*250000:(length(data)); 
t = 6.3*25000:(1/Fs)*250000:7.3*25000; 
% refer to http://www.mathworks.com/support/tech-notes/1700/1702.html
x=data(t,1);
% Use next highest power of 2 greater than or equal to length(x) to calculate FFT.
nfft= 2^(nextpow2(length(x))); 
% Take fft, padding with zeros so that length(fftx) is equal to nfft 
fftx = fft(x,nfft);
% Calculate the numberof unique points
NumUniquePts = ceil((nfft+1)/2); 
% FFT is symmetric, throw away second half 
fftx = fftx(1:NumUniquePts); 
% Take the magnitude of fft of x and scale the fft so that it is not a function of the length of x
mx = abs(fftx)/length(x); 
% Take the square of the magnitude of fft of x. 
mx = mx.^2; 
% Since we dropped half the FFT, we multiply mx by 2 to keep the same energy.
% The DC component and Nyquist component, if it exists, are unique and
% should not be multiplied by 2.
if rem(nfft, 2) % odd nfft excludes Nyquist point
  mx(2:end) = mx(2:end)*2;
else
  mx(2:end -1) = mx(2:end -1)*2;
end
% This is an evenly spaced frequency vector with NumUniquePts points. 
f = (0:NumUniquePts-1)*Fs/nfft; 
% Generate the plot, title and labels. 
plot(f,mx); 
xlabel('Frequency (Hz)'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

peak_threshold=-0.05;
min_num_pts=25;
spikes_position=find_spikes(data(:,1), peak_threshold, min_num_pts)+260000;
diff_spikes=diff(spikes_position);
period=[];
period=[period;diff_spikes];

sort_freq=sort(1./period)*25000;
binWidth=5;
binCtrs=0:binWidth:max(sort_freq);
counts=hist(sort_freq,binCtrs);
bar(binCtrs, counts/sum(counts));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=data(5*25000:10*25000,1);
N = length(x); %% number of points
T= length(x)/25000; %% time
f = data(:,1); %%define function, 10 Hz sine wave
p = abs(fft(f))/(N/2); %% absolute value of the fft
p = p(1:N/2).^2; %% take the power of positve freq. half
freq = [0:N/2-1]/T; %% find the corresponding frequency in Hz
semilogy(freq,p); %% plot on semilog scale
%axis([0 20 0 1]); %% zoom in

%%%%%%%%%%%%%%%%%%%%%





%% distribution of seizure length june-aug 09
seizure_length=[ 20
14.6
18.1
20.2
11.5
22
21.3
21.6
22.3
22.4
25.8
21.3];
% seizure_length_1=[
% 22.8
% 18.87
% 20.1
% 21.6
% 21.4
% 20.4
% 21.1
% 21.1
% 20.8
% 21.5
% 20.6
% 21.9
% ;];
% 
seizure_length_2=[24.5
25.28
15.82
17.22
13.84
15.29
18.82
23.01
28.02
26.7
35.6
35.3
];
seizure_length_3=[
5.7801
6.7245
7.523
6.89
8.6394
7.7218
10.033
9.1346
8.53
7.62
11.332
13.184];

% seizure_length_4=[
% 6.4568
% 7.666
% 8.4989
% 8.1064
% 8.77 
% 9.05 
% 7.9
% 10.77
% 9.8621
% 10.805
% 10.635
% 11.331
% ];
% 
% seizure_length_5=[
% 18.198
% 19.679
% 20.286
% 20.554
% 23.132
% 16.064
% 21.607
% 28.964
% 25.869
% 25.24
% 20.564
% 35.9];
% 
x=1:12;
plot(x,seizure_length,'*')
hold on 
% plot(x,seizure_length_1,'*')
plot(x,seizure_length_2,'*')
plot(x,seizure_length_3,'*')
% plot(x,seizure_length_4,'*')
% plot(x,seizure_length_5,'*')
%% july_16_seizure_4
% refer to seizure_2
% four channel

Hd=bpf1;
data(:,1)=filter(Hd, data(:,1));
data(:,2)=filter(Hd, data(:,2));
data(:,3)=filter(Hd, data(:,3));
data(:,4)=filter(Hd, data(:,4));

peak_threshold=-0.05;%threshold to pick the spike
period=[];
min_num_pts=5;
spikes_position_1=find_spikes(-abs(data(:,1)), peak_threshold,min_num_pts);
spikes_position_2=find_spikes(-abs(data(:,2)), peak_threshold,min_num_pts);
spikes_position_3=find_spikes(-abs(data(:,3)), peak_threshold,min_num_pts);
spikes_position_4=find_spikes(-abs(data(:,4)), peak_threshold,min_num_pts);

diff_spikes_1=diff(spikes_position_1);
diff_spikes_2=diff(spikes_position_2);
diff_spikes_3=diff(spikes_position_3);
diff_spikes_4=diff(spikes_position_4);

window=[];
for i=1:length(diff_spikes_4)
if diff_spikes_4(i)>3000
    window=[window,spikes_position_4(i)];
end
end
window=window';


coeff=[];
 for i=1:length(window)-1
     data_window=[];
     data_window=[(window(i)+5000:window(i+1)+5000)' data(window(i)+5000:window(i+1)+5000,:)];
     coeff1=corrcoef(data_window(:,2),data_window(:,3)); % ch1 ch2
     coeff2=corrcoef(data_window(:,2),data_window(:,4)); % ch1 ch3
     coeff3=corrcoef(data_window(:,2),data_window(:,5)); % ch1 ch4
     coeff4=corrcoef(data_window(:,3),data_window(:,4)); % ch2 ch3
     coeff5=corrcoef(data_window(:,3),data_window(:,5)); % ch2 ch4
     coeff6=corrcoef(data_window(:,4),data_window(:,5)); % ch3 ch4
     coeff=[coeff; coeff1(1,2) coeff2(1,2) coeff3(1,2) coeff4(1,2) coeff5(1,2) coeff6(1,2)]; 
 end 
 
coeff_spikes=[];
 for i=1:length(window)-1
     data_window=[];
     % get bursting window
     data_window=[(window(i)+5000:window(i+1)+5000)' data(window(i)+5000:window(i+1)+5000,:)];
     % set threshold for spike finding
     peak_threshold=-0.05;%threshold to pick the spike
     period=[];
     min_num_pts=5;
     % use a channel to pick spikes
     spikes_position=find_spikes(-abs(data_window(:,5)), peak_threshold,min_num_pts);
     spikes_position=data_window(spikes_position,1);
     for j=1:length(spikes_position)
     spike_window=[(spikes_position(j)-50:spikes_position(j)+50)' data((spikes_position(j)-50:spikes_position(j)+50),1:4)];
     coeff1=corrcoef(spike_window(:,2),spike_window(:,3)); % ch1 ch2
     coeff2=corrcoef(spike_window(:,2),spike_window(:,4)); % ch1 ch3
     coeff3=corrcoef(spike_window(:,2),spike_window(:,5)); % ch1 ch4
     coeff4=corrcoef(spike_window(:,3),spike_window(:,4)); % ch2 ch3
     coeff5=corrcoef(spike_window(:,3),spike_window(:,5)); % ch2 ch4
     coeff6=corrcoef(spike_window(:,4),spike_window(:,3)); % ch3 ch4
     coeff_spikes=[coeff_spikes; j  coeff1(1,2) coeff2(1,2) coeff3(1,2) coeff4(1,2) coeff5(1,2) coeff6(1,2)]; 
     end 
 end 
 



%%%%%%%

A=mean(coeff(:,1));
B=mean(coeff(:,2));
C=mean(coeff(:,3));
D=mean(coeff(:,4));
E=mean(coeff(:,5));
F=mean(coeff(:,6));
coeff_demo=[C E F 1; B D 1 F; A 1 D E ; 1 A B C];
imagesc(coeff_demo)
colormap jet
axis image


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot( (1:length(data))/25000,data(:,1)+0.5,'k')
hold on
plot( (1:length(data))/25000,data(:,2),'b')
plot( (1:length(data))/25000,data(:,3)-0.5,'r')
plot( (1:length(data))/25000,data(:,4)-1,'color',[0 0.5 0])
h = legend('ch1','ch2','ch3','ch4',2);
set(h,'Interpreter','none','FontSize',20)
set(gca,'FontSize',20)
xlabel('Time(second)','FontSize',20)
ylabel('Voltage(mV)','FontSize',20)

% ch2-ch4
plot( (1:length(data))/25000,data(:,2)+0.5,'b')
hold on
plot( (1:length(data))/25000,data(:,4),'color',[0 0.5 0])
h = legend('ch2','ch4',2);
set(h,'Interpreter','none','FontSize',20)
set(gca,'FontSize',20)
xlabel('Time(second)','FontSize',20)
ylabel('Voltage(mV)','FontSize',20)

% coeff ch2-4

plot(window(2:length(window))/25000, coeff(:,5),'linewidth',4)
set(gca,'FontSize',20)
xlabel('Time(second)','FontSize',20)
ylabel('correlation coefficient','FontSize',20)

%%%%
binWidth=0.05;
subplot(1,3,1)
sort_coeff_spikes_1=sort(coeff_spikes((spikeburstno(5):spikeburstno(6)-1)',5));
binCtrs=0:binWidth:1;
counts=hist(sort_coeff_spikes_1,binCtrs);
bar(binCtrs, counts/sum(counts));
xlim([0 1.2]);
ylim([0 0.5]);
set(gca,'FontSize',20)

subplot(1,3,2)
sort_coeff_spikes_2=sort(coeff_spikes((spikeburstno(18):spikeburstno(19)-1)',5));
binCtrs=0:binWidth:1;
counts=hist(sort_coeff_spikes_2,binCtrs);
bar(binCtrs, counts/sum(counts));
xlim([0 1.2]);
ylim([0 0.5]);
set(gca,'FontSize',20)
xlabel('correlation coefficient','FontSize',20)

subplot(1,3,3)
sort_coeff_spikes_3=sort(coeff_spikes((spikeburstno(20):spikeburstno(21)-1)',5));
binCtrs=0:binWidth:1;
counts=hist(sort_coeff_spikes_3,binCtrs);
bar(binCtrs, counts/sum(counts));
xlim([0 1.2]);
ylim([0 0.5]);
set(gca,'FontSize',20)
% 

%% july_16_10 seizure 2
% theta rythme
data(:,5)=data(:,4);

Hd=bpf1;
data(:,1)=filter(Hd, data(:,1));
data(:,2)=filter(Hd, data(:,2));
data(:,3)=filter(Hd, data(:,3));
data(:,4)=filter(Hd, data(:,4));


theta_data=xin_seizure_2_MultiCH_Data(252500:length(xin_seizure_2_MultiCH_Data),4);
% theta_data=xin_seizure_10_MultiCH_Data(10.2*25000:50*25000,4);
Hd=mytheta;
theta_data=filter(Hd,theta_data); %butterworth

plot( (1:length(data))/25000,data(:,5)+2)
hold on
plot( (1:length(data))/25000,theta_data,'k')

h = legend('EEG','Theta Oscillation',2);
set(h,'Interpreter','none','FontSize',20)
set(gca,'FontSize',20)
xlabel('Time(second)','FontSize',20)
ylabel('Voltage (mV)','FontSize',20)

% find spikes
peak_threshold=-0.05;%threshold to pick the spike
period=[];
min_num_pts=5;
spikes_position_1=find_spikes(-abs(data(:,1)), peak_threshold,min_num_pts);
spikes_position_2=find_spikes(-abs(data(:,2)), peak_threshold,min_num_pts);
spikes_position_3=find_spikes(-abs(data(:,3)), peak_threshold,min_num_pts);
spikes_position_4=find_spikes(-abs(data(:,4)), peak_threshold,min_num_pts);

 
 plot( (1:length(data))/25000,data(:,1)+0.5)
hold on
plot( (1:length(data))/25000,data(:,2))
plot( (1:length(data))/25000,data(:,3)-0.5)
plot( (1:length(data))/25000,data(:,4)-1)
ylim([-2 2])
set(gca,'FontSize',20)
xlabel('Time(second)','FontSize',20)
ylabel('Amplitude (mV)','FontSize',20)

plot(spikes_position_1/25000,data(spikes_position_1,1)+0.5,'*','MarkerEdgeColor','k')
hold on
plot(spikes_position_2/25000,data(spikes_position_2,2),'*','MarkerEdgeColor','r')
plot(spikes_position_3/25000,data(spikes_position_3,3)-0.5,'*','MarkerEdgeColor','g')
plot(spikes_position_4/25000,data(spikes_position_4,4)-1,'*','MarkerEdgeColor','m')

% raster plot
plot(spikes_position_1/25000,-0.5,'.','MarkerEdgeColor','k')
hold on
plot(spikes_position_2/25000,-1,'.','MarkerEdgeColor','r')
plot(spikes_position_3/25000,-1.5,'.','MarkerEdgeColor','g')
plot(spikes_position_4/25000,-2,'.','MarkerEdgeColor','m')


%spikes vs theta
plot( (1:length(data))/25000,theta_data)
hold on

plot(spikes_position_1/25000,theta_data(spikes_position_1,1),'*','MarkerEdgeColor','k')
hold on
plot(spikes_position_2/25000,theta_data(spikes_position_2,1),'*','MarkerEdgeColor','r')
plot(spikes_position_3/25000,theta_data(spikes_position_3,1),'*','MarkerEdgeColor','g')
plot(spikes_position_4/25000,theta_data(spikes_position_4,1),'*','MarkerEdgeColor','m')

% plot(spikes_position_1/25000,0.4,'*','MarkerEdgeColor','k')
% hold on
% plot(spikes_position_2/25000,0.4,'*','MarkerEdgeColor','r')
% plot(spikes_position_3/25000,0.4,'*','MarkerEdgeColor','g')
% plot(spikes_position_4/25000,0.4,'*','MarkerEdgeColor','m')

h = legend('Theta Oscillation','ch1_spike','ch2_spike','ch3_spike','ch4_spike',2);
set(h,'Interpreter','none','FontSize',20)
set(gca,'FontSize',20)
xlabel('Time(second)','FontSize',20)
ylabel('Amplitude (mV)','FontSize',20)

% find spike for theta
plot( (1:length(data))/25000,theta_data)
hold on
peak_threshold=-0.05;%threshold to pick the spike
period=[];
min_num_pts=1;
theta_spikes_position=find_spikes((theta_data(:,1)), peak_threshold,min_num_pts);
plot(theta_spikes_position/25000,theta_data(theta_spikes_position,1),'*','MarkerEdgeColor','r')


%define interval
phase_data=[];
phase_data(:,1)=(1:length(data))';% 1st time
%2nd give each time point a phase value
for k=1:theta_spikes_position(1)
phase_data(k,2)=2*pi/theta_spikes_position(1)*k;
end
for i=1:(length(theta_spikes_position)-1)
    for j=theta_spikes_position(i):theta_spikes_position(i+1)
    phase_data(j,2)=2*pi/(theta_spikes_position(i+1)-theta_spikes_position(i))*(j-theta_spikes_position(i));
    end
end
for z=theta_spikes_position(length(theta_spikes_position)):length(data)
    phase_data(z,2)=2*pi/(length(data)-theta_spikes_position(length(theta_spikes_position)))*(z-theta_spikes_position(length(theta_spikes_position)));
end

plot( (1:length(data))/25000,theta_data)
hold on
plot(phase_data(:,1)/25000,phase_data(:,2),'k')
hold on
plot(theta_spikes_position/25000,theta_data(theta_spikes_position,1),'*','MarkerEdgeColor','r')
h = legend('Theta oscillation','Theta Phase',2);
set(h,'Interpreter','none','FontSize',20)
set(gca,'FontSize',20)
xlabel('Time(second)','FontSize',20)


% give spike-time a phase value
phase_ch1=phase_data(spikes_position_1,2);
phase_ch2=phase_data(spikes_position_2,2);
phase_ch3=phase_data(spikes_position_3,2);
phase_ch4=phase_data(spikes_position_4,2);

plot(phase_data(:,1)/25000,phase_data(:,2),'k')
hold on;plot(spikes_position_1/25000,phase_ch1,'*','MarkerEdgeColor','k')
hold on; plot(spikes_position_2/25000,phase_ch2,'*','MarkerEdgeColor','r')
hold on; plot(spikes_position_3/25000,phase_ch3,'*','MarkerEdgeColor','g')
hold on; plot(spikes_position_4/25000,phase_ch4,'*','MarkerEdgeColor','m')

h = legend('Theta Phase','ch1_spike_phase','ch2_spike_phase','ch3_spike_phase','ch4_spike_phase',2);
set(h,'Interpreter','none','FontSize',20)
set(gca,'FontSize',20)
xlabel('Time(second)','FontSize',20)
ylabel('phase','FontSize',20)

%%%%%%%%%%%%%%%
binWidth=0.2;

subplot(2,2,1)
sort_phase_ch1=sort(phase_ch1);
binCtrs=0:binWidth:2*pi;
counts1=hist(sort_phase_ch1,binCtrs);
bar(binCtrs, counts1/sum(counts1));
xlim([0 2*pi]);
ylim([0 0.15]);
set(gca,'FontSize',20)
set(gca,'XTick',[0 pi 2*pi])


subplot(2,2,2)
sort_phase_ch2=sort(phase_ch2);
binCtrs=0:binWidth:2*pi;
counts2=hist(sort_phase_ch2,binCtrs);
bar(binCtrs, counts2/sum(counts2));
xlim([0 2*pi]);
ylim([0 0.15]);
set(gca,'FontSize',20)
set(gca,'XTick',[0 pi 2*pi])



subplot(2,2,3)
sort_phase_ch3=sort(phase_ch3);
binCtrs=0:binWidth:2*pi;
counts3=hist(sort_phase_ch3,binCtrs);
bar(binCtrs, counts3/sum(counts3));
xlim([0 2*pi]);
ylim([0 0.15]);
set(gca,'FontSize',20)
set(gca,'XTick',[0 pi 2*pi])

subplot(2,2,4)
sort_phase_ch4=sort(phase_ch4);
binCtrs=0:binWidth:2*pi;
counts4=hist(sort_phase_ch4,binCtrs);
bar(binCtrs, counts4/sum(counts4));
xlim([0 2*pi]);
ylim([0 0.15]);
set(gca,'FontSize',20)
set(gca,'XTick',[0 pi 2*pi])    

%%%%%%%%%
% phase difference

window_size=25;
new_spikes_position_1=[];
new_spikes_position_2=[];
new_spikes_position_3=[];
new_spikes_position_4=[];
for i=1:length(spikes_position_1)
    new_spikes_position2=spikes_position_2(spikes_position_2>=spikes_position_1(i)-window_size & spikes_position_2<=spikes_position_1(i)+window_size);
    new_spikes_position3=spikes_position_3(spikes_position_3>=spikes_position_1(i)-window_size& spikes_position_3<=spikes_position_1(i)+window_size);
    new_spikes_position4=spikes_position_4(spikes_position_4>=spikes_position_1(i)-window_size& spikes_position_4<=spikes_position_1(i)+window_size);
    if isempty(new_spikes_position2)|isempty(new_spikes_position3)|isempty(new_spikes_position4)
    else  
        if new_spikes_position2>=spikes_position_1(i)
            new_spikes_position2=min(new_spikes_position2);
        else
            new_spikes_position2=max(new_spikes_position2);
        end
        if new_spikes_position3>=spikes_position_1(i)
            new_spikes_position3=min(new_spikes_position3);
        else
            new_spikes_position3=max(new_spikes_position3);
        end
        if new_spikes_position4>=spikes_position_1(i)
            new_spikes_position4=min(new_spikes_position4);
        else
            new_spikes_position4=max(new_spikes_position4);
        end
    new_spikes_position_1=[new_spikes_position_1 spikes_position_1(i)];
    new_spikes_position_2=[new_spikes_position_2 new_spikes_position2];
    new_spikes_position_3=[new_spikes_position_3 new_spikes_position3];
    new_spikes_position_4=[new_spikes_position_4 new_spikes_position4];
    end
    clear new_spikes_position2;
    clear new_spikes_position3;
    clear new_spikes_position3;
end

plot((1:length(data))/25000,data(:,1)+1)
hold on
plot((1:length(data))/25000,data(:,2)+0.5)
plot((1:length(data))/25000,data(:,3))
plot((1:length(data))/25000,data(:,4)-0.5)

plot(new_spikes_position_1/25000,data(new_spikes_position_1,1)+1,'*','MarkerEdgeColor','k')
hold on
plot(new_spikes_position_2/25000,data(new_spikes_position_2,2)+0.5,'*','MarkerEdgeColor','r')
plot(new_spikes_position_3/25000,data(new_spikes_position_3,3),'*','MarkerEdgeColor','g')
plot(new_spikes_position_4/25000,data(new_spikes_position_4,4)-0.5,'*','MarkerEdgeColor','m')

%%%%% find theta phase difference
window_size=50;
new_spikes_position_12=[];
new_spikes_position_21=[];
new_spikes_position_13=[];
new_spikes_position_31=[];
new_spikes_position_14=[];
new_spikes_position_41=[];
new_spikes_position_23=[];
new_spikes_position_32=[];
new_spikes_position_24=[];
new_spikes_position_42=[];
new_spikes_position_34=[];
new_spikes_position_43=[];

for i=1:length(spikes_position_1)
    new_spikes_position2=spikes_position_2(spikes_position_2>=spikes_position_1(i)-window_size & spikes_position_2<=spikes_position_1(i)+window_size);
    new_spikes_position3=spikes_position_3(spikes_position_3>=spikes_position_1(i)-window_size& spikes_position_3<=spikes_position_1(i)+window_size);
    new_spikes_position4=spikes_position_4(spikes_position_4>=spikes_position_1(i)-window_size& spikes_position_4<=spikes_position_1(i)+window_size);
%     if isempty(new_spikes_position2)|isempty(new_spikes_position3)|isempty(new_spikes_position4)
    if isempty(new_spikes_position2)
    else  
        if new_spikes_position2>=spikes_position_1(i)
            new_spikes_position2=min(new_spikes_position2);
        else
            new_spikes_position2=max(new_spikes_position2);
        end
    new_spikes_position_12=[new_spikes_position_12 spikes_position_1(i)];
    new_spikes_position_21=[new_spikes_position_21 new_spikes_position2];
    end
    
   if isempty(new_spikes_position3)
   else
       if new_spikes_position3>=spikes_position_1(i)
            new_spikes_position3=min(new_spikes_position3);
         else
             new_spikes_position3=max(new_spikes_position3);
       end
        new_spikes_position_13=[new_spikes_position_13 spikes_position_1(i)];
        new_spikes_position_31=[new_spikes_position_31 new_spikes_position3];
   end
    
   if isempty(new_spikes_position4)
   else
        if new_spikes_position4>=spikes_position_1(i)
            new_spikes_position4=min(new_spikes_position4);
        else
            new_spikes_position4=max(new_spikes_position4);
        end
        new_spikes_position_14=[new_spikes_position_14 spikes_position_1(i)];
        new_spikes_position_41=[new_spikes_position_41 new_spikes_position4];
    end
    clear new_spikes_position2;
    clear new_spikes_position3;
    clear new_spikes_position4;
end

for i=1:length(spikes_position_2)
    new_spikes_position3=spikes_position_3(spikes_position_3>=spikes_position_2(i)-window_size& spikes_position_3<=spikes_position_2(i)+window_size);
    new_spikes_position4=spikes_position_4(spikes_position_4>=spikes_position_2(i)-window_size& spikes_position_4<=spikes_position_2(i)+window_size);
    
   if isempty(new_spikes_position3)
   else
       if new_spikes_position3>=spikes_position_2(i)
            new_spikes_position3=min(new_spikes_position3);
         else
             new_spikes_position3=max(new_spikes_position3);
       end
        new_spikes_position_23=[new_spikes_position_23 spikes_position_2(i)];
        new_spikes_position_32=[new_spikes_position_32 new_spikes_position3];
   end
    
   if isempty(new_spikes_position4)
   else
        if new_spikes_position4>=spikes_position_2(i)
            new_spikes_position4=min(new_spikes_position4);
        else
            new_spikes_position4=max(new_spikes_position4);
        end
        new_spikes_position_24=[new_spikes_position_24 spikes_position_2(i)];
        new_spikes_position_42=[new_spikes_position_42 new_spikes_position4];
    end
    clear new_spikes_position3;
    clear new_spikes_position4;
end


for i=1:length(spikes_position_3)
    new_spikes_position4=spikes_position_4(spikes_position_4>=spikes_position_3(i)-window_size& spikes_position_4<=spikes_position_3(i)+window_size);       
   if isempty(new_spikes_position4)
   else
        if new_spikes_position4>=spikes_position_3(i)
            new_spikes_position4=min(new_spikes_position4);
        else
            new_spikes_position4=max(new_spikes_position4);
        end
        new_spikes_position_34=[new_spikes_position_34 spikes_position_3(i)];
        new_spikes_position_43=[new_spikes_position_43 new_spikes_position4];
    end
    clear new_spikes_position4;
end

common_phase_12=phase_data(new_spikes_position_12,2);
common_phase_21=phase_data(new_spikes_position_21,2);
common_phase_13=phase_data(new_spikes_position_13,2);
common_phase_31=phase_data(new_spikes_position_31,2);
common_phase_14=phase_data(new_spikes_position_14,2);
common_phase_41=phase_data(new_spikes_position_41,2);
common_phase_23=phase_data(new_spikes_position_23,2);
common_phase_32=phase_data(new_spikes_position_32,2);
common_phase_24=phase_data(new_spikes_position_24,2);
common_phase_42=phase_data(new_spikes_position_42,2);
common_phase_34=phase_data(new_spikes_position_34,2);
common_phase_43=phase_data(new_spikes_position_43,2);



binWidth=0.005;

subplot(2,3,1)
sort_diff_1=common_phase_12-common_phase_21;
binCtrs=-0.3:binWidth: 0.3;
counts1=hist(sort_diff_1,binCtrs);
bar(binCtrs, counts1/sum(counts1));
xlim([-0.15  0.15]);
ylim([0 0.2]);
set(gca,'FontSize',20)
xlabel(' theta phase difference','FontSize',20)

subplot(2,3,2)
sort_diff_2=common_phase_13-common_phase_31;
binCtrs=-0.3:binWidth: 0.3;
sort_diff_2=sort_diff_2(abs(sort_diff_2)<=1);
counts2=hist(sort_diff_2,binCtrs);
bar(binCtrs, counts2/sum(counts2));
xlim([-0.15  0.15]);
ylim([0 0.2]);
set(gca,'FontSize',20)
xlabel(' theta phase difference','FontSize',20)

subplot(2,3,3)
sort_diff_3=common_phase_14-common_phase_41;
sort_diff_3=sort_diff_3(abs(sort_diff_3)<=1);
binCtrs=-0.3:binWidth: 0.3;
counts3=hist(sort_diff_3,binCtrs);
bar(binCtrs, counts3/sum(counts3));
xlim([-0.15  0.15]);
ylim([0 0.2]);
set(gca,'FontSize',20)
xlabel(' theta phase difference','FontSize',20)

subplot(2,3,4)
sort_diff_4=common_phase_23-common_phase_32;
sort_diff_4=sort_diff_4(abs(sort_diff_4)<=1);
binCtrs=-0.3:binWidth: 0.3;
counts4=hist(sort_diff_4,binCtrs);
bar(binCtrs, counts4/sum(counts4));
xlim([-0.15  0.15]);
ylim([0 0.2]);
set(gca,'FontSize',20)
xlabel(' theta phase difference','FontSize',20)


subplot(2,3,5)
sort_diff_5=common_phase_24-common_phase_42;
sort_diff_5=sort_diff_5(abs(sort_diff_5)<=1);
binCtrs=-0.3:binWidth: 0.3;
counts5=hist(sort_diff_5,binCtrs);
bar(binCtrs, counts4/sum(counts4));
xlim([-0.15  0.15]);
ylim([0 0.2]);
set(gca,'FontSize',20)
xlabel(' theta phase difference','FontSize',20)


subplot(2,3,6)
sort_diff_6=common_phase_34-common_phase_43;
sort_diff_6=sort_diff_6(abs(sort_diff_6)<=1);
binCtrs=-0.3:binWidth: 0.3;
counts6=hist(sort_diff_6,binCtrs);
bar(binCtrs, counts4/sum(counts4));
xlim([-0.15  0.15]);
ylim([0 0.2]);
set(gca,'FontSize',20)
xlabel(' theta phase difference','FontSize',20)



%% june 24 baseline_1
plot((1:length(data))/25000,data(:,1))
set(gca,'FontSize',20)
xlabel('Time(second)','FontSize',20)
ylabel('Voltage(mV)','FontSize',20)




%% dec_02_10


%theta rythme

%% feb_22_11 seizure 2
% theta rythme
data=xin_seizure_2_MultiCH_Data(:,:);
data=data(10.2*25000:length(data),:);

data(:,5)=data(:,4);

Hd=bpf1;
data(:,4)=filter(Hd, data(:,4));


theta_data=data(:,5);
% theta_data=xin_seizure_10_MultiCH_Data(10.2*25000:50*25000,4);
Hd=mytheta;
theta_data=filter(Hd,theta_data); %butterworth

plot( (1:length(data))/25000,data(:,5)+2)
hold on
plot( (1:length(data))/25000,theta_data,'k')

h = legend('EEG','Theta Oscillation',2);
set(h,'Interpreter','none','FontSize',20)
set(gca,'FontSize',20)
xlabel('Time(second)','FontSize',20)
ylabel('Voltage (mV)','FontSize',20)


 plot( (1:length(data))/25000,data(:,1)+0.5)
hold on
plot( (1:length(data))/25000,data(:,2))
plot( (1:length(data))/25000,data(:,3)-0.5)
plot( (1:length(data))/25000,data(:,4)-1)
plot( (1:length(data))/25000,data(:,5)-1.5,'k')
ylim([-2 1])
set(gca,'FontSize',20)
xlabel('Time(second)','FontSize',20)
ylabel('Amplitude (mV)','FontSize',20)

% find spikes
peak_threshold=-0.05;%threshold to pick the spike
period=[];
min_num_pts=5;
spikes_position_1=find_spikes(-abs(data(:,1)), peak_threshold,min_num_pts);
spikes_position_2=find_spikes(-abs(data(:,2)), peak_threshold,min_num_pts);
spikes_position_3=find_spikes(-abs(data(:,3)), peak_threshold,min_num_pts);
spikes_position_4=find_spikes(-abs(data(:,4)), peak_threshold,min_num_pts);

 


% find spike for theta
plot( (1:length(data))/25000,theta_data)
hold on
peak_threshold=-0.05;%threshold to pick the spike
period=[];
min_num_pts=1;
theta_spikes_position=find_spikes((theta_data(:,1)), peak_threshold,min_num_pts);
plot(theta_spikes_position/25000,theta_data(theta_spikes_position,1),'*','MarkerEdgeColor','r')


%define interval
phase_data=[];
phase_data(:,1)=(1:length(data))';% 1st time
%2nd give each time point a phase value
for k=1:theta_spikes_position(1)
phase_data(k,2)=2*pi/theta_spikes_position(1)*k;
end
for i=1:(length(theta_spikes_position)-1)
    for j=theta_spikes_position(i):theta_spikes_position(i+1)
    phase_data(j,2)=2*pi/(theta_spikes_position(i+1)-theta_spikes_position(i))*(j-theta_spikes_position(i));
    end
end
for z=theta_spikes_position(length(theta_spikes_position)):length(data)
    phase_data(z,2)=2*pi/(length(data)-theta_spikes_position(length(theta_spikes_position)))*(z-theta_spikes_position(length(theta_spikes_position)));
end

plot( (1:length(data))/25000,theta_data)
hold on
plot(phase_data(:,1)/25000,phase_data(:,2),'k')
hold on
plot(theta_spikes_position/25000,theta_data(theta_spikes_position,1),'*','MarkerEdgeColor','r')
h = legend('Theta oscillation','Theta Phase',2);
set(h,'Interpreter','none','FontSize',20)
set(gca,'FontSize',20)
xlabel('Time(second)','FontSize',20)


% give spike-time a phase value
phase_ch1=phase_data(spikes_position_1,2);
phase_ch2=phase_data(spikes_position_2,2);
phase_ch3=phase_data(spikes_position_3,2);
phase_ch4=phase_data(spikes_position_4,2);

plot(phase_data(:,1)/25000,phase_data(:,2),'k')
hold on;plot(spikes_position_1/25000,phase_ch1,'*','MarkerEdgeColor','k')
hold on; plot(spikes_position_2/25000,phase_ch2,'*','MarkerEdgeColor','r')
hold on; plot(spikes_position_3/25000,phase_ch3,'*','MarkerEdgeColor','g')
hold on; plot(spikes_position_4/25000,phase_ch4,'*','MarkerEdgeColor','m')

h = legend('Theta Phase','ch1_spike_phase','ch2_spike_phase','ch3_spike_phase','ch4_spike_phase',2);
set(h,'Interpreter','none','FontSize',20)
set(gca,'FontSize',20)
xlabel('Time(second)','FontSize',20)
ylabel('phase','FontSize',20)

%%%%%%%%%%%%%%%
binWidth=0.2;

subplot(2,2,1)
sort_phase_ch1=sort(phase_ch1);
binCtrs=0:binWidth:2*pi;
counts1=hist(sort_phase_ch1,binCtrs);
bar(binCtrs, counts1/sum(counts1));
xlim([0 2*pi]);
ylim([0 0.15]);
set(gca,'FontSize',20)
set(gca,'XTick',[0 pi 2*pi])


subplot(2,2,2)
sort_phase_ch2=sort(phase_ch2);
binCtrs=0:binWidth:2*pi;
counts2=hist(sort_phase_ch2,binCtrs);
bar(binCtrs, counts2/sum(counts2));
xlim([0 2*pi]);
ylim([0 0.15]);
set(gca,'FontSize',20)
set(gca,'XTick',[0 pi 2*pi])

subplot(2,2,3)
sort_phase_ch3=sort(phase_ch3);
binCtrs=0:binWidth:2*pi;
counts3=hist(sort_phase_ch3,binCtrs);
bar(binCtrs, counts3/sum(counts3));
xlim([0 2*pi]);
ylim([0 0.15]);
set(gca,'FontSize',20)
set(gca,'XTick',[0 pi 2*pi])

subplot(2,2,4)
sort_phase_ch4=sort(phase_ch4);
binCtrs=0:binWidth:2*pi;
counts4=hist(sort_phase_ch4,binCtrs);
bar(binCtrs, counts4/sum(counts4));
xlim([0 2*pi]);
ylim([0 0.15]);
set(gca,'FontSize',20)
set(gca,'XTick',[0 pi 2*pi])    

%%%%%%%%%
% phase difference

window_size=25;
new_spikes_position_1=[];
new_spikes_position_2=[];
new_spikes_position_3=[];
new_spikes_position_4=[];
for i=1:length(spikes_position_1)
    new_spikes_position2=spikes_position_2(spikes_position_2>=spikes_position_1(i)-window_size & spikes_position_2<=spikes_position_1(i)+window_size);
    new_spikes_position3=spikes_position_3(spikes_position_3>=spikes_position_1(i)-window_size& spikes_position_3<=spikes_position_1(i)+window_size);
    new_spikes_position4=spikes_position_4(spikes_position_4>=spikes_position_1(i)-window_size& spikes_position_4<=spikes_position_1(i)+window_size);
    if isempty(new_spikes_position2)|isempty(new_spikes_position3)|isempty(new_spikes_position4)
    else  
        if new_spikes_position2>=spikes_position_1(i)
            new_spikes_position2=min(new_spikes_position2);
        else
            new_spikes_position2=max(new_spikes_position2);
        end
        if new_spikes_position3>=spikes_position_1(i)
            new_spikes_position3=min(new_spikes_position3);
        else
            new_spikes_position3=max(new_spikes_position3);
        end
        if new_spikes_position4>=spikes_position_1(i)
            new_spikes_position4=min(new_spikes_position4);
        else
            new_spikes_position4=max(new_spikes_position4);
        end
    new_spikes_position_1=[new_spikes_position_1 spikes_position_1(i)];
    new_spikes_position_2=[new_spikes_position_2 new_spikes_position2];
    new_spikes_position_3=[new_spikes_position_3 new_spikes_position3];
    new_spikes_position_4=[new_spikes_position_4 new_spikes_position4];
    end
    clear new_spikes_position2;
    clear new_spikes_position3;
    clear new_spikes_position3;
end

plot((1:length(data))/25000,data(:,1)+1)
hold on
plot((1:length(data))/25000,data(:,2)+0.5)
plot((1:length(data))/25000,data(:,3))
plot((1:length(data))/25000,data(:,4)-0.5)

plot(new_spikes_position_1/25000,data(new_spikes_position_1,1)+1,'*','MarkerEdgeColor','k')
hold on
plot(new_spikes_position_2/25000,data(new_spikes_position_2,2)+0.5,'*','MarkerEdgeColor','r')
plot(new_spikes_position_3/25000,data(new_spikes_position_3,3),'*','MarkerEdgeColor','g')
plot(new_spikes_position_4/25000,data(new_spikes_position_4,4)-0.5,'*','MarkerEdgeColor','m')

%%%%% find theta phase difference

window_size=50;
new_spikes_position_12=[];
new_spikes_position_21=[];
new_spikes_position_13=[];
new_spikes_position_31=[];
new_spikes_position_14=[];
new_spikes_position_41=[];
new_spikes_position_23=[];
new_spikes_position_32=[];
new_spikes_position_24=[];
new_spikes_position_42=[];
new_spikes_position_34=[];
new_spikes_position_43=[];

for i=1:length(spikes_position_1)
    new_spikes_position2=spikes_position_2(spikes_position_2>=spikes_position_1(i)-window_size & spikes_position_2<=spikes_position_1(i)+window_size);
    new_spikes_position3=spikes_position_3(spikes_position_3>=spikes_position_1(i)-window_size& spikes_position_3<=spikes_position_1(i)+window_size);
    new_spikes_position4=spikes_position_4(spikes_position_4>=spikes_position_1(i)-window_size& spikes_position_4<=spikes_position_1(i)+window_size);
%     if isempty(new_spikes_position2)|isempty(new_spikes_position3)|isempty(new_spikes_position4)
    if isempty(new_spikes_position2)
    else  
        if new_spikes_position2>=spikes_position_1(i)
            new_spikes_position2=min(new_spikes_position2);
        else
            new_spikes_position2=max(new_spikes_position2);
        end
    new_spikes_position_12=[new_spikes_position_12 spikes_position_1(i)];
    new_spikes_position_21=[new_spikes_position_21 new_spikes_position2];
    end
    
   if isempty(new_spikes_position3)
   else
       if new_spikes_position3>=spikes_position_1(i)
            new_spikes_position3=min(new_spikes_position3);
         else
             new_spikes_position3=max(new_spikes_position3);
       end
        new_spikes_position_13=[new_spikes_position_13 spikes_position_1(i)];
        new_spikes_position_31=[new_spikes_position_31 new_spikes_position3];
   end
    
   if isempty(new_spikes_position4)
   else
        if new_spikes_position4>=spikes_position_1(i)
            new_spikes_position4=min(new_spikes_position4);
        else
            new_spikes_position4=max(new_spikes_position4);
        end
        new_spikes_position_14=[new_spikes_position_14 spikes_position_1(i)];
        new_spikes_position_41=[new_spikes_position_41 new_spikes_position4];
    end
    clear new_spikes_position2;
    clear new_spikes_position3;
    clear new_spikes_position4;
end

for i=1:length(spikes_position_2)
    new_spikes_position3=spikes_position_3(spikes_position_3>=spikes_position_2(i)-window_size& spikes_position_3<=spikes_position_2(i)+window_size);
    new_spikes_position4=spikes_position_4(spikes_position_4>=spikes_position_2(i)-window_size& spikes_position_4<=spikes_position_2(i)+window_size);
    
   if isempty(new_spikes_position3)
   else
       if new_spikes_position3>=spikes_position_2(i)
            new_spikes_position3=min(new_spikes_position3);
         else
             new_spikes_position3=max(new_spikes_position3);
       end
        new_spikes_position_23=[new_spikes_position_23 spikes_position_2(i)];
        new_spikes_position_32=[new_spikes_position_32 new_spikes_position3];
   end
    
   if isempty(new_spikes_position4)
   else
        if new_spikes_position4>=spikes_position_2(i)
            new_spikes_position4=min(new_spikes_position4);
        else
            new_spikes_position4=max(new_spikes_position4);
        end
        new_spikes_position_24=[new_spikes_position_24 spikes_position_2(i)];
        new_spikes_position_42=[new_spikes_position_42 new_spikes_position4];
    end
    clear new_spikes_position3;
    clear new_spikes_position4;
end


for i=1:length(spikes_position_3)
    new_spikes_position4=spikes_position_4(spikes_position_4>=spikes_position_3(i)-window_size& spikes_position_4<=spikes_position_3(i)+window_size);       
   if isempty(new_spikes_position4)
   else
        if new_spikes_position4>=spikes_position_3(i)
            new_spikes_position4=min(new_spikes_position4);
        else
            new_spikes_position4=max(new_spikes_position4);
        end
        new_spikes_position_34=[new_spikes_position_34 spikes_position_3(i)];
        new_spikes_position_43=[new_spikes_position_43 new_spikes_position4];
    end
    clear new_spikes_position4;
end

common_phase_12=phase_data(new_spikes_position_12,2);
common_phase_21=phase_data(new_spikes_position_21,2);
common_phase_13=phase_data(new_spikes_position_13,2);
common_phase_31=phase_data(new_spikes_position_31,2);
common_phase_14=phase_data(new_spikes_position_14,2);
common_phase_41=phase_data(new_spikes_position_41,2);
common_phase_23=phase_data(new_spikes_position_23,2);
common_phase_32=phase_data(new_spikes_position_32,2);
common_phase_24=phase_data(new_spikes_position_24,2);
common_phase_42=phase_data(new_spikes_position_42,2);
common_phase_34=phase_data(new_spikes_position_34,2);
common_phase_43=phase_data(new_spikes_position_43,2);



binWidth=0.005;

subplot(2,3,1)
sort_diff_1=common_phase_12-common_phase_21;
binCtrs=-0.3:binWidth: 0.3;
counts1=hist(sort_diff_1,binCtrs);
bar(binCtrs, counts1/sum(counts1));
xlim([-0.15  0.15]);
ylim([0 0.2]);
set(gca,'FontSize',20)
xlabel(' theta phase difference','FontSize',20)

subplot(2,3,2)
sort_diff_2=common_phase_13-common_phase_31;
binCtrs=-0.3:binWidth: 0.3;
sort_diff_2=sort_diff_2(abs(sort_diff_2)<=1);
counts2=hist(sort_diff_2,binCtrs);
bar(binCtrs, counts2/sum(counts2));
xlim([-0.15  0.15]);
ylim([0 0.2]);
set(gca,'FontSize',20)
xlabel(' theta phase difference','FontSize',20)

subplot(2,3,3)
sort_diff_3=common_phase_14-common_phase_41;
sort_diff_3=sort_diff_3(abs(sort_diff_3)<=1);
binCtrs=-0.3:binWidth: 0.3;
counts3=hist(sort_diff_3,binCtrs);
bar(binCtrs, counts3/sum(counts3));
xlim([-0.15  0.15]);
ylim([0 0.2]);
set(gca,'FontSize',20)
xlabel(' theta phase difference','FontSize',20)

subplot(2,3,4)
sort_diff_4=common_phase_23-common_phase_32;
sort_diff_4=sort_diff_4(abs(sort_diff_4)<=1);
binCtrs=-0.3:binWidth: 0.3;
counts4=hist(sort_diff_4,binCtrs);
bar(binCtrs, counts4/sum(counts4));
xlim([-0.15  0.15]);
ylim([0 0.2]);
set(gca,'FontSize',20)
xlabel(' theta phase difference','FontSize',20)


subplot(2,3,5)
sort_diff_5=common_phase_24-common_phase_42;
sort_diff_5=sort_diff_5(abs(sort_diff_5)<=1);
binCtrs=-0.3:binWidth: 0.3;
counts5=hist(sort_diff_5,binCtrs);
bar(binCtrs, counts4/sum(counts4));
xlim([-0.15  0.15]);
ylim([0 0.2]);
set(gca,'FontSize',20)
xlabel(' theta phase difference','FontSize',20)


subplot(2,3,6)
sort_diff_6=common_phase_34-common_phase_43;
sort_diff_6=sort_diff_6(abs(sort_diff_6)<=1);
binCtrs=-0.3:binWidth: 0.3;
counts6=hist(sort_diff_6,binCtrs);
bar(binCtrs, counts4/sum(counts4));
xlim([-0.15  0.15]);
ylim([0 0.2]);
set(gca,'FontSize',20)
xlabel(' theta phase difference','FontSize',20)




%% mar_23_11 
% baseline 0

data(:,5)=data(:,4);
Hd=bpf1;
data(:,4)=filter(Hd, data(:,4));
baseline_data=data;

plot((1:2500000)/25000,data(1:2500000,1)+0.5)
hold on
plot((1:2500000)/25000,data(1:2500000,2))
plot((1:2500000)/25000,data(1:2500000,3)-0.5)
plot((1:2500000)/25000,data(1:2500000,4)-1)

plot((2500000:5000000)/25000,data(2500000:5000000,1)+0.5)
hold on
plot((2500000:5000000)/25000,data(2500000:5000000,2))
plot((2500000:5000000)/25000,data(2500000:5000000,3)-0.5)
plot((2500000:5000000)/25000,data(2500000:5000000,4)-1)

plot((5000000:7500000)/25000,data(5000000:7500000,1)+0.5)
hold on
plot((5000000:7500000)/25000,data(5000000:7500000,2))
plot((5000000:7500000)/25000,data(5000000:7500000,3)-0.5)
plot((5000000:7500000)/25000,data(5000000:7500000,4)-1)


sample_data_1=baseline_data(1:250000,1);
sample_data_2=baseline_data(1:250000,2);
sample_data_3=baseline_data(1:250000,3);
sample_data_4=baseline_data(1:250000,4);
hilbert_sample_1=hilbert(sample_data_1);
hilbert_sample_2=hilbert(sample_data_2);
hilbert_sample_3=hilbert(sample_data_3);
hilbert_sample_4=hilbert(sample_data_4);

r1=sqrt(real(hilbert_sample_1).^2+imag(hilbert_sample_1).^2);
r2=sqrt(real(hilbert_sample_2).^2+imag(hilbert_sample_2).^2);
r3=sqrt(real(hilbert_sample_3).^2+imag(hilbert_sample_3).^2);
r4=sqrt(real(hilbert_sample_4).^2+imag(hilbert_sample_4).^2);
mean_r=[mean(r1),mean(r2),mean(r3),mean(r4)];
r=min(mean_r);
rms=sqrt(mean(r.^2)); %rms


% seizure 4

data=xin_seizure_4_MultiCH_Data(:,:);
data=data(10.2*25000:length(data),:);

data(:,5)=data(:,1);
data(:,6)=data(:,2);
data(:,7)=data(:,3);
data(:,8)=data(:,4);

Hd=bpf1;
data(:,1)=filter(Hd, data(:,1));
data(:,2)=filter(Hd, data(:,2));
data(:,3)=filter(Hd, data(:,3));
data(:,4)=filter(Hd, data(:,4));

hilbert_1=hilbert(data(:,1));
hilbert_2=hilbert(data(:,2));
hilbert_3=hilbert(data(:,3));
hilbert_4=hilbert(data(:,4));

% plot(hilbert_4,'*','MarkerEdgeColor','r')
% hold on
% plot(hilbert_sample_1,'*')

r1=sqrt(real(hilbert_1).^2 + imag(hilbert_1).^2);
r2=sqrt(real(hilbert_2).^2 + imag(hilbert_2).^2);
r3=sqrt(real(hilbert_3).^2 + imag(hilbert_3).^2);
r4=sqrt(real(hilbert_4).^2 + imag(hilbert_4).^2);
hilbert_1( r1< rms ) = NaN;
hilbert_2( r2< rms)=NaN;
hilbert_3( r3< rms ) = NaN;
hilbert_4( r4< rms)=NaN;

ksdensity(diff(unwrap(mod((angle(hilbert_sample_1)+pi),2*pi)-pi)))
hold on

ksdensity(diff(unwrap(mod((angle(hilbert_1)+pi),2*pi)-pi)))
hold all
ksdensity(diff(unwrap(mod((angle(hilbert_2)+pi),2*pi)-pi)))
ksdensity(diff(unwrap(mod((angle(hilbert_3)+pi),2*pi)-pi)))
ksdensity(diff(unwrap(mod((angle(hilbert_4)+pi),2*pi)-pi)))
h = legend('baseline','ch1','ch2','ch3','ch4',2);
set(h,'Interpreter','none','FontSize',20)
set(gca,'FontSize',20)
xlabel('phase velocity')
ylabel('Density')