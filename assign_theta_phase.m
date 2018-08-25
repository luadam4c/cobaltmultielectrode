function [phase_data,phase_ch1,phase_ch2,phase_ch3,phase_ch4]=theta_dist(binWidth,new,theta_spikes_position,spikes_position_1,spikes_position_2,spikes_position_3,spikes_position_4)
% Assign theta phase value to all the peaks.
% new--filtered data(band pass: 300-4000Hz)
% theta_spikes_position--define theta cycle
% spikes_position_1,2,3,4--all the peaks
binWidth=0.2;
phase_data=[];
phase_data(:,1)=(1:length(new))'; 
%give each time point a phase value
for k=1:theta_spikes_position(1) 
    phase_data(k,2)=2*pi/theta_spikes_position(1)*k;
end
for i=1:(length(theta_spikes_position)-1)
    for j=theta_spikes_position(i):theta_spikes_position(i+1) phase_data(j,2)=2*pi/(theta_spikes_position(i+1)-theta_spikes_position(i))*(j-theta_spikes_position(i));
    end
end

for z=theta_spikes_position(length(theta_spikes_position)):length(new)
    phase_data(z,2)=2*pi/(length(new)-theta_spikes_position(length(theta_spikes_position)))*(z-theta_spikes_position(length(theta_spikes_position)));
end
phase_ch1=phase_data(spikes_position_1,2);
phase_ch2=phase_data(spikes_position_2,2);
phase_ch3=phase_data(spikes_position_3,2);
phase_ch4=phase_data(spikes_position_4,2);

%plot theta phase distribution
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
