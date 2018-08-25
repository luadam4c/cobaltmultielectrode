window_size=30*25; % window length
slide=125; % sliding 

bin_number = 30; %set bin number
number_bins = linspace(-pi,pi,bin_number); % set range
ShanEnt=[];

% i=1000;
% 
% sample_window=delta_phi(1+slide*(i-1):window_size+slide*(i-1));
%     P= hist(sample_window,number_bins)/window_size;
%         index_zeros = find(P==0);
%         P(index_zeros) = [];
%     ShanEnt(i)=-bin_number*sum(P.*log(P));
%     
%     hist(sample_window,number_bins)
%     
    
for i=1:round((length(delta_phi)-window_size)/slide)
    sample_window=delta_phi(1+slide*(i-1):window_size+slide*(i-1));
    P= hist(sample_window,number_bins)/window_size;
        index_zeros = find(P==0);
        P(index_zeros) = [];
    ShanEnt(i)=-sum(P.*log(P));
end
% 


time=(1:slide:round((length(sample_data_2)-window_size)/slide)*slide)+window_size/2;
ax(1)=subplot(3,1,1);
plot((1:length(sample_data_2))/25000,sample_data_2)
hold on
plot((1:length(sample_data_4))/25000,sample_data_4-1)
hold off
ylabel('Voltage(V)','FontSize',20)

ax(2)=subplot(3,1,2);
plot((1:length(sample_data_2))/25000,delta_phi)
ylabel('delta\_phi','FontSize',20)

ax(3)=subplot(3,1,3);
plot(time/25000, ShanEnt)
linkaxes(ax,'x')
xlabel('time(second)','FontSize',20)
ylabel('ShanEnt','FontSize',20)


