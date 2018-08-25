
function [mean_1, mean_2, mean_3,mean_4] = mean_ShanEnt(phase_ch1,phase_ch2,phase_ch3,phase_ch4,window_size,slide,number_bins)
% Calculate shannon entropy for each channel
% phase_ch1,2,3,4--theta phase value for all the peaks
% Shannon Entropy=-sum(P.*log(P)); P--probability of each bin

ShanEnt=[];
for i=1:round((length(phase_ch1)-window_size)/slide) 
    sample_window=phase_ch1(1+slide*(i-1):window_size+slide*(i-1));
    P = hist(sample_window,number_bins)/window_size;
        index_zeros = -nd(P==0);
        P(index_zeros) = [];
    ShanEnt(i)=-sum(P.*log(P));
end
mean_1=mean(ShanEnt);
ShanEnt=[];

for i=1:round((length(phase_ch2)-window_size)/slide)
    sample_window=phase_ch2(1+slide*(i-1):window_size+slide*(i-1));
    P= hist(sample_window,number_bins)/window_size;
        index_zeros = -nd(P==0);
        P(index_zeros) = [];
    ShanEnt(i)=-sum(P.*log(P));
end
mean_2=mean(ShanEnt);
ShanEnt=[];


for i=1:round((length(phase_ch3)-window_size)/slide)
    sample_window=phase_ch3(1+slide*(i-1):window_size+slide*(i-1)); 
    P= hist(sample_window,number_bins)/window_size;
        index_zeros = -nd(P==0);
        P(index_zeros) = [];
    ShanEnt(i)=-sum(P.*log(P));
end
mean_3=mean(ShanEnt);
ShanEnt=[];


for i=1:round((length(phase_ch4)-window_size)/slide) 
    sample_window=phase_ch4(1+slide*(i-1):window_size+slide*(i-1));
    P= hist(sample_window,number_bins)/window_size;
        index_zeros = -nd(P==0);
        P(index_zeros) = [];
    ShanEnt(i)=-sum(P.*log(P));
end
mean_4=mean(ShanEnt);