function [ output_args ] = find_spikes(fft_data, peak_threshold, min_num_pts)

indeces = [1:1:length(fft_data)];
    y = fft_data(:,1);
    max_data = [ ] ;
    spike_indeces=[ ];
    peak_data_indeces = y<peak_threshold;
    peak_y = y(peak_data_indeces);
    peak_indeces = indeces(peak_data_indeces);
    group_break_indeces = find(diff(peak_indeces)>1);
    
    for group_counter = 2:length(group_break_indeces)
        
        block_indeces = [group_break_indeces(group_counter-1)+1:group_break_indeces(group_counter)];
        data_block_y = peak_y(block_indeces);
        data_block_x = peak_indeces(block_indeces)';
        
        if length(data_block_y)>min_num_pts
            max_indeces=find(data_block_y==min(data_block_y));
            spike_indeces_temp=data_block_x(max_indeces);
            spike_indeces=[spike_indeces; spike_indeces_temp(1)];
        end 
    end
    
    output_args = spike_indeces;
    
