% Event Synchronization
function [sum_J] = event_synch(action_potentials_1,action_potentials_2,lag_index)
% action_potentials_1 , action_potentials_2 --peaks in two time series

sum_J = 0; %Initialize variable

for counter = 1:length(action_potentials_1) %Look at every action potential
    %Verify that lag_index is smaller than half the interval between 2 action potentials
    if counter <= length(action_potentials_1)-1
        min_lag = min([lag_index,(action_potentials_1(1,counter+1)-action_potentials_1(1,counter))/2]);
    else
        min_lag = lag_index;
    end
    
    %Find a action potential of 2 that is between the action potential in 1 and the
        lag_index
    ind = find(action_potentials_2 (1,:) <= action_potentials_1(1,counter)+min_lag & action_potentials_2(1,:) > action_potentials_1(1,counter), 1);
    if isempty(ind)
        %See if the action potentials in 1 and 2 have the same index
        ind = find(action_potentials_2 (1,:) == action_potentials_1(1,counter), 1);
        if isempty(ind)
            J = 0; %No action potential in 2 after action potential in 1
            sum J = sum J+J;
        else
            J = .5; %Action potentials occur at the same time
            sum J = sum J+J;
        end
    else
        %Recalulate lag_index
if ind>=2 && counter>=2 && counter~=length(action_potentials_1) && ind~= length(action_potentials_2)
    min_lag = min([lag_index,(action_potentials_1(1,counter+1)- ...
    action_potentials_1(1,counter))/2,(action_potentials_2(1,ind+1)- ...
    action_potentials_2(1,ind))/2,(action_potentials_1(1,counter)- ...
    action_potentials_1(1,counter-1))/2,(action_potentials_2(1,ind)- ...
    action_potentials_2(1,ind-1))/2]); 
elseif ind>=2 && counter==length(action_potentials_1) && ind~=length(action_potentials_2)
    min_lag = min([lag_index,(action_potentials_2(1,ind+1)- ...
    action_potentials_2 (1,ind))/2,(action_potentials_1(1,counter)- ...
    action_potentials_1(1,counter -1))/2,(action_potentials_2(1,ind)- ...
    action_potentials_2(1,ind-1))/2]); 
elseif counter>=2 && counter~=length(action_potentials_1) && ind==length(action_potentials_2)
    min_lag = min([lag_index,(action_potentials_1(1,counter+1)- ...
    action_potentials_1(1,counter))/2,(action_potentials_1(1,counter)- ...
    action_potentials_1(1,counter-1))/2,(action_potentials_2(1,ind)- ...
    action_potentials_2(1,ind-1))/2]);
elseif ind<2 && counter<2
    min_lag = min([lag_index,(action_potentials_1(1,counter+1)- ...
    action_potentials_1(1,counter))/2,(action_potentials_2(1,ind+1)- ...
    action_potentials_2(1,ind))/2]);
elseif ind<2 && counter~=length(action_potentials_1)
    min_lag = min([lag_index,(action_potentials_1(1,counter+1)- ...
    action_potentials_1(1,counter))/2,(action_potentials_2(1,ind+1)- ...
    action_potentials_2(1,ind))/2,(action_potentials_1(1,counter)- ...
    action_potentials_1(1,counter-1))/2]);
elseif counter<2 && ind~=length(action_potentials_2)
    min_lag = min([lag_index,(action_potentials_1(1,counter+1)- ...
    action_potentials_1(1,counter))/2,(action_potentials_2(1,ind+1)- ...
    action_potentials_2(1,ind))/2,(action_potentials_2(1,ind)- ...
    action_potentials_2(1,ind-1))/2]);
elseif counter==length(action_potentials_1) && ind==length(action_potentials_2)
    min_lag = min([lag_index,(action_potentials_1(1,counter)- ...
    action_potentials_1(1,counter-1))/2,(action_potentials_2(1,ind)- ...
    action_potentials_2(1,ind-1))/2]); 
elseif counter<2 && ind==length(action_potentials_2)
    min_lag = min([lag_index,(action_potentials_1(1,counter+1)- ...
    action_potentials_1(1,counter))/2,(action_potentials_2(1,ind)- ...
    action_potentials_2(1,ind-1))/2]);
else %counter at end and ind<2
    min_lag = min([lag_index,(action_potentials_2(1,ind+1)- ...
    action_potentials_2 (1,ind))/2,(action_potentials_1(1,counter)- ...
    action_potentials_1(1,counter -1))/2]);
end

    %Find a action potential of 2 that is between the action potential in 1 and the
        lag_index
        ind = find(action_potentials_2 (1,:) <= action_potentials_1(1,counter)+min_lag & action_potentials_2(1,:) > action_potentials_1(1,counter), 1);
        if isempty(ind) J = 0; %No action potential in 2 after action potential in 1
            sum J = sum J+J;
        else
            J = 1; %Action potential in 2 is within the lag_index of 1
            sum J = sum J+J;
        end
    
    end
end