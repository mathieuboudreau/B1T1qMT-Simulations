function [sled_flips sled_offsets sled_pulse_durations sled_trs sled_excite_flips] = create_protocol_sled(inputArray)

%% Create Protocol Struct
% 
% LAST EDITED BY: Mathieu Boudreau (_MJB)
% LAST MODIFIED ON: August 10nth 2012

sled_offsets = unique(inputArray(:,2))';
sled_pulse_durations = unique(inputArray(:,3));
sled_trs = unique(inputArray(:,4));
sled_excite_flips = unique(inputArray(:,5));

sled_flips = cell(length(sled_trs),1);

for trIndex = 1:length(sled_trs)
    count = 1;
    for forIndex = 1:length(inputArray)
        if inputArray(forIndex,4) == sled_trs(trIndex)
            temp_sled_flips(count) = inputArray(forIndex,1);
            count = count+1;
        end
    end
   temp_sled_flips = unique(temp_sled_flips);
   sled_flips{trIndex} = {temp_sled_flips};
end
