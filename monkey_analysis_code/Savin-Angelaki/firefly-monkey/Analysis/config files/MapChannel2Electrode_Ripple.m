function [channel_id ,electrode_id, brain_area, channels_per_area, electrode_type] = MapChannel2Electrode_Ripple(MetaTags, prs)
channels_per_area = [];
for a = 1:size(prs.area, 2)
    
    probe_in_area = prs.electrode_type{a};
    
    if strcmp(probe_in_area, 'linearprobe12')
        channels_per_area(a) = 12;
    elseif strcmp(probe_in_area, 'linearprobe16')
        channels_per_area(a) = 16;
    elseif strcmp(probe_in_area, 'linearprobe24')
        channels_per_area(a) = 24;
    elseif strcmp(probe_in_area, 'linearprobe32')
        channels_per_area(a) = 32;
    elseif strcmp(probe_in_area, 'linearprobe64')
        channels_per_area(a) = 64;
    elseif strcmp(probe_in_area, 'utah48')
        channels_per_area(a) = 48;
    elseif strcmp(probe_in_area, 'utah96')
        channels_per_area(a) = 96;
    elseif strcmp(probe_in_area, 'utah128')
        channels_per_area(a) = 128;
    else
        disp('WARNING; dont recognize the probe type')
        break
    end
end

transition_points = [0 cumsum(channels_per_area)];
transition_points(end) = [];
for a = 1:sum(channels_per_area)
    channel_id(a) = MetaTags.ChannelID(a);
    area_index = double(find(a>transition_points, 1, 'last'));
    brain_area(a) = prs.area(area_index);
    electrode_type(a) = prs.electrode_type(area_index);
end

electrode_id = []; 
for a = 1:size(prs.area, 2)
    mapping = load([prs.filepath_conf, [prs.electrode_type{a}, '.mat']]);
    electrode_id = [electrode_id, mapping.map.channel(mapping.map.electrode)];
end

%% Correction for MT/MST border
if ismember('MST', brain_area)
    if isfield(prs,'MT_border')
        if ~isempty(prs.MT_border)
            
            
            MST_counter = 0;
            for a = 1:size(brain_area, 2)
                if strcmp(brain_area{a},'MST')
                    MST_counter = MST_counter + 1;
                    if MST_counter >= prs.MT_border
                        brain_area{a} = 'MT';
                    end
                end   
            end
           
            
            
        end
    end
end


end

