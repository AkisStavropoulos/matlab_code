function [ch_id, electrode_id, brain_area] = MapChannel2Electrode2(total_channels, electrode_type, area)
        % total_channels should be an array with number of electrodes per
    % channels
        % electrode type should tell us the electrode type
        % area, same
        % prs, any other info we may need
               
        n_areas = size(total_channels, 1); 
        for i = 1:n_areas
            electrode = electrode_type{i};
            [temp_chnl_indx,temp_elec_indx] = MapChannel2Electrode(electrode);
            temp_area = cell(1, length(temp_chnl_indx));
            temp_area(:) = area(i);
        end
        
        ch_id = temp_chnl_indx; 
        electrode_id = temp_elec_indx;
        brain_area = temp_area; 
        
        
        
end