classdef unit < handle
    %%
    properties
        cluster_id
        channel_id
        electrode_id
        uQ
        cR
        isiV
        electrode_type
        brain_area
        spkwf
        spkwidth
        type
        trials
        stats
    end
    %%
    methods
        %% class constructor
        function this = unit(unittype,unit,Fs)
            this.cluster_id = unit.cluster_id;
            this.channel_id = unit.channel_id;
            try
                this.uQ = unit.uQ;
                this.cR= unit.cR;
                this.isiV = unit.isiV;
            catch
                this.uQ = nan;
                this.cR= nan;
                this.isiV = nan;
            end
            this.electrode_id = unit.electrode_id;
            this.electrode_type = unit.electrode_type;
            try
                this.brain_area = unit.brain_area;
            catch
            end
            this.spkwf = unit.spkwf; %mean spike-waveform;
            this.spkwidth = Compute_SpikeWidth(unit.spkwf,Fs);
            this.type = unittype;
        end
    end
end