%% add units
function AddUnits(this,prs)

cd(prs.filepath_neur); % cd(fullfile(prs.filepath_neur,'Sorted'));
% determine type of electrode
linearprobe_type = []; utaharray_type = [];
for k=1:length(prs.electrode_type)
    linearprobe_type = [linearprobe_type find(cellfun(@(electrode_type) strcmp(prs.electrode_type{k},electrode_type), prs.linearprobe.types),1)];
    utaharray_type = [utaharray_type find(cellfun(@(electrode_type) strcmp(prs.electrode_type{k},electrode_type), prs.utaharray.types),1)];
end

if isfield(prs, 'isRipple')
    if prs.isRipple
        
        prs.neur_filetype = 'nev';
        % Directory changed back to array data
        SortedFolder = dir(prs.filepath_neur);
        ind = arrayfun(@(x) contains(x.name,'Sorted'),SortedFolder);
        SortedFolder = {SortedFolder(ind).name};
        Nfolders = numel(SortedFolder);

        [electrode_id, brain_area, ~, electrode_type] = MapChannel2Electrode_Ripple_Units(prs);

        Nch = numel(electrode_id);
        ch_per_fld = Nch/Nfolders;
        ch_ind = 0:ch_per_fld:Nch;
        
        sua = []; mua = [];
        for k = 1:numel(SortedFolder)
            cd(fullfile(prs.filepath_neur,SortedFolder{k}));
%             brain_area = prs.area{strcmp(prs.electrode_type,prs.linearprobe.types{linearprobe_type})};
            [sua_tmp, mua_tmp] = GetUnits_phy_Ripple('spike_times.npy', 'spike_clusters.npy', 'cluster_info.tsv',electrode_id(ch_ind(k)+1:ch_ind(k+1)), electrode_type(ch_ind(k)+1:ch_ind(k+1)), brain_area(ch_ind(k)+1:ch_ind(k+1)));
            % correction for breaking data in parts
            for jj = 1:numel(sua_tmp)
                sua_tmp(jj).channel_id = sua_tmp(jj).channel_id + ch_ind(k);
                sua_tmp(jj).cluster_id = sua_tmp(jj).cluster_id + ch_ind(k);
            end
            for jj = 1:numel(mua_tmp)
                mua_tmp(jj).channel_id = mua_tmp(jj).channel_id + ch_ind(k);
                mua_tmp(jj).cluster_id = mua_tmp(jj).cluster_id + ch_ind(k);
            end
            sua = [sua sua_tmp];
            mua = [mua mua_tmp];
        end
        
        %brain_area = prs.area{strcmp(prs.electrode_type,prs.utaharray.types{utaharray_type})};
        
        % Get Events
        cd(prs.filepath_neur); 
        file_nev=dir('*.nev'); 
        fprintf(['... reading events from ' file_nev.name '\n']);
        [events_nev,prs] = GetEvents_nev(file_nev.name,prs); % requires package from Blackrock Microsystems: https://github.com/BlackrockMicrosystems/NPMK
        t_end_nev = [events_nev.t_end];
        t_beg_nev = [events_nev.t_beg];
        if length(t_beg_nev) > length(t_end_nev)
            if all(t_end_nev - t_beg_nev(1:length(t_end_nev)) > 0)
                events_nev.t_beg = t_beg_nev(1:length(t_end_nev));
            elseif all(t_end_nev - t_beg_nev(length(t_beg_nev)-length(t_end_nev)+1:end) > 0)
                events_nev.t_beg = t_beg_nev(length(t_beg_nev)-length(t_end_nev)+1:end);
            end
        end
        %prs.fs_spk = 30000;
        if length(this.behaviours.trials)~=length(events_nev.t_end)
            events_nev = FixEvents_nev(events_nev,this.behaviours.trials);
        end
        if length(this.behaviours.trials)==length(events_nev.t_end)
            if ~isempty(sua)
                for i=1:length(sua)
                    %fetch singleunit
                    this.units(end+1) = unit('singleunit',sua(i),prs.fs_spk);
%                     if strcmp({prs.utaharray.types{utaharray_type}},'utah96')
%                         this.units(end).brain_area = brain_area;
%                     else
%                         this.units(end).brain_area = prs.MapDualArray2BrainArea(brain_area, this.units(end).electrode_id);
%                     end
                    this.units(end).AddTrials(sua(i).tspk,events_nev,this.behaviours,prs);
                end
            end
            if ~isempty(mua)
                for i=1:length(mua)
                    %fetch multiunit
                    this.units(end+1) = unit('multiunit',mua(i),prs.fs_spk);
%                     if strcmp(prs.utaharray.types{utaharray_type},'utah96')
%                         this.units(end).brain_area = brain_area;
%                     else
%                         this.units(end).brain_area = prs.MapDualArray2BrainArea(brain_area, this.units(end).electrode_id);
%                     end
                    this.units(end).AddTrials(mua(i).tspk,events_nev,this.behaviours,prs);
                end
            end
        else
            fprintf('Cannot segment spikes: Trial counts in smr and nev files do not match \n');
            fprintf(['Trial end events: NEV file - ' num2str(length(events_nev.t_end)) ...
                ' , SMR file - ' num2str(length(this.behaviours.trials)) '\n']);
            fprintf('Debug and try again! \n');
        end
    else
        if ~isempty(linearprobe_type) % assume linearprobe is recorded using Plexon
            % Need to change directory to Plexon data
            try
                cd(fullfile(prs.filepath_neur,'PLEXON FILES','Sorted'));
            catch
                cd(fullfile(prs.filepath_neur,'Sorted'));
            end
            brain_area = prs.area{strcmp(prs.electrode_type,prs.linearprobe.types{linearprobe_type})};
            file_ead=dir('*_ead.plx'); prs.neur_filetype = 'plx';
            fprintf(['... reading ' file_ead.name '\n']);
            if ~isempty(file_ead)
                [events_plx] = GetEvents_plx(file_ead.name);
            else
                file_plx = dir([prs.filepath_neur, '*.plx']);
                [events_plx] = GetEvents_plx(fullfile(prs.filepath_neur, file_plx.name));
            end
            
            prs.fs_spk = 20000;
            file_plx=dir('*.dat');
            fprintf(['... reading ' file_plx.name '\n']);
            [sua, mua] = GetUnits_phy('spike_times.npy', 'spike_clusters.npy', 'cluster_info.tsv',prs.linearprobe.types{linearprobe_type});
            
            if ~isempty(sua)
                for i=1:length(sua)
                    %fetch singleunit
                    if numel(sua(i).tspk)/numel(events_plx.t_beg) > prs.minspk
                        this.units(end+1) = unit('singleunit',sua(i),prs.fs_spk);
                        this.units(end).brain_area = brain_area;
                        this.units(end).AddTrials(sua(i).tspk,events_plx,this.behaviours,prs);
                    end
                end
            end
            if ~isempty(mua)
                for i=1:length(mua)
                    %fetch multiunit
                    if numel(mua(i).tspk)/numel(events_plx.t_beg) > prs.minspk
                        this.units(end+1) = unit('multiunit',mua(i),prs.fs_spk);
                        this.units(end).brain_area = brain_area;
                        this.units(end).AddTrials(mua(i).tspk,events_plx,this.behaviours,prs);
                    end
                end
            end
        end
        
        
        if ~isempty(utaharray_type) % assume utaharray is recorded using Cereplex or Ripple
            % Directory changed back to array data
            SortedFolder = dir(prs.filepath_neur);
            ind = arrayfun(@(x) contains(x.name,'Sorted'),SortedFolder);
            SortedFolder = {SortedFolder(ind).name};
            
            sua = []; mua = [];
            for k = 1:numel(SortedFolder)
                cd(fullfile(prs.filepath_neur,SortedFolder{k}));
                [sua_tmp, mua_tmp] = GetUnits_phy('spike_times.npy', 'spike_clusters.npy', 'cluster_info.tsv',prs.utaharray.types{utaharray_type});%     [sua, mua] = GetUnits_phy('spike_times.npy', 'spike_clusters.npy', 'cluster_group.tsv','cluster_location.xls',prs.utaharray.types{utaharray_type}); % requires npy-matlab package: https://github.com/kwikteam/npy-matlab
                sua = [sua sua_tmp];
                mua = [mua mua_tmp];
            end
            
            % Get events
            cd(prs.filepath_neur);
            file_nev=dir('*.nev'); prs.neur_filetype = 'nev';
            brain_area = prs.area{strcmp(prs.electrode_type,{prs.utaharray.types{utaharray_type}})};
            fprintf(['... reading events from ' file_nev.name '\n']);
            [events_nev,prs] = GetEvents_nev(file_nev.name,prs); % requires package from Blackrock Microsystems: https://github.com/BlackrockMicrosystems/NPMK
            t_end_nev = [events_nev.t_end];
            t_beg_nev = [events_nev.t_beg];
            if length(t_beg_nev) > length(t_end_nev)
                if all(t_end_nev - t_beg_nev(1:length(t_end_nev)) > 0)
                    events_nev.t_beg = t_beg_nev(1:length(t_end_nev));
                elseif all(t_end_nev - t_beg_nev(length(t_beg_nev)-length(t_end_nev)+1:end) > 0)
                    events_nev.t_beg = t_beg_nev(length(t_beg_nev)-length(t_end_nev)+1:end);
                end
            end
            prs.fs_spk = 30000;
            if length(this.behaviours.trials)~=length(events_nev.t_end)
                events_nev = FixEvents_nev(events_nev,this.behaviours.trials);
            end
            
            if length(this.behaviours.trials)==length(events_nev.t_end)
                if ~isempty(sua)
                    for i=1:length(sua)
                        %fetch singleunit
                        this.units(end+1) = unit('singleunit',sua(i),prs.fs_spk);
                        if strcmp(prs.utaharray.types{utaharray_type},'utah96')
                            this.units(end).brain_area = brain_area;
                        else
                            this.units(end).brain_area = prs.MapDualArray2BrainArea(brain_area, this.units(end).electrode_id);
                        end
                        this.units(end).AddTrials(sua(i).tspk,events_nev,this.behaviours,prs);
                    end
                end
                if ~isempty(mua)
                    for i=1:length(mua)
                        %fetch multiunit
                        this.units(end+1) = unit('multiunit',mua(i),prs.fs_spk);
                        if strcmp(prs.utaharray.types{utaharray_type},'utah96')
                            this.units(end).brain_area = brain_area;
                        else
                            this.units(end).brain_area = prs.MapDualArray2BrainArea(brain_area, this.units(end).electrode_id);
                        end
                        this.units(end).AddTrials(mua(i).tspk,events_nev,this.behaviours,prs);
                    end
                end
            else
                fprintf('Cannot segment spikes: Trial counts in smr and nev files do not match \n');
                fprintf(['Trial end events: NEV file - ' num2str(length(events_nev.t_end)) ...
                    ' , SMR file - ' num2str(length(this.behaviours.trials)) '\n']);
                fprintf('Debug and try again! \n');
            end
        else
            fprintf('No cereplex neural data files in the specified path \n');
        end
    end
end
