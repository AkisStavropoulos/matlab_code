function [sua, mua] = GetUnits_phy(f_spiketimes, f_spikeclusters, f_clusterinfo, electrode_type)

sua = []; mua = [];

[~,electrode_id] = MapChannel2Electrode(electrode_type);
spiketimes = readNPY(f_spiketimes);
cluster_ids = readNPY(f_spikeclusters);
if exist(f_clusterinfo,'file') % remove if clause once we have these files for all recording sessions
    cluster_info = tdfread('cluster_info.tsv');
end
try
    load 'waveForms.mat'; load 'QualityMetr.mat';
catch
end
sua_indx = find(strtrim(string(cluster_info.group))=='good');
for i = 1:length(sua_indx)
    sua(i).tspk = spiketimes(cluster_ids == cluster_info.id(sua_indx(i)));
    sua(i).cluster_id = cluster_info.id(sua_indx(i));
    if ~isempty(cluster_info)
        try
            sua(i).channel_id = cluster_info.channel(sua_indx(i))+1; %K2, Phy2 returns channels from 0
        catch
            sua(i).channel_id = cluster_info.ch(sua_indx(i))+1; %K2, Phy2 returns channels from 0
        end
        sua(i).electrode_id = electrode_id(sua(i).channel_id);
        sua(i).electrode_type = electrode_type;
        %         sua(i).spkwf = squeeze(mean(waveForms(str2double({clusters.id}) == str2double(clusters(sua_indx(i)).id),:,:),2));
        try
            sua(i).spkwf = waveFormsMean(sua_indx(i),:);
            sua(i).uQ = uQ(sua_indx(i)); sua(i).cR = cR(sua_indx(i)); sua(i).isiV = isiV(sua_indx(i));
        catch
        end
    else
        sua(i).channel_id = [];
        sua(i).electrode_id = [];
        sua(i).electrode_type = electrode_type;
        sua(i).spkwf = [];
    end
end

mua_indx = find(strtrim(string(cluster_info.group))=='mua');
for i = 1:length(mua_indx)
    mua(i).tspk = spiketimes(cluster_ids == cluster_info.id(mua_indx(i)));
    mua(i).cluster_id = cluster_info.id(mua_indx(i));
    if ~isempty(cluster_info)
        try
            mua(i).channel_id = cluster_info.channel(mua_indx(i))+1; %K2, Phy2 returns channels from 0
        catch
            mua(i).channel_id = cluster_info.ch(mua_indx(i))+1; %K2, Phy2 returns channels from 0
        end
        mua(i).electrode_id = electrode_id(mua(i).channel_id);
        mua(i).electrode_type = electrode_type;
        %         mua(i).spkwf = squeeze(mean(waveForms(str2double({clusters.id}) == str2double(clusters(mua_indx(i)).id),:,:),2));
        try
            mua(i).spkwf = waveFormsMean(mua_indx(i),:);
            mua(i).uQ = uQ(mua_indx(i)); mua(i).cR = cR(mua_indx(i)); mua(i).isiV = isiV(mua_indx(i));
        catch
        end
    else
        mua(i).channel_id = [];
        mua(i).electrode_id = [];
        mua(i).electrode_type = electrode_type;
        mua(i).spkwf = [];
    end
end