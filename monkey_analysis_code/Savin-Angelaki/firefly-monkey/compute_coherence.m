%% list all file in a folder
%listing = dir('/Volumes/WD Edo/firefly_analysis/DATASET/PPC+PFC+MST/');
listing = dir('D:\Savin-Angelaki\saved\');
expression = '^m[0-9]+s[0-9]+.mat$';
for i = 1:length(listing)
    % check that file names matches with the regexp
    if isempty(regexp(listing(i).name,expression))
        continue
    end

    splt = split(listing(i).name,'s');
    monkey_id = str2num(splt{1}(2:end));
    
    strng = splt{2};
    splt = split(strng,'.');
    session_id = str2num(splt{1});
    % load, set flags and analyze
    load(fullfile(listing(i).folder,listing(i).name))
    prs.fitGAM_coupled = 0;
    prs.compute_canoncorr = 0;
    prs.regress_popreadout = 0;
	prs.simulate_population = 0;
	prs.compute_coherencyLFP = 1;
	prs.corr_neuronbehverr = 0;
    stats_lfp = AnalysePopulation(lfps,trials_behv,behv_stats,lfps,prs);

    pop_lfps.stats = stats_lfp;
    lead_electrode = plotCoherence(lfps,pop_lfps,prs,monkey_id,session_id,save_lfd);
    %save
    save(fullfile('/Volumes/WD Edo/firefly_analysis/LFP_band/processed_data/LFP_coherence',strcat('LFP_coherence_',listing(i).name)),'lead_electrode','stats_lfp')
    
    % plot 
    pop_lfps.stats = stats_lfp;
    
end