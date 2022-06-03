function lead_electrode = plotCoherence(lfps,pop_lfps,prs,monk_id,session_id,save_fld)
% lfps : object of class lfps, fields are channel_id, electrode_id...
% pop_lfps : stats variable containing the coherence info
nlfps = length(lfps);
theta_peak = prs.lfp_theta_peak;
beta_peak = prs.lfp_beta_peak;

C = pop_lfps.stats.crosslfp.coher;
f = pop_lfps.stats.crosslfp.freq;

% create a struct containing area ID
brain_area_cell = {lfps.brain_area};
unq = unique(brain_area_cell);
for jj = 1:length(unq)
    
    brain_area.(unq{jj}) = jj;
end

lead_electrode = {};

len_pairs = size(C,2);
vec_ba_i = zeros(len_pairs,1);
vec_ba_j = zeros(len_pairs,1);
electrode_i = zeros(len_pairs,1);
electrode_j = zeros(len_pairs,1);
electrode_type_i = {};

% create pairing vecs based on how C is constructed in
% "coherencyc_unequal_length_trials"
cnt = 1;
for i=1:nlfps
    for j=1:i-1
        electrode_i(cnt) = lfps(i).electrode_id;
        electrode_j(cnt) = lfps(j).electrode_id;
        electrode_type_i{cnt} = lfps(i).electrode_type;
        vec_ba_i(cnt) = brain_area.(lfps(i).brain_area);
        vec_ba_j(cnt) = brain_area.(lfps(j).brain_area);
        cnt = cnt + 1;
    end
end

% id_map = [];
% unique(lfps.
% for i=1:nlfps
%     if lfp
%     brain_area_map
% end


for jj = 1:length(unique(brain_area_cell))
    
    % select the brain area
    filt = (vec_ba_i == jj) & (vec_ba_j == jj);
    ele_filt_i = electrode_i(filt);
    ele_filt_j = electrode_j(filt);
    C_filt = C(:, filt);
    % set the electrode type (assumption: same electrode type x brain
    % area)
    electrode = electrode_type_i{find(filt,1)};
    if contains(electrode,'utah')
        [xloc,yloc] = map_utaharray([],electrode);
    else
        xloc = 1:max(ele_filt_i);
        yloc = ones(size(xloc));
    end
    
    spatial_coher = []; spatial_dist = [];
    
    for ch = 1:length(ele_filt_i)
        spatial_dist(end+1) = sqrt((xloc(ele_filt_i(ch)) - xloc(ele_filt_j(ch)))^2 +...
                                   (yloc(ele_filt_i(ch)) - yloc(ele_filt_j(ch)))^2);
        spatial_coher(end+1,:) = C_filt(:, ch);
        
    end
    
    
    spatial_distances = unique(spatial_dist);
    spatial_coher_mu = zeros(length(spatial_distances), length(f));
    spatial_coher_sem = zeros(length(spatial_distances), length(f));

    for ii = 1:length(spatial_distances)
        spatial_coher_mu(ii,:) = mean(spatial_coher(spatial_dist==spatial_distances(ii),:));
        spatial_coher_sem(ii,:) = std(spatial_coher(spatial_dist==spatial_distances(ii),:))/sqrt(sum(spatial_dist==spatial_distances(ii)));
    end
    
    
    fig = figure; hold on; subplot(1,2,1); hold on;
    suptitle(unq{jj})
    cmap = gray(numel(spatial_distances));
    for i=1:numel(spatial_distances), plot(f,spatial_coher_mu(i,:),'Color',cmap(i,:)); end
    axis([2 80 0.65 1]); xlabel('Frequency (Hz)'); ylabel('Magnitude of coherence between LFPs'); set(gca,'Fontsize',10);
    [~,theta] = min(abs(f - theta_peak)); [~,beta] = min(abs(f - beta_peak)); spatial_multiplier = prs.electrodespacing;
    subplot(1,2,2); hold on; %errorbar(spatial_multiplier*spatial_distances,spatial_coher_mu(:,theta),spatial_coher_sem(:,theta),'ok','MarkerFaceColor','r','Capsize',0);
    hold on; errorbar(spatial_multiplier*spatial_distances,spatial_coher_mu(:,beta),spatial_coher_sem(:,beta),'dk','Capsize',0);
    spatialprs = fmincon(@(x) sum([spatial_coher_mu(:,beta)' - (1 - x(1)*(1 - exp(-(spatial_multiplier*spatial_distances)/x(2))))].^2),[1 1],[],[]);
    plot(linspace(0,5,100),(1 - spatialprs(1)*(1 - exp(-(linspace(0,5,100))/spatialprs(2)))),'k');
    axis([0 5 0.7 1]); xlabel('Distance between electrodes (mm)'); ylabel('Magnitude of coherence between LFPs');
    %                     legend('\theta (8.5 Hz)','\beta (18.5 Hz)'); 
    set(gca,'Fontsize',10);
    saveas(fig,fullfile(save_fld,strcat(sprintf('%s_coher_x_spat_dist_%dm%d.png',unq{jj},monk_id,session_id))))
    if contains(electrode, 'utah') 
        % filter the electrode used in a specific brain area
        electrode_id = [];
        for ii =1:length(lfps)
            if strcmp(lfps(ii).brain_area, unq{jj})
                electrode_id = [electrode_id, lfps(ii).electrode_id];
            end
        end
        % create the heatmap of coherance wrt central electrode
        nflps_filt = length(unique(ele_filt_i));
        
        % extract the xloc and yloc only for the electrodes used in the
        % brain area
        xloc_brain_area = [];
        yloc_brain_area = [];
        electrode_x_brain_area = [];
        for ele = electrode_id
            xloc_brain_area = [xloc_brain_area; xloc(ele)];
            yloc_brain_area = [yloc_brain_area; yloc(ele)];

        end
        
        central_loc = [round(median(xloc_brain_area)),round(median(yloc_brain_area))];    
        % take the cental location
        [~,mmx]=min(abs(xloc_brain_area-central_loc(1)));
        [~,mmy]=min(abs(yloc_brain_area-central_loc(2)));
        
        % select the central ele in original coordingate
        central_ele = electrode_id(find((xloc_brain_area == xloc_brain_area(mmx)) & (yloc_brain_area == yloc_brain_area(mmy)),1));
   
        nele_x = length(unique(xloc_brain_area));
        nele_y = length(unique(yloc_brain_area));

        heatmap_matrix = nan(nele_x,nele_y);
        
        
        for ele = unique(ele_filt_i)'

            if ele == central_ele
                xx = xloc(central_ele);
                yy = yloc(central_ele);
                loc_idx = find((xloc_brain_area == xx) & (yloc_brain_area == yy));
                % start indexing from 1
                xx = xloc_brain_area(loc_idx) - min(xloc_brain_area) + 1;
                yy = yloc_brain_area(loc_idx) - min(yloc_brain_area) + 1;
                heatmap_matrix(xx,yy) = 1;
                continue
            end
            idx = (ele_filt_i ==  central_ele) & (ele_filt_j == ele);
            if ~any(idx)
                idx = (ele_filt_j == central_ele) & (ele_filt_i == ele);
                % electrode loc in origina coordintes
                xx = xloc(ele_filt_i(idx));
                yy = yloc(ele_filt_i(idx));
                
            else
                xx = xloc(ele_filt_j(idx));
                yy = yloc(ele_filt_j(idx));
            end
            % get the idx of the location 
            loc_idx = find((xloc_brain_area == xx) & (yloc_brain_area == yy));
            % start indexing from 1
            xx = xloc_brain_area(loc_idx) - min(xloc_brain_area) + 1;
            yy = yloc_brain_area(loc_idx) - min(yloc_brain_area) + 1;
            heatmap_matrix(xx,yy) = mean(spatial_coher(idx,:));
        end
        fig = figure;
        suptitle(unq{jj})
        % interpolate nan
        for ii = find(isnan(heatmap_matrix))
            [idx1,idx2] = ind2sub(size(heatmap_matrix),ii);
            M = zeros(size(heatmap_matrix));
            M(idx1,idx2) = 1; % location
            heatmap_matrix(idx1,idx2) = mean(heatmap_matrix(conv2(M,[1,1,1;1,0,1;1,1,1],'same')>0));
        end

        [M,I] = sort(max(heatmap_matrix));
        imagesc(imresize(heatmap_matrix,10),[min(min(heatmap_matrix)),M(I(end-1))]); 
        xlabel('y pos')
        ylabel('x pos')
        title('heatmap coherence')
        c = colorbar('Fontsize',14);
        c.Label.String = 'coherence';
        saveas(fig,fullfile(save_fld,strcat(sprintf('%s_coherence_x_heatmap_%dm%d.png',unq{jj},monk_id,session_id))))

        
        %% plot other phase distance
        spatial_phasediff = zeros(nflps_filt,nflps_filt,numel(f)); 
        spatial_dist = zeros(nflps_filt,nflps_filt);
        phi = pop_lfps.stats.crosslfp.phase(:,filt);
        ind2row = @(i,j) min(i,j) + (max(i,j)-1)*(max(i,j)-2)/2; % to read the output of "coherencyc_unequal_length_trials" function from Chronux
        %chan2elec = @(i,j) [electrode_id(i) electrode_id(j)];
        cnt = 1;
        for i=1:nflps_filt
            for j=1:nflps_filt
                if i==j, spatial_phasediff(i,j,:) = zeros(numel(f),1);  % zero phase-lag with itself
%                 elseif i>j, spatial_phasediff(i,j,:) = phi(:,ind2row(i,j));
%                 elseif i<j, spatial_phasediff(i,j,:) = -phi(:,ind2row(i,j)); %phase x rel. to y = -phase y rel. to x
%                 end

                elseif j < i
                    spatial_phasediff(i,j,:) = -phi(:,cnt);
                    spatial_phasediff(j,i,:) = phi(:,cnt);
                    ele_i = electrode_i(cnt);
                    ele_j = electrode_j(cnt);
                    spatial_dist(i,j) = sqrt((xloc(ele_i)-xloc(ele_j))^2 + (yloc(ele_i)-yloc(ele_j))^2);
                    
                    cnt = cnt + 1;
                end
                  
            end
        end
        spatial_dist = spatial_dist + spatial_dist';
        %% theta phase map
        %             [~,theta] = min(abs(f - theta_peak)); 
        %% beta phase map
        [~,beta] = min(abs(f - beta_peak));
        

        phasediffs = mean(squeeze(spatial_phasediff(:,:,beta)),2);
        [phasediffs_sorted,phaseorder] = sort(phasediffs);
        electrode_id_sorted = electrode_id(phaseorder);
        lead_electrode(end+1,:) = {unq{jj},electrode_id_sorted(1)};
        cmap = cool(nflps_filt);
        fig = figure; hold on;
        suptitle(unq{jj})
        subplot(1,2,1); hold on;
        for i=1:nflps_filt
            plot(f,(180/pi)*squeeze(spatial_phasediff(phaseorder(1),phaseorder(i),:)),'Color',cmap(i,:)); % plot leader vs everyone else
        end
        axis([1 80 -50 50]); hline(0,'k'); xlabel('Frequency (Hz)'); ylabel('Phase of coherence between LFPs (deg)');
        subplot(1,2,2); hold on;
        
        centered_xloc = xloc_brain_area - min(xloc_brain_area) + 1;
        centered_yloc = yloc_brain_area - min(yloc_brain_area) + 1;

        zloc = nan(max(centered_xloc),max(centered_yloc));
        for i = 1:nflps_filt
            xx = find(electrode_id == electrode_id_sorted(i),1);
            zloc(centered_xloc(xx),centered_yloc(xx)) = i;
            plot(centered_xloc(xx),centered_yloc(xx),'o','Color',cmap(i,:),'MarkerFaceColor',cmap(i,:));
        end
        % fill in the edges of the array (for smoothing)
        zloc(1,1) = 0.5*(zloc(1,2) + zloc(2,1)); zloc(1,end) = 0.5*(zloc(1,end-1) + zloc(2,end));
        zloc(end,1) = 0.5*(zloc(end-1,1) + zloc(end,2)); zloc(end,end) = 0.5*(zloc(end,end-1) + zloc(end-1,end));
        
        % interpolate nan
        for ii = find(isnan(zloc))
            [idx1,idx2] = ind2sub(size(zloc),ii);
            M = zeros(size(zloc));
            M(idx1,idx2) = 1; % location
            zloc(idx1,idx2) = mean(zloc(conv2(M,[1,1,1;1,0,1;1,1,1],'same')>0));
        end
        
        axis([0 11 0 11]); axis off; colormap(cool); colorbar;
        maxtimelag = ((phasediffs_sorted(end) - phasediffs_sorted(1))/(2*pi))*(1/beta_peak)*1e3;
        colorbar('Ticks',[0,1],'TickLabels',{[num2str(0) ' ms'],[num2str(round(maxtimelag*10)/10) ' ms']},'Fontsize',14);
        %             DrawPhaseArrows([xloc yloc],electrode_id_sorted);
        saveas(fig,fullfile(save_fld,strcat(sprintf('%s_delta_phase_x_spat_dist_%dm%d.png',unq{jj},monk_id,session_id))))

        fig = figure; imagesc(imresize(zloc,10)); colormap(cool);
        suptitle(unq{jj})
        c = colorbar('Fontsize',14);
        c.Label.String = 'delta phase - beta';
        saveas(fig,fullfile(save_fld,strcat(sprintf('%s_delta_phase_heatmap_%dm%d.png',unq{jj},monk_id,session_id))))

        %% velocity
        spatial_phasediff_beta = squeeze(spatial_phasediff(:,:,beta)); spatial_phasediff_beta = spatial_phasediff_beta(:);
        spatial_dists = unique(spatial_dist);
       
        for i=2:numel(spatial_dists)
           
            median_phasediff_beta(i) = median(abs(spatial_phasediff_beta(spatial_dist(:)==spatial_dists(i)))); 
            sem_phasediff_beta(i) = 1.2533*std(abs(spatial_phasediff_beta(spatial_dist(:)==spatial_dists(i))))/sqrt(sum(spatial_dist(:)==spatial_dists(i)));
            [F,X] = ecdf(abs(spatial_phasediff_beta(spatial_dist(:)==spatial_dists(i)))); 
            if length(unique(X)) < 3 
                ecdf_phasediff_beta(i,:) = nan;
                continue
            end
            [X,indx] = unique(X); F = F(indx);
            ecdf_phasediff_beta(i,:) = interp1(X,F,linspace(0,max(X),101));
        end
        %figure; surf(spatial_dists',linspace(0,0.5,101)',ecdf_phasediff_beta'); colormap(goodcolormap('bwr',100)); set(gca,'YDir','reverse')
        %hold on; plot(spatial_dists,median_phasediff_beta,'dk'); axis([0 9.5 0 0.2906]); 
    end
    

end

end