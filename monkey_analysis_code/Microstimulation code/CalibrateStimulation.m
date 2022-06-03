%% specify which eyes were tracked
prs.eyechannels = [2 2]; % example: use [1 1] when tracking both eyes with coil, [0 2] for only right-eye with eye-tracker etc.
prs.saccade_thresh = 50; % saccade threshold for eye velocity (deg/s)
prs.blink_thresh = 50; %deg
prs.nanpadding = 5; % samples on either side
prs.filtwidth = 2; % samples
prs.filtsize = 2*prs.filtwidth; % samples
prs.min_intersaccade = 0.1; % (s) minimum inter-saccade interval
prs.post_saccade = 0.2; % time after saccade onset to plot
prs.saccade_duration = 0.15; % duration of saccade to compute direction

%% read log file
file_log = dir('*.log');
fid = fopen(file_log.name, 'r');
eof=0; newline = 'nothingnew'; count=0;
while ~strcmp(newline(1:6),'Monkey'), newline = fgetl(fid); end
stimlog.animalname = newline(14:end);
while ~strcmp(newline(1:12),'DataFileName'), newline = fgetl(fid); end
stimlog.datafilename = newline(15:end);
while ~strcmp(newline(1:14),'Pulse Duration'), newline = fgetl(fid); end
stimlog.pulseduration = str2double(newline(21:end));
while ~strcmp(newline(1:15),'Pulse Frequency'), newline = fgetl(fid); end
stimlog.stimfrequency = str2double(newline(17:end));
while ~strcmp(newline(1:14),'Train Duration'), newline = fgetl(fid); end
stimlog.trainduration = str2double(newline(22:end));
while ~strcmp(newline(1:14),'Train Interval'), newline = fgetl(fid); end
stimlog.intertraininterval = str2double(newline(22:end));
while ~strcmp(newline(1:14),'Pulse Polarity'), newline = fgetl(fid); end
stimlog.pulsepolarity = str2double(newline(17:end));
while ~strcmp(newline(1:10),'Repetition'), newline = fgetl(fid); end
stimlog.nreps = str2double(newline(13:end));
ntrls = 0;
while newline ~= -1
    while ~strcmp(newline(1:9),'Trial Num')
        newline = fgetl(fid);
        if newline == -1, break; end
    end
    ntrls = ntrls + 1;
    newline = fgetl(fid); if newline == -1, break; end
    stimlog.pulseamplitude(ntrls) = str2double(newline(27:end));
end
fclose(fid);

%% load smr file
file_smr = dir('*.smr');
data = ImportSMR(file_smr.name);

%% channel titles
nch = length(data);
ch_title = cell(1,nch);
hdr = {data.hdr};
for i=1:nch
    if ~isempty(hdr{i})
        ch_title{i} = hdr{i}.title;
    else
        ch_title{i} = 'nan';
    end
end

%% channel number
chnum.pulse = find(strcmp(ch_title,'Pulse'));
chnum.yle = find(strcmp(ch_title,'LDy')); chnum.zle = find(strcmp(ch_title,'LDz'));
chnum.yre = find(strcmp(ch_title,'RDy')); chnum.zre = find(strcmp(ch_title,'RDz'));
chnum.marker = find(strcmp(ch_title,'marker'));

%% parameters for converting from volts
scaling.t = data(chnum.marker).hdr.tim.Scale*data(chnum.marker).hdr.tim.Units;
scaling.pulse = data(chnum.pulse).hdr.adc.Scale; offset.pulse = data(chnum.pulse).hdr.adc.DC;
scaling.yle = data(chnum.yle).hdr.adc.Scale; offset.yle = data(chnum.yle).hdr.adc.DC;
scaling.yre = data(chnum.yre).hdr.adc.Scale; offset.yre = data(chnum.yre).hdr.adc.DC; 
scaling.zle = data(chnum.zle).hdr.adc.Scale; offset.zle = data(chnum.zle).hdr.adc.DC; 
scaling.zre = data(chnum.zre).hdr.adc.Scale; offset.zre = data(chnum.zre).hdr.adc.DC;

%% read channel data
chnames = fieldnames(chnum); MAX_LENGTH = inf; dt = [];
for i=1:length(chnames)-1
    ch.(chnames{i}) = double(data(chnum.(chnames{i})).imp.adc)*scaling.(chnames{i}) + offset.(chnames{i});
    dt.(chnames{i}) = prod(data(chnum.(chnames{i})).hdr.adc.SampleInterval);
    ts.(chnames{i}) = dt.(chnames{i}):dt.(chnames{i}):length(ch.(chnames{i}))*dt.(chnames{i});
    MAX_LENGTH = min(length(ch.(chnames{i})),MAX_LENGTH);
end
ch.yle = ch.yle(1:MAX_LENGTH); ts.yle = ts.yle(1:MAX_LENGTH);
ch.yre = ch.yre(1:MAX_LENGTH); ts.yre = ts.yre(1:MAX_LENGTH);
ch.zle = ch.zle(1:MAX_LENGTH); ts.zle = ts.zle(1:MAX_LENGTH);
ch.zre = ch.zre(1:MAX_LENGTH); ts.zre = ts.zre(1:MAX_LENGTH);

%% event markers
markers = data(chnum.marker).imp.mrk(:,1);
t.events = double(data(chnum.marker).imp.tim)*scaling.t;
t.stimbeg = t.events(markers == 32); %
t.stimend = t.stimbeg + stimlog.trainduration*1e-3; % ms to s

%% remove eye blinks if using eye tracker and smooth
X = [ch.zle ch.zre ch.yle ch.yre];
X = ReplaceWithNans(X, prs.blink_thresh, prs.nanpadding);
ch.zle = X(:,1); ch.zre = X(:,2); ch.yle = X(:,3); ch.yre = X(:,4);
sig = 1*prs.filtwidth; %filter width
sz = 1*prs.filtsize; %filter size
t2 = linspace(-sz/2, sz/2, sz);
h = exp(-t2.^2/(2*sig^2));
h = h/sum(h); % normalise filter to ensure area under the graph of the data is not altered
ch.zle = conv(ch.zle,h,'same'); ch.zre = conv(ch.zre,h,'same'); 
ch.yle = conv(ch.yle,h,'same'); ch.yre = conv(ch.yre,h,'same');

%% take derivative of eye position to get eye velocity
if prs.eyechannels(2) == 0 % use the eye with a working eye coil
    dze = diff(ch.zle);
    dye = diff(ch.yle);
elseif prs.eyechannels(1) == 0
    dze = diff(ch.zre);
    dye = diff(ch.yre);
else
    dze = [diff(ch.zle) diff(ch.zre)];
    dye = [diff(ch.yle) diff(ch.yre)];
end
eyevel_ver = dze/dt.zle; 
eyevel_hor = dye/dt.yle;
de = sqrt(dze.^2 + dye.^2); % speed of eye movement
sig = 5*prs.filtwidth; %filter width
sz = 5*prs.filtsize; %filter size
t2 = linspace(-sz/2, sz/2, sz);
h = exp(-t2.^2/(2*sig^2));
h = h/sum(h); % normalise filter to ensure area under the graph of the data is not altered
for i=1:size(de,2), de(:,i) = conv(de(:,i),h,'same'); end

%% apply threshold on eye speed
saccade_thresh = prs.saccade_thresh;
thresh = saccade_thresh*(dt.zle); % threshold in units of deg/sample
if size(de,2) == 1, indx_thresh = (de>thresh & de<20*thresh); else, indx_thresh = all((de>thresh & de<20*thresh)')'; end % upper limit to remove blink
dindx_thresh = diff(indx_thresh); 
t_saccade = find(dindx_thresh>0)*(dt.zle);
% remove duplicates by applying a saccade refractory period
min_isi = prs.min_intersaccade;
t_saccade(diff(t_saccade)<min_isi) = [];
t.saccade = t_saccade;

%% interpolate nans
% nanx = isnan(ch.zle); t1 = 1:numel(ch.zle); ch.zle(nanx) = interp1(t1(~nanx), ch.zle(~nanx), t1(nanx), 'pchip');
% nanx = isnan(ch.zre); t1 = 1:numel(ch.zle); ch.zre(nanx) = interp1(t1(~nanx), ch.zre(~nanx), t1(nanx), 'pchip');
% nanx = isnan(ch.yle); t1 = 1:numel(ch.yle); ch.yle(nanx) = interp1(t1(~nanx), ch.yle(~nanx), t1(nanx), 'pchip');
% nanx = isnan(ch.yre); t1 = 1:numel(ch.yre); ch.yre(nanx) = interp1(t1(~nanx), ch.yre(~nanx), t1(nanx), 'pchip');

%% plot stimulation and eye movements
k = 21; % skip first and last k points when plotting because of edge artifacts
figure; hold on;
subplot(3,1,1); hold on; plot(ts.pulse,ch.pulse);
xlabel('Time (s)'); ylabel({'Current' ; 'amplitude (V)'});
subplot(3,1,2); hold on; plot(ts.yle(k:end-k),ch.yle(k:end-k));
xlabel('Time (s)'); ylabel({'Eye horizontal (deg)'});
subplot(3,1,2); hold on; plot(ts.yre(k:end-k),ch.yre(k:end-k));
ylim([-40 40]);
subplot(3,1,3); hold on; plot(ts.zle(k:end-k),ch.zle(k:end-k));
xlabel('Time (s)'); ylabel({'Eye Vertical (deg)'});
subplot(3,1,3); hold on; plot(ts.zre(k:end-k),ch.zre(k:end-k));
ylim([-40 40]);
saveas(gcf,['Fig_' num2str(i), '.epsc']);

%% plot stimulation and eye movements
k = 21; % skip first and last k points when plotting because of edge artifacts
figure; hold on;
subplot(3,1,1); hold on; plot(ts.pulse,ch.pulse);
xlabel('Time (s)'); ylabel({'Current' ; 'amplitude (V)'});
subplot(3,1,2); hold on; plot(ts.yle(k:end-k),ch.yle(k:end-k));
xlabel('Time (s)'); ylabel({'Eye horizontal (deg)'});
subplot(3,1,2); hold on; plot(ts.yre(k:end-k),ch.yre(k:end-k)); vline(t.stimbeg,'--k');
ylim([-40 40]);
subplot(3,1,3); hold on; plot(ts.zle(k:end-k),ch.zle(k:end-k));
xlabel('Time (s)'); ylabel({'Eye Vertical (deg)'});
subplot(3,1,3); hold on; plot(ts.zre(k:end-k),ch.zre(k:end-k)); vline(t.stimbeg,'--k');
ylim([-40 40]);

%% compute saccade direction
pulseamplitudes = unique(stimlog.pulseamplitude); namplitudes = length(pulseamplitudes);
nsaccades = zeros(1,namplitudes);
nt_saccade = round(prs.saccade_duration/dt.yle);
for i=1:namplitudes
    % saccade count
    nsaccades(i) = sum(t.saccade > t.stimbeg(stimlog.pulseamplitude == pulseamplitudes(i))' & ...
        t.saccade < t.stimend(stimlog.pulseamplitude == pulseamplitudes(i))', 'all');
    % saccade times
    saccade_time{i} = t.saccade(logical(sum(t.saccade > t.stimbeg(stimlog.pulseamplitude == pulseamplitudes(i))' & ...
        t.saccade < t.stimend(stimlog.pulseamplitude == pulseamplitudes(i))',2)));
    % stimulation times
    stim_onset{i} = t.stimbeg(stimlog.pulseamplitude == pulseamplitudes(i));
    % saccade direction aligned to saccade onset
    for j=1:nsaccades(i)
        post_saccade_period = (ts.yle > saccade_time{i}(j) & ts.yle < saccade_time{i}(j)+prs.post_saccade);
        horpos{i}{j} = 0.5*(ch.yle(post_saccade_period) + ch.yre(post_saccade_period)); horpos{i}{j} = horpos{i}{j} - horpos{i}{j}(1);
        verpos{i}{j} = 0.5*(ch.zle(post_saccade_period) + ch.zre(post_saccade_period)); verpos{i}{j} = verpos{i}{j} - verpos{i}{j}(1);
        saccade_dir{i}(:,j) = [horpos{i}{j}(nt_saccade) - horpos{i}{j}(1) ;  verpos{i}{j}(nt_saccade) - verpos{i}{j}(1)];
    end
    minlgth = min(cellfun(@(x) numel(x),horpos{i}));
    horpos{i} = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],horpos{i},'un',0));
    verpos{i} = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],verpos{i},'un',0));
    
    % saccade direction aligned to stimulation onset
    for j=1:stimlog.nreps
        post_stim_period = (ts.yle > stim_onset{i}(j) & ts.yle < stim_onset{i}(j)+stimlog.trainduration*1e-3); % doesn't work for extremely short trains as saccade is not completed
        post_stim_period = (ts.yle > stim_onset{i}(j) & ts.yle < stim_onset{i}(j)+prs.post_saccade);
        left_horpos_stimevoked{i}{j} = (ch.yle(post_stim_period)); left_horpos_stimevoked{i}{j} = left_horpos_stimevoked{i}{j} - left_horpos_stimevoked{i}{j}(1);
        right_horpos_stimevoked{i}{j} = (ch.yre(post_stim_period)); right_horpos_stimevoked{i}{j} = right_horpos_stimevoked{i}{j} - right_horpos_stimevoked{i}{j}(1);
        left_verpos_stimevoked{i}{j} = (ch.zle(post_stim_period)); left_verpos_stimevoked{i}{j} = left_verpos_stimevoked{i}{j} - left_verpos_stimevoked{i}{j}(1);
        right_verpos_stimevoked{i}{j} = (ch.zre(post_stim_period)); right_verpos_stimevoked{i}{j} = right_verpos_stimevoked{i}{j} - right_verpos_stimevoked{i}{j}(1);
        hordir = 0.5*((left_horpos_stimevoked{i}{j}(nt_saccade) - left_horpos_stimevoked{i}{j}(1)) + (right_horpos_stimevoked{i}{j}(nt_saccade) - right_horpos_stimevoked{i}{j}(1)));
        verdir = 0.5*((left_verpos_stimevoked{i}{j}(nt_saccade) - left_verpos_stimevoked{i}{j}(1)) + (right_verpos_stimevoked{i}{j}(nt_saccade) - right_verpos_stimevoked{i}{j}(1)));
        saccade_dir_stimevoked{i}(:,j) = [hordir ; verdir];
    end
    minlgth = min(cellfun(@(x) numel(x),left_horpos_stimevoked{i}));
    left_horpos_stimevoked{i} = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],left_horpos_stimevoked{i},'un',0));
    right_horpos_stimevoked{i} = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],right_horpos_stimevoked{i},'un',0));
    left_verpos_stimevoked{i} = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],left_verpos_stimevoked{i},'un',0));
    right_verpos_stimevoked{i} = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],right_verpos_stimevoked{i},'un',0));

end

% baseline saccade count
nsaccades_baseline = sum(t.saccade > (t.stimbeg - stimlog.trainduration*1e-3)' & t.saccade < t.stimbeg', 'all');
% baseline saccade times
saccade_time_baseline = t.saccade(logical(sum(t.saccade > (t.stimbeg - stimlog.trainduration*1e-3)' & t.saccade < t.stimbeg',2)));
% baseline saccade direction aligned to saccade onset
for j=1:nsaccades_baseline
    post_saccade_period = (ts.yle > saccade_time_baseline(j) & ts.yle < saccade_time_baseline(j)+prs.post_saccade);
    horpos_baseline{j} = 0.5*(ch.yle(post_saccade_period) + ch.yre(post_saccade_period)); horpos_baseline{j} = horpos_baseline{j} - horpos_baseline{j}(1);
    verpos_baseline{j} = 0.5*(ch.zle(post_saccade_period) + ch.zre(post_saccade_period)); verpos_baseline{j} = verpos_baseline{j} - verpos_baseline{j}(1);
    saccade_dir_baseline(:,j) = [horpos_baseline{j}(nt_saccade) - horpos_baseline{j}(1) ;  verpos_baseline{j}(nt_saccade) - verpos_baseline{j}(1)];
end
minlgth = min(cellfun(@(x) numel(x),horpos_baseline));
horpos_baseline = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],horpos_baseline,'un',0));
verpos_baseline = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],verpos_baseline,'un',0));

% baseline saccade direction aligned to 200ms before stimulation
stim_onset_all = cell2mat(stim_onset');
for j=1:stimlog.nreps*namplitudes
%     pre_stim_period = (ts.yle > (stim_onset_all(j)-stimlog.trainduration*1e-3) & ts.yle < stim_onset_all(j)); % doesn't work for extremely short trains as saccade is not completed
    pre_stim_period = (ts.yle > (stim_onset_all(j)-prs.post_saccade) & ts.yle < stim_onset_all(j));
    left_horpos_baseline_prestim{j} = (ch.yle(pre_stim_period)); left_horpos_baseline_prestim{j} = left_horpos_baseline_prestim{j} - left_horpos_baseline_prestim{j}(1);
    right_horpos_baseline_prestim{j} = (ch.yre(pre_stim_period)); right_horpos_baseline_prestim{j} = right_horpos_baseline_prestim{j} - right_horpos_baseline_prestim{j}(1);
    left_verpos_baseline_prestim{j} = (ch.zle(pre_stim_period)); left_verpos_baseline_prestim{j} = left_verpos_baseline_prestim{j} - left_verpos_baseline_prestim{j}(1);
    right_verpos_baseline_prestim{j} = (ch.zre(pre_stim_period)); right_verpos_baseline_prestim{j} = right_verpos_baseline_prestim{j} - right_verpos_baseline_prestim{j}(1);
    hordir = 0.5*((left_horpos_baseline_prestim{j}(nt_saccade) - left_horpos_baseline_prestim{j}(1)) + (right_horpos_baseline_prestim{j}(nt_saccade) - right_horpos_baseline_prestim{j}(1)));
    verdir = 0.5*((left_verpos_baseline_prestim{j}(nt_saccade) - left_verpos_baseline_prestim{j}(1)) + (right_verpos_baseline_prestim{j}(nt_saccade) - right_verpos_baseline_prestim{j}(1)));
    saccade_dir_baseline_prestim(:,j) = [hordir ;  verdir];
end
minlgth = min(cellfun(@(x) numel(x),left_horpos_baseline_prestim));
left_horpos_baseline_prestim = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],left_horpos_baseline_prestim,'un',0));
right_horpos_baseline_prestim = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],right_horpos_baseline_prestim,'un',0));
left_verpos_baseline_prestim = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],left_verpos_baseline_prestim,'un',0));
right_verpos_baseline_prestim = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],right_verpos_baseline_prestim,'un',0));

%% compute saccade probability
psaccade.mu = nsaccades/stimlog.nreps; 
psaccade.sem = sqrt(((psaccade.mu).*(1 - (psaccade.mu)))/stimlog.nreps);
psaccade_baseline.mu = nsaccades_baseline/(namplitudes*stimlog.nreps); 
psaccade_baseline.sem = sqrt(((psaccade_baseline.mu).*(1 - (psaccade_baseline.mu)))/(namplitudes*stimlog.nreps));

%% plot psychometric function
figure; hold on;
shadedErrorBar(pulseamplitudes, repmat(psaccade_baseline.mu,[1 namplitudes]), repmat(psaccade_baseline.sem,[1 namplitudes]),'lineprops','--k')
h = errorbar(pulseamplitudes,psaccade.mu,psaccade.sem,'or','MarkerFaceColor','r'); h.CapSize = 0;
axis([0 max(pulseamplitudes) + min(pulseamplitudes) 0 1]);
xlabel('Current amplitude (\muA)'); ylabel('Number of saccades per stimulation');
saveas(gcf,['Fig_' num2str(i), '.epsc']);

%% plot direction of saccade aligned to saccade onset
tt = dt.yle:dt.yle:size(horpos_baseline,1)*dt.yle;
figure; hold on;
subplot(namplitudes+1,3,1); hold on; plot(tt,horpos_baseline,'k'); ylim([-40 40]);
subplot(namplitudes+1,3,2); hold on; plot(tt,verpos_baseline,'k'); ylim([-40 40]);
subplot(namplitudes+1,3,3); hold on; quiver(zeros(1,nsaccades_baseline),zeros(1,nsaccades_baseline),saccade_dir_baseline(1,:),saccade_dir_baseline(2,:),'k'); axis([-40 40 -40 40]);

for i=1:namplitudes
    tt = dt.yle:dt.yle:size(horpos{i},1)*dt.yle;
    subplot(namplitudes+1,3,3*i+1); hold on; plot(tt,horpos{i},'r'); ylim([-40 40]);
    subplot(namplitudes+1,3,3*i+2); hold on; plot(tt,verpos{i},'r'); ylim([-40 40]);
    subplot(namplitudes+1,3,3*i+3); hold on; quiver(zeros(1,nsaccades(i)),zeros(1,nsaccades(i)),saccade_dir{i}(1,:),saccade_dir{i}(2,:),'r'); axis([-40 40 -40 40]);
end

subplot(namplitudes+1,3,1); hold on; xlabel('Time from saccade onset (s)'); ylabel('Horizontal eye position (deg)'); 
subplot(namplitudes+1,3,2); hold on; xlabel('Time from saccade onset (s)'); ylabel('Vertical eye position (deg)');
subplot(namplitudes+1,3,3); hold on; title('Saccade amplitude and direction'); 
xlabel('Horizontal eye displacement (deg)'); ylabel('Vertical eye displacement (deg)');

%% plot direction of saccade aligned to stimulation onset
figure(21); hold on; figure(22); hold on;

tt = dt.yle:dt.yle:size(left_horpos_baseline_prestim,1)*dt.yle;
figure(21);
subplot(namplitudes+1,4,1); hold on; plot(tt,left_horpos_baseline_prestim,'k'); ylim([-40 40]);
subplot(namplitudes+1,4,2); hold on; plot(tt,right_horpos_baseline_prestim,'k'); ylim([-40 40]);
subplot(namplitudes+1,4,3); hold on; plot(tt,left_verpos_baseline_prestim,'k'); ylim([-40 40]);
subplot(namplitudes+1,4,4); hold on; plot(tt,right_verpos_baseline_prestim,'k'); ylim([-40 40]);
figure(22);
subplot(namplitudes+1,1,1); hold on; quiver(zeros(1,stimlog.nreps*namplitudes),zeros(1,stimlog.nreps*namplitudes),saccade_dir_baseline_prestim(1,:),saccade_dir_baseline_prestim(2,:),'k'); axis([-40 40 -40 40]);

for i=1:namplitudes
    tt = dt.yle:dt.yle:size(left_horpos_stimevoked{i},1)*dt.yle;
    figure(21);
    subplot(namplitudes+1,4,4*i+1); hold on; plot(tt,left_horpos_stimevoked{i},'r'); ylim([-40 40]);
    subplot(namplitudes+1,4,4*i+2); hold on; plot(tt,right_horpos_stimevoked{i},'r'); ylim([-40 40]);
    subplot(namplitudes+1,4,4*i+3); hold on; plot(tt,left_verpos_stimevoked{i},'r'); ylim([-40 40]);
    subplot(namplitudes+1,4,4*i+4); hold on; plot(tt,right_verpos_stimevoked{i},'r'); ylim([-40 40]);
    figure(22);
    subplot(namplitudes+1,1,1*i+1); hold on; quiver(zeros(1,stimlog.nreps),zeros(1,stimlog.nreps),saccade_dir_stimevoked{i}(1,:),saccade_dir_stimevoked{i}(2,:),'r'); axis([-40 40 -40 40]);
    
end

figure(21); hold on;
subplot(namplitudes+1,4,1); hold on; xlabel('Time (s)'); ylabel('Left Horizontal (deg)');
subplot(namplitudes+1,4,2); hold on; xlabel('Time (s)'); ylabel('Right Horizontal (deg)');
subplot(namplitudes+1,4,3); hold on; xlabel('Time (s)'); ylabel('Left Vertical (deg)');
subplot(namplitudes+1,4,4); hold on; xlabel('Time (s)'); ylabel('Right Vertical (deg)');
saveas(gcf,['Fig_21' num2str(i), '.epsc']);

figure(22); hold on;
subplot(namplitudes+1,1,1); hold on; title('Saccade amplitude and direction'); 
xlabel('Horizontal eye displacement (deg)'); ylabel('Vertical eye displacement (deg)');
saveas(gcf,['Fig_22' num2str(i), '.epsc']);

%% Compute saccade endpoints
pulseamplitudes = unique(stimlog.pulseamplitude); namplitudes = length(pulseamplitudes);
nsaccades = zeros(1,namplitudes);
nt_saccade = round(prs.saccade_duration/dt.yle);
for i=1:namplitudes
    % saccade count
    nsaccades(i) = sum(t.saccade > t.stimbeg(stimlog.pulseamplitude == pulseamplitudes(i))' & ...
        t.saccade < t.stimend(stimlog.pulseamplitude == pulseamplitudes(i))', 'all');
    % saccade times
    saccade_time{i} = t.saccade(logical(sum(t.saccade > t.stimbeg(stimlog.pulseamplitude == pulseamplitudes(i))' & ...
        t.saccade < t.stimend(stimlog.pulseamplitude == pulseamplitudes(i))',2)));
    % stimulation times
    stim_onset{i} = t.stimbeg(stimlog.pulseamplitude == pulseamplitudes(i));
    % saccade endpoint aligned to saccade onset
    for j=1:nsaccades(i)
        post_saccade_period = (ts.yle > saccade_time{i}(j) & ts.yle < saccade_time{i}(j)+prs.post_saccade);
        horpos_temp{i}{j} = 0.5*(ch.yle(post_saccade_period) + ch.yre(post_saccade_period)); 
        verpos_temp{i}{j} = 0.5*(ch.zle(post_saccade_period) + ch.zre(post_saccade_period)); 
        saccade_endpoint{i}(:,j) = [horpos_temp{i}{j}(nt_saccade) ;  verpos_temp{i}{j}(nt_saccade)];
    end
    minlgth = min(cellfun(@(x) numel(x),horpos_temp{i}));
    horpos{i} = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],horpos_temp{i},'un',0));
    verpos{i} = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],verpos_temp{i},'un',0));

    % saccade endpoint aligned to stimulation onset
    for j=1:stimlog.nreps
        post_stim_period = (ts.yle > stim_onset{i}(j) & ts.yle < stim_onset{i}(j)+stimlog.trainduration*1e-3); % doesn't work for extremely short trains as saccade is not completed
        post_stim_period = (ts.yle > stim_onset{i}(j) & ts.yle < stim_onset{i}(j)+prs.post_saccade);
        left_horpos_stimevoked_temp{i}{j} = (ch.yle(post_stim_period)); 
        right_horpos_stimevoked_temp{i}{j} = (ch.yre(post_stim_period));
        left_verpos_stimevoked_temp{i}{j} = (ch.zle(post_stim_period)); 
        right_verpos_stimevoked_temp{i}{j} = (ch.zre(post_stim_period));
        horendpoint = 0.5*(left_horpos_stimevoked_temp{i}{j}(nt_saccade) + right_horpos_stimevoked_temp{i}{j}(nt_saccade));
        verendpoint = 0.5*(left_verpos_stimevoked_temp{i}{j}(nt_saccade) + right_verpos_stimevoked_temp{i}{j}(nt_saccade));
        saccade_endpoint_stimevoked{i}(:,j) = [horendpoint ; verendpoint];
    end
    minlgth = min(cellfun(@(x) numel(x),left_horpos_stimevoked_temp{i}));
    left_horpos_stimevoked{i} = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],left_horpos_stimevoked_temp{i},'un',0));
    right_horpos_stimevoked{i} = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],right_horpos_stimevoked_temp{i},'un',0));
    left_verpos_stimevoked{i} = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],left_verpos_stimevoked_temp{i},'un',0));
    right_verpos_stimevoked{i} = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],right_verpos_stimevoked_temp{i},'un',0));
end

% baseline saccade count
nsaccades_baseline = sum(t.saccade > (t.stimbeg - stimlog.trainduration*1e-3)' & t.saccade < t.stimbeg', 'all');
% baseline saccade times
saccade_time_baseline = t.saccade(logical(sum(t.saccade > (t.stimbeg - stimlog.trainduration*1e-3)' & t.saccade < t.stimbeg',2)));
% baseline saccade endpoints aligned to saccade onset
for j=1:nsaccades_baseline
    post_saccade_period = (ts.yle > saccade_time_baseline(j) & ts.yle < saccade_time_baseline(j)+prs.post_saccade);
    horpos_baseline_temp{j} = 0.5*(ch.yle(post_saccade_period) + ch.yre(post_saccade_period));
    verpos_baseline_temp{j} = 0.5*(ch.zle(post_saccade_period) + ch.zre(post_saccade_period));
    saccade_endpoint_baseline(:,j) = [horpos_baseline_temp{j}(nt_saccade) ;  verpos_baseline_temp{j}(nt_saccade)];
end
minlgth = min(cellfun(@(x) numel(x),horpos_baseline_temp));
horpos_baseline = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],horpos_baseline_temp,'un',0));
verpos_baseline = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],verpos_baseline_temp,'un',0));

% baseline saccade endpoints aligned to 200ms before stimulation
stim_onset_all = cell2mat(stim_onset');
for j=1:stimlog.nreps*namplitudes
%     pre_stim_period = (ts.yle > (stim_onset_all(j)-stimlog.trainduration*1e-3) & ts.yle < stim_onset_all(j)); % doesn't work for extremely short trains as saccade is not completed
    pre_stim_period = (ts.yle > (stim_onset_all(j)-prs.post_saccade) & ts.yle < stim_onset_all(j));
    left_horpos_baseline_prestim_temp{j} = (ch.yle(pre_stim_period)); 
    right_horpos_baseline_prestim_temp{j} = (ch.yre(pre_stim_period));
    left_verpos_baseline_prestim_temp{j} = (ch.zle(pre_stim_period)); 
    right_verpos_baseline_prestim_temp{j} = (ch.zre(pre_stim_period));
    horendpoint = 0.5*(left_horpos_baseline_prestim_temp{j}(nt_saccade) + right_horpos_baseline_prestim_temp{j}(nt_saccade));
    verendpoint = 0.5*(left_verpos_baseline_prestim_temp{j}(nt_saccade) + right_verpos_baseline_prestim_temp{j}(nt_saccade));
    saccade_endpoint_baseline_prestim(:,j) = [horendpoint ;  verendpoint];
end
minlgth = min(cellfun(@(x) numel(x),left_horpos_baseline_prestim_temp));
left_horpos_baseline_prestim = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],left_horpos_baseline_prestim_temp,'un',0));
right_horpos_baseline_prestim = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],right_horpos_baseline_prestim_temp,'un',0));
left_verpos_baseline_prestim = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],left_verpos_baseline_prestim_temp,'un',0));
right_verpos_baseline_prestim = cell2mat(cellfun(@(x) [x(1:minlgth) ; nan(minlgth-length(x),1)],right_verpos_baseline_prestim_temp,'un',0));

%% plot distribution of saccade endpoints
N = 50;
nx = linspace(-40,40,N+1);

figure(32);
subplot(namplitudes+1,2,1); hold on; histogram(saccade_endpoint_baseline_prestim(1,:),nx,'facecolor','k'); xlim([-40 40]);
subplot(namplitudes+1,2,2); hold on; histogram(saccade_endpoint_baseline_prestim(2,:),nx,'facecolor','k'); xlim([-40 40]);


for i=1:namplitudes
    figure(32);
    subplot(namplitudes+1,2,2*i+1); hold on; histogram(saccade_endpoint_stimevoked{i}(1,:),nx); xlim([-40 40]);
    subplot(namplitudes+1,2,2*i+2); hold on; histogram(saccade_endpoint_stimevoked{i}(2,:),nx); xlim([-40 40]);
end

figure(32); hold on;
subplot(namplitudes+1,2,1); hold on; title('Saccade endpoint (hor)'); xlabel('Horizontal eye displacement (deg)'); 
subplot(namplitudes+1,2,2); hold on; title('Saccade endpoint (ver)'); xlabel('Vertical eye displacement (deg)');
saveas(gcf,['Fig_22' num2str(i), '.epsc']);

% Test significance of difference between baseline and stimulation
[h,p] = arrayfun(@(i) ttest2(saccade_endpoint_baseline_prestim',saccade_endpoint_stimevoked{i}'),1:namplitudes,'un',0);
h = cell2mat(h');
p = cell2mat(p');

