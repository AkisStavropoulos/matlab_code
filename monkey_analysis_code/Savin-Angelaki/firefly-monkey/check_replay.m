%% set path
% if contains(computer,'MAC')
%     separ = '/';
% else
%     separ = '\';
% end
% this_path = pwd;
% path_split = split(this_path,separ);
% base_fold = join(path_split(1:end-1),separ);
% if ~contains(path,base_fold{1})
%     addpath(fullfile(base_fold{1},'genpath2'))
%     addpath(genpath2(base_fold{1},{'.git','genpath2'}))
%     
% end
%% load and plot
%smr_active = ImportSMR('Z:\Data\Monkey2_newzdrive\Schro\Utah Array\Mar 02 2018\behavioural data\m53c0325.smr');
%smr_passive = ImportSMR('Z:\Data\Monkey2_newzdrive\Schro\Utah Array\Mar 02 2018\behavioural data\m53c0326.smr');
% smr_active = ImportSMR('/Volumes/server/Data/Monkey2_newzdrive/Schro/Utah Array/Mar 02 2018/behavioural data/m53c0325.smr');
% ?smr_passive = ImportSMR('/Volumes/server/Data/Monkey2_newzdrive/Schro/Utah Array/Mar 02 2018/behavioural data/m53c0326.smr');
%% ciao
data = smr_active;
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
chno.xmp = find(strcmp(ch_title,'MonkeyX'));
chno.v = find(strcmp(ch_title,'ForwardV')); chno.w = find(strcmp(ch_title,'AngularV'));

chnames = fieldnames(chno); MAX_LENGTH = inf; dt = [];
scaling.xmp = data(chno.xmp).hdr.adc.Scale;
scaling.v = data(chno.v).hdr.adc.Scale;
offset.v = data(chno.v).hdr.adc.DC;
r_vel_active = double(data(chno.v).imp.adc)*scaling.v + offset.v;
scaling.w = data(chno.w).hdr.adc.Scale; offset.w = data(chno.w).hdr.adc.DC;
theta_vel_active = double(data(chno.w).imp.adc)*scaling.w + offset.w;
offset.xmp = data(chno.xmp).hdr.adc.DC;
xmp_active = double(data(chno.xmp).imp.adc)*scaling.xmp + offset.xmp;
chno.mrk = find(strcmp(ch_title,'marker'));
scaling.t = data(chno.mrk).hdr.tim.Scale*data(chno.mrk).hdr.tim.Units;
markers = data(chno.mrk).imp.mrk(:,1);
t.events = double(data(chno.mrk).imp.tim)*scaling.t;
t.beg_active = t.events(markers ==2);
t.end_active = t.events(markers ==3);
dt = prod(data(chno.v).hdr.adc.SampleInterval);

data = smr_passive;
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
chno.xmp = find(strcmp(ch_title,'MonkeyX'));
scaling.xmp = data(chno.xmp).hdr.adc.Scale;
offset.xmp = data(chno.xmp).hdr.adc.DC;
xmp_passive = double(data(chno.xmp).imp.adc)*scaling.xmp + offset.xmp;
chno.mrk = find(strcmp(ch_title,'marker'));
scaling.t = data(chno.mrk).hdr.tim.Scale*data(chno.mrk).hdr.tim.Units;
scaling.v = data(chno.v).hdr.adc.Scale;
offset.v = data(chno.v).hdr.adc.DC;
r_vel_passive = double(data(chno.v).imp.adc)*scaling.v + offset.v;
scaling.w = data(chno.w).hdr.adc.Scale; offset.w = data(chno.w).hdr.adc.DC;
theta_vel_passive = double(data(chno.w).imp.adc)*scaling.w + offset.w;
markers = data(chno.mrk).imp.mrk(:,1);
t.events = double(data(chno.mrk).imp.tim)*scaling.t;
t.beg_passive = t.events(markers ==2);
t.end_passive = t.events(markers ==3);
if t.beg_active(end)> t.end_active(end)
    t.beg_active = t.beg_active(1:end-1);
end
if t.beg_passive(end)> t.end_passive(end)
    t.beg_passive = t.beg_passive(1:end-1);
    
end
diff_active =  t.end_active -  t.beg_active;
diff_passive = t.end_passive - t.beg_passive;
% plot(diff_active(1:100))
% hold on
% plot(diff_passive(1:100))
% title('trial duration')

%% get the movement start
prs = default_prs(53,47);
v_thresh = prs.v_thresh;
w_thresh = prs.w_thresh;
v_time2thresh = prs.v_time2thresh;
indx = find(r_vel_active > v_thresh, 1);
orign_indx = 1:length(r_vel_active);
idx_start_active = indx;
idx_stop_active = [];
vel=r_vel_active;
avel = theta_vel_active;
ii0 = 0;
while ~isempty(indx)
    % do something
    indx_end = indx + find(abs(vel(indx:end)) < v_thresh & abs(avel(indx:end)) < w_thresh,1); % first downward threshold-crossing
    if isempty(indx_end)
        indx_end = length(vel);
        
    end
    idx_stop_active = [idx_stop_active; ii0+indx_end];
    vel = vel(indx_end+1:end);
    avel = avel(indx_end+1:end);
    ii0 = ii0 + indx_end;
    indx = find(vel > v_thresh, 1);
    idx_start_active = [idx_start_active; ii0+indx];

end

%% get the passive

indx = find(r_vel_passive > v_thresh, 1);
orign_indx = 1:length(r_vel_passive);
idx_start_passive = indx;
idx_stop_passive = [];
vel=r_vel_passive;
avel = theta_vel_passive;
ii0 = 0;
while ~isempty(indx)
    % do something
    indx_end = indx + find(abs(vel(indx:end)) < v_thresh & abs(avel(indx:end)) < w_thresh,1); % first downward threshold-crossing
    if isempty(indx_end)
        indx_end = length(vel);
        
    end
    idx_stop_passive = [idx_stop_passive; ii0+indx_end];
    vel = vel(indx_end+1:end);
    avel = avel(indx_end+1:end);
    ii0 = ii0 + indx_end;
    indx = find(vel > v_thresh, 1);
    idx_start_passive = [idx_start_passive; ii0+indx];

end
 