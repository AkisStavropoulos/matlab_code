%% Ripple markers testing for Unity PROCESS NOISE
% Directory changed back to array data
neural_path = 'Z:\Users\Akis\Process noise tests\Monkey setup\neural data'; % 'Z:\Users\Akis\Ripple Test\Moog 1\Nov 15 2021\neural data\PM'; % neural_path =  'Z:\Data\Monkey2_newzdrive\Jimmy\U-probe\Sep 16 2021\neural data'; 
% neural_path = 'Y:\Projects\Monkey_Unity_VR\data\testData\20210907\neural data'; % 'Z:\Data\Monkey2_newzdrive\Jimmy\U-probe\Jul 14 2021\neural data'; % 'Z:\Data\Monkey2_newzdrive\Jimmy\U-probe\Jul 02 2021\neural data'; % 'Z:\Users\Akis\Ripple Test\Moog 3\Jul 08 2021\neural data'; % 'D:\ripple marker test moog1\test 17';
cd(neural_path) 
file_nev=dir('*.nev'); prs.neur_filetype = 'nev';

fprintf(['... reading events from ' file_nev.name '\n']);
[events_nev,prs] = GetEvents_nev(file_nev.name,prs); % requires package from Blackrock Microsystems: https://github.com/BlackrockMicrosystems/NPMK

behv_path = 'Z:\Users\Akis\Process noise tests\Monkey setup\behavioural data'; % 'Z:\Users\Akis\Ripple Test\Moog 1\Nov 15 2021\behavioural data\PM2'; % behv_path = 'Z:\Data\Monkey2_newzdrive\Jimmy\U-probe\Sep 16 2021\behavioural data'; 
% behv_path = 'Y:\Projects\Monkey_Unity_VR\data\testData\20210907\behavioural data\001'; % 'Z:\Data\Monkey2_newzdrive\Jimmy\U-probe\Jul 29 2021\behavioural data';   % 'Z:\Users\Akis\Ripple Test\Moog 3\Jul 08 2021\behavioural data';
cd(behv_path)
file_behv = dir('discontinuous*.txt');
if numel(file_behv) > 1
    events_vr.t_end = []; events_vr.t_beg = []; events_vr.t_rew = [];
    for n = 1:numel(file_behv)
        data1 = table2array(readtable(file_behv(n).name));
        events_vr.t_end = [events_vr.t_end data1(:,end-15)'];
        events_vr.t_beg = [events_vr.t_beg data1(:,end-18)'];
        events_vr.t_rew = [events_vr.t_rew data1(:,end-16)'];
    end
else
    data1 = table2array(readtable(file_behv(1).name));
    events_vr.t_end = data1(:,end-15)';
    events_vr.t_beg = data1(:,end-18)';
    events_vr.t_rew = data1(:,end-16)';
end

% compare markers
for n = 1:numel(events_nev.t_start)
    fldnames = fields(events_nev);
    for k = 1:numel(fldnames)
        if n==1; ripple_mrk.(fldnames{k}) = []; end
        if any(strcmp(fldnames{k},{'t_rew','t_beg','t_end'}))
            indx = find(events_nev.(fldnames{k}) > events_nev.t_start(n) & events_nev.(fldnames{k}) < events_nev.t_stop(n));
            ripple_mrk.(fldnames{k}) = [ripple_mrk.(fldnames{k}) events_nev.(fldnames{k})(indx) - events_nev.t_start(n)];
        end
    end
end

N = numel(events_vr.t_end);
figure; hold on;
subplot(3,1,1); hold on;
hist(events_vr.t_beg - ripple_mrk.t_beg,100); title('t_{beg}: Unity - Ripple');
subplot(3,1,2); hold on;
hist(events_vr.t_end - ripple_mrk.t_end,100); title('t_{end}: Unity - Ripple');
subplot(3,1,3); hold on;
hist(events_vr.t_rew(events_vr.t_rew>0) - ripple_mrk.t_rew,100); title('t_{reward}: Unity - Ripple'); xlabel('time [s]');

% check drift
figure;plot(events_vr.t_beg ,events_vr.t_beg - ripple_mrk.t_beg); hold on;
plot(events_vr.t_end, events_vr.t_end - ripple_mrk.t_end)
xlabel('unity time [s]'); ylabel('Unity-Ripple time [s]');

% slope of time drift between ripple and unity
b_dt = regress([events_vr.t_beg - ripple_mrk.t_beg]',events_vr.t_beg');

%% Check markers between behavioura and neural data are the same
events_vr
rewards = sum(events_vr.t_rew > 0)
events_nev

%%
% t_end_nev = [events_nev.t_end];
% t_beg_nev = [events_nev.t_beg];
% if length(t_beg_nev) > length(t_end_nev)
%     if all(t_end_nev - t_beg_nev(1:length(t_end_nev)) > 0)
%         events_nev.t_beg = t_beg_nev(1:length(t_end_nev));
%     elseif all(t_end_nev - t_beg_nev(length(t_beg_nev)-length(t_end_nev)+1:end) > 0)
%         events_nev.t_beg = t_beg_nev(length(t_beg_nev)-length(t_end_nev)+1:end);
%     end
% end
