%% Ripple markers testing
% Directory changed back to array data
cd('D:\ripple marker test moog1\test 17') % cd('Z:\Data\Monkey2_newzdrive\Viktor\Utah Array\May 06 2021\neural data');
file_nev=dir('*.nev'); prs.neur_filetype = 'nev';

fprintf(['... reading events from ' file_nev.name '\n']);
[events_nev,prs] = GetEvents_nev(file_nev.name,prs); % requires package from Blackrock Microsystems: https://github.com/BlackrockMicrosystems/NPMK

file_smr=dir('*.smr');
prs = default_prs(00,1);
data_smr = ImportSMR(file_smr.name);
trials_smr = AddSMRData(data_smr,prs);

events_smr.t_beg = arrayfun(@(x) x.events.t_beg, trials_smr);
events_smr.t_end = arrayfun(@(x) x.events.t_end + x.events.t_beg, trials_smr);
events_smr.t_rew = arrayfun(@(x) x.events.t_rew + x.events.t_beg, trials_smr);
events_smr.t_beg_correction = arrayfun(@(x) x.events.t_beg_correction, trials_smr);

% compare markers

figure; hold on;
subplot(3,1,1); hold on;
hist(diff(events_smr.t_beg-events_smr.t_beg_correction - events_nev.t_beg),100); title('t_{beg}: Spike2 - Ripple');
subplot(3,1,2); hold on;
hist(diff(events_smr.t_end - events_nev.t_end),100); title('t_{end}: Spike2 - Ripple');
subplot(3,1,3); hold on;
hist(diff(events_smr.t_rew(~isnan(events_smr.t_rew)) - events_nev.t_rew),100); title('t_{reward}: Spike2 - Ripple');

%%
t_end_nev = [events_nev.t_end];
t_beg_nev = [events_nev.t_beg];
if length(t_beg_nev) > length(t_end_nev)
    if all(t_end_nev - t_beg_nev(1:length(t_end_nev)) > 0)
        events_nev.t_beg = t_beg_nev(1:length(t_end_nev));
    elseif all(t_end_nev - t_beg_nev(length(t_beg_nev)-length(t_end_nev)+1:end) > 0)
        events_nev.t_beg = t_beg_nev(length(t_beg_nev)-length(t_end_nev)+1:end);
    end
end
