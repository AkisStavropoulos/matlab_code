% fpath = 'Z:\Users\Akis\Ripple Test\Moog 3\May 26 2021';
% fname = 'moog3_saline_test_3.ns5';
% fpath = 'Z:\Data\Monkey2_newzdrive\Marco\U-probe\Mar 29 2021\neural data';
% fname = 'marco_03_29_2021.plx';
% fpath = 'Z:\Data\Monkey2_newzdrive\Viktor\U-probe\Mar 30 2021\neural data';
% fname = 'm71s387.ns5';
% fpath = 'Z:\Data\Monkey2_newzdrive\Jimmy\U-probe\Jul 29 2021\neural data';
% fname = 'm73s14.ns5';
fpath = 'Z:\Users\Akis\Ripple Test\Moog 3\Sep 30 2021_Weird Noise';
fname = 'wied_noise_1.ns5';

cd(fpath);
fullname = fullfile(fpath,fname);
[fpath,fname,ftype] = fileparts(fullname);


MakeDatFile(fullname);

nch = 32;
nsamples = 60*30000; % plexon: 20000, ripple: 30000
precision = '*int16';
startsample = 10;
[data,startsample] = ReadDatFile([fname '.dat'],nch,nsamples,precision,startsample);
data = double(data);

t = 1:size(data,2);
for n = 1:nch
    nanx = isnan(data(n,:));
    data(n,nanx) = interp1(t(~nanx),data(n,~nanx),t(nanx));
    
end

Fs = 30000;            % Sampling frequency                    
T = 1/Fs;              % Sampling period       
L = size(data,2);      % Length of signal
t = (0:L-1)*T;         % Time vector

figure; hold on;
for n = 1:nch
Y = fft(data(n,:));

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

plot(f,P1)
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([0 120])
title(fname);
end
