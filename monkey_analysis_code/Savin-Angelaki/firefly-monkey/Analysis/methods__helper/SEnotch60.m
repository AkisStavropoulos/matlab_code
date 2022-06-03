
% Power line noise filter (60Hz)
function lfp_f = SEnotch60(lfp,fs)

d = designfilt('bandstopiir','FilterOrder',2, ...
'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
'DesignMethod','butter','SampleRate',fs);
lfp_f = filtfilt(d,lfp);