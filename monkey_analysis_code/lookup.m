%% create a file called lookup.m with the following:
sessions(1).foldername = 'C:\Users\stavropo\Documents\Firefly analysis\Data\Behavioural training\Jan 17 2018'
sessions(2).foldername = 'C:\...'
...
%% your analysis code
lookup;
for i=1:length(sessions)
    cd(sessions(i).foldername)
    dir('..')
    ...
    your_analysis_code_here
    ...
end
