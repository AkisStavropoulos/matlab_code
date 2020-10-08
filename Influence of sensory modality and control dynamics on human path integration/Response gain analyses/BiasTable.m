function [BiasTbl,Bias] = BiasTable(mu,sem,legend_input)
%% Create Table with response bias information
plusminus = sprintf('\x00b1'); % unicode characters
tau = sprintf('\x03c4');

mods = unique(cellfun(@(x) strtok(x,' '),legend_input,'uniformoutput',false),'stable');
taugrps = unique(cellfun(@(x) x(strfind(x,'['):end),legend_input,'uniformoutput',false),'stable');
taugrps = cellfun(@(x) [tau ': ' x],taugrps,'uniformoutput',false);

for s = 1:size(mu.r,1)
    for t = 1:size(mu.r,2)
    Bias.r{t,s} = [num2str(mu.r(s,t),'%4.3f') ' ' plusminus ' ' num2str(sem.r(s,t),'%4.3f')];
    Bias.th{t,s} = [num2str(mu.th(s,t),'%4.3f') ' ' plusminus ' ' num2str(sem.th(s,t),'%4.3f')];
    BiasTbl.r(t,s) = table(Bias.r(t,s));
    BiasTbl.th(t,s) = table(Bias.th(t,s));

    end
end

BiasTbl.r.Properties.RowNames = taugrps;    BiasTbl.r.Properties.VariableNames = mods;
BiasTbl.th.Properties.RowNames = taugrps;    BiasTbl.th.Properties.VariableNames = mods;