function [CorrTbl,Correlation] = CorrTable(rho,pval,legend_input)
%% Create Table with response bias information
plusminus = sprintf('\x00b1'); % unicode characters
middledot = sprintf('\x00b7'); % unicode characters

mods = unique(cellfun(@(x) strtok(x,' '),legend_input,'uniformoutput',false),'stable');
subgrps = cellfun(@(x) ['Subject ' num2str(x)], num2cell(1:size(rho.r,1)),'uniformoutput',false);

for i = 1:size(rho.r,1)
    for s = 1:size(rho.r,2)
        if pval.r(i,s) < 10^-3; p_r = num2str(pval.r(i,s),'%.1e'); p_r = [strtok(p_r,'e') middledot '10^_{' p_r(strfind(p_r,'-'):end) '}'];
        else; p_r = num2str(pval.r(i,s),'%4.3f'); end
        
        if pval.th(i,s) < 10^-3; p_th = num2str(pval.th(i,s),'%.1e'); p_th = [strtok(p_th,'e') middledot '10^_{' p_th(strfind(p_th,'-'):end) '}'];
        else; p_th = num2str(pval.th(i,s),'%4.3f'); end
        
    Correlation.r{i,s} = ['r = ' num2str(rho.r(i,s),'%4.3f') ', p = ' p_r];
    Correlation.th{i,s} = ['r = ' num2str(rho.th(i,s),'%4.3f') ', p = ' p_th];
    CorrTbl.r(i,s) = table(Correlation.r(i,s));
    CorrTbl.th(i,s) = table(Correlation.th(i,s));

    end
end

CorrTbl.r.Properties.RowNames = subgrps;    CorrTbl.r.Properties.VariableNames = mods;
CorrTbl.th.Properties.RowNames = subgrps;    CorrTbl.th.Properties.VariableNames = mods;
