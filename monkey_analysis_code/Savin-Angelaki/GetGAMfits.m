function pgam = GetGAMfits(fits_path,name_pattern)

% get the python GAM fit results

cd(fits_path);
fitfile = dir(fullfile(fits_path,name_pattern));
sessnum = arrayfun(@(x) strsplit(x.name,{'s','_'}),fitfile,'un',0); sessnum = cellfun(@(x) str2double(x{2}),sessnum);
[~,I] = sort(sessnum);
fitfile = fitfile(I);
Nsessions = numel(fitfile);

for n = 1:Nsessions
    temp = importdata(fitfile(n).name);
    units = unique([temp.neuron]);
    Nunits = numel(units);
    Nvar = numel(temp)/Nunits;

    for k = 1:Nunits
        pgam(n).units(k).info = [];
        indx = arrayfun(@(x) x.neuron == units(k), temp);
        
        pgam(n).units(k).info = [pgam(n).units(k).info temp(indx)];
    end
end