function extractMCvariables(mcvariables)
% extracts variable names from MCvariables struct and evaluates them
varnames = fieldnames(mcvariables);

for i = 1:length(varnames)
    assignin('caller',varnames{i},mcvariables.(varnames{i}));
end
    