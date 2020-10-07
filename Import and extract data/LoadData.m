function [subject,DataTable] = LoadData(data_folder,subject_name)
%% load data

% to load certain people from subject_name list, choose the indices of
% subject_name before you execute the function, e.g. subject name([1 3 7])
disp('...... Loading Data');

clear subject
for i = 1:length(subject_name)
    cd([data_folder subject_name{i}]);
    subject(i) = load(subject_name{i});
end
for i = 1:length(subject)
    subject(i).name = subject(i).trials(1).prs.subject;
end

% CreateDataTable(subject);
cd(data_folder);
if 0
DataTable = load('DataTable');
else
    DataTable = [];
end
% DataTable = DataTable.data;

disp('...Data Loaded')
