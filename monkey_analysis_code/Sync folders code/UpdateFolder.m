function [FolderNames] = UpdateFolder(PathName,UpdPath,UpdSubFolders)
% PathName: Up-to-date folder
% UpdPath: Folder to be updated
% UpdSubF: 1 if you want to update subfolders too, 0 if not
%
% Outputs the subfolder so you can run the same procedure for the
% subfolders

% PathName = 'G:\My Drive\SF w MC\Paper'; 
% UpdPath = 'Z:\Users\Akis\Paper';

% Make sure a folder exists on destination path
if ~isfolder(UpdPath)
    mkdir(UpdPath)
    disp(['Folder "' UpdPath '" created.'])
end
%% Update Folder content

% Up-to-date folder
Fcontent = dir(PathName);
Fcontent = rmdotfiles(Fcontent);
SubFolders = [];
for i = 1:length(Fcontent)
    if Fcontent(i).isdir == 1
        SubFolders = [SubFolders i];
    end
end
FolderNames = {Fcontent(SubFolders).name};
FileIndx = setdiff(1:length(Fcontent),SubFolders);
FileNames = {Fcontent(FileIndx).name};


% Folder to be updated
UFcontent = dir(UpdPath);
UFcontent = rmdotfiles(UFcontent);
USubFolders = [];
for i = 1:length(UFcontent)
    if UFcontent(i).isdir == 1
        USubFolders = [USubFolders i];
    end
end
UFolderNames = {UFcontent(USubFolders).name};
UFileIndx = setdiff(1:length(UFcontent),USubFolders);
UFileNames = {UFcontent(UFileIndx).name};

% Find new files
new = setdiff(FileNames,UFileNames);
if ~isempty(new)
    Nblank = 30;
    status = [];    msg = [];   report_temp = [];
    for i = 1:length(new)
        [status(i),msg{i}] = copyfile([PathName '\' new{i}],UpdPath) ;
        
        ind = length(new{i});
        if ind >= Nblank
            ind = Nblank;
        end
        report_temp{i} = [num2str(i) ') ' new{i}(1:ind) blanks(Nblank-ind) ...
            ': status ' num2str(status(i)) ', message: ' msg{i} '\n'];
    end
    fprintf(['\n New files: \n' [report_temp{:}] '\n']);
end

% Update all
status = [];    msg = [];   report_temp = [];
for i = 1:length(FileNames)
   [status(i),msg{i}] = copyfile([PathName '\' FileNames{i}],UpdPath) ;
   if status(i) == 0
      fprintf(['Error: ' FileNames{i} '\n' msg{i} '\n']) 
   end
end

%% Update Subfolders
if UpdSubFolders
status = [];    msg = [];  
for i = 1:length(FolderNames)
   newPath = [PathName,'\',FolderNames{i}];
   newDest = [UpdPath,'\',FolderNames{i}];

   [status(i),msg{i}] = copyfile(newPath,newDest) ;
   
   if status(i) == 0
      fprintf(['Error: ' FileNames{i} '\n' msg{i} '\n']) 
   end
end
end

