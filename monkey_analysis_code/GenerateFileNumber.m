function [fname,fileNum] = GenerateFileNumber(FileName)

[fPath, fName, fExt] = fileparts(FileName);
[~,~,~,ind] = sscanf(flip(fName),'%d');
startname = fName(1:end-ind+1);

fld_contents = dir(fPath); 

while any(arrayfun(@(x) contains(x.name,fName), fld_contents)) %exist(FileName, 'file')
  [fPath, fName, fExt] = fileparts(FileName);
  [~,~,~,ind] = sscanf(flip(fName),'%d');
  startname = fName(1:end-ind+1);
  % Get number of files:
  fDir     = dir(fullfile(fPath, [fName, '*', fExt]));
  fStr     = sprintf('%s*', fDir.name); %lower(sprintf('%s*', fDir.name));
  fNum     = sscanf(fStr, [startname, '%d', fExt, '*']);
  newNum   = max(fNum) + 1;
  FileName = fullfile(fPath, [startname, sprintf('%03d', newNum), fExt]);

  [fPath, fName, fExt] = fileparts(FileName);
end

% output file name
[~,fName,fExt] = fileparts(FileName);
fname = [fName,  fExt];
fileNum = strsplit(fName,'_'); fileNum = fileNum{end};

