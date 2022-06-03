function fname = SaveFileWithNumber(FileName, Data)
[fPath, fName, fExt] = fileparts(FileName);
if isempty(fExt)  % No '.mat' in FileName
  fExt     = '.mat';
  FileName = fullfile(fPath, [fName, fExt]);
end
[~,~,~,ind] = sscanf(flip(fName),'%d');
startname = fName(1:end-ind+1);

while exist(FileName, 'file')
  [fPath, fName, fExt] = fileparts(FileName);
  [~,~,~,ind] = sscanf(flip(fName),'%d');
  startname = fName(1:end-ind+1);
  % Get number of files:
  fDir     = dir(fullfile(fPath, [fName, '*', fExt]));
  fStr     = sprintf('%s*', fDir.name); %lower(sprintf('%s*', fDir.name));
  fNum     = sscanf(fStr, [startname, '%d', fExt, '*']);
  newNum   = max(fNum) + 1;
  FileName = fullfile(fPath, [startname, sprintf('%03d', newNum), fExt]);
end
save(FileName, 'Data');

% output file name
[~,fName,fExt] = fileparts(FileName);
fname = [fName,  fExt];

