%function updateprog(mice); %Input ['2' '3' '4']

mice={'02' '03' '04'};

for ii=1:length(mice)

Lname1='BCIT';
Lname2='.mat';
Lmouse=strcat(Lname1,mice{ii},Lname2); %create name of file to load
load(Lmouse);

%% choose dir
Dname1='/Volumes/server/Data/Mice/JeanPaul/VR/Subjects/TRAINING_PHASE_2/BCIT';
Dname2='/smr';
Dmouse=strcat(Dname1,mice{ii},Dname2);

d=dir([Dmouse,'/*.smr']);

%% Find Names of all files in direcory
NewNames={};
for w=1:size(d,1)
    NewNames{w}=d(w).name(1:17);
end

%% Find names of already saved files
OldNames={};
for f=1:size(allSMR,2)
    OldNames{f}=allSMR(f).name(1:17);
end

%% Find The Paths of the files to be imported
new=setdiff(NewNames,OldNames); %find different names

Paths=[]; Names=[];
for j=1:size(new,2)                                 
    idx=strcmp(new(j),NewNames);                    
    
    %the case that we have a single .smr file
    if sum(idx)==1
       pth=strcat(d(idx).folder,'/',d(idx).name);
       Paths=[Paths {pth}]; Names=[Names {d(idx).name}];
    
    %the case we have two .smr files
    else
        
        pos=find(idx==1);
        for g=1:sum(idx)
            pth=strcat(d(pos(g)).folder,'/',d(pos(g)).name);
            Paths=[Paths {pth}]; Names=[Namesd {d(pos(g)).name}];
        end  
        
    end
    
end

%% Import Files

for x=1:length(Paths)
    data=ImportSMR(Paths{x});
     

%check channel headers
nch = length(data);
ch_title = cell(1,nch);
hdr = {data.hdr};
for i=1:nch
    if ~isempty(hdr{i})
        ch_title{i} = hdr{i}.title;
    else
        ch_title{i} = 'nan';
    end
end

chno.mrk = find(strcmp(ch_title,'marker'));
markers = data(chno.mrk).imp.mrk(:,1);

sz=length(allSMR); %calculate size to find last row with it
allSMR(sz+1).name=Names{x}(1:17); %get names
allSMR(sz+1).trls=sum(markers==4); %get trials
allSMR(sz+1).rwd=sum(markers==6); %get rwd
    
end
%% Save mouse
Sname1='BCIT';
Smouse=strcat(Sname1,mice{ii});
save(Smouse,'allSMR')

%% Display
show=strcat('BCIT',mice{ii},' **Upload Completed**');
disp(show)

end