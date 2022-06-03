function MakeDatFile(fullname)

[fpath,fname,ftype] = fileparts(fullname);
disp(['Reading file ',fname,ftype]);

switch ftype
    
    case '.plx'
        channels = (32:55);
        [~, npoints, ~, ~] = plx_ad_gap_info([fpath,'\',fname,ftype],0);        
        ad=zeros(length(channels),npoints);
        
        for cc=1:length(channels)
            [~, ~, ~, ~, ad(cc,:)] = plx_ad([fpath,'\',fname,ftype],channels(cc));
        end
        ad = int16(ad);
        fidP = fopen([fname,'.dat'],'w');
        disp('Writing dat file ...');
        fwrite(fidP,ad,'int16');
        fclose(fidP);
        fclose all;
        clear ad
        
    case '.ns5'
        
%       NS6 = openNSx([fpath,'\',fname,ftype],'noread');
%       chn = NS6.MetaTags.ChannelCount;
        NS6 = openNSx([fpath,'\',fname,ftype]);
        pause (1);
        fidN = fopen([fname,'.dat'], 'a');
        disp('Writing dat file ...');
        fwrite(fidN, NS6.Data,'int16');
        fclose(fidN);
        pause (1);
        clear NS6;
end

fclose all; % close all;

if ~strcmp(pwd,fpath)
    movefile ([fname,'.dat'],[fpath,'\',fname,'.dat']);
end

disp('Binary file created');
