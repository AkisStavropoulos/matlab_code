function[ChanList]=AddEventCh(fhand)
% SONCHANLIST returns a structure with details of active channels in a SMR file
% 

if fhand<1
    ChanList=[];
    return;
end

AcChan=0;

for i=1:CEDS64MaxChan(fhand)
    kind=CEDS64ChanType(fhand,i);
    if(kind>0)                   % Only look at channels that are active
        AcChan=AcChan+1;
        ChanList(AcChan).number=i;
        ChanList(AcChan).kind=kind;
        [~, tt]=CEDS64ChanTitle(fhand,i);
        ChanList(AcChan).title= tt;
        [~, cc]=CEDS64ChanComment(fhand,i);
        ChanList(AcChan).comment=cc;
%         ChanList(AcChan).phyChan=c.phyChan;
    end
end

            
            
            
