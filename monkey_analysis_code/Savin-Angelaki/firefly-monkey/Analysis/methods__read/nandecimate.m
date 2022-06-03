% Downsample usig Matlab decimate function ingore NANs
function decimdata = nandecimate(data,dw)

if ~isnan(data)
    decimdata = decimate(data,dw);
else
    decimdata = downsample(data,dw);
    %     data(isnan(data))=0;
%     data=fillmissing(data,'linear');
%     decimdata = decimate(data,dw);
%     disp('NANdecimate');
end
    