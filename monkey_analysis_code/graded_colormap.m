function newcolr = graded_colormap(colr,Ngrades)
% Create graded colormap
% If existing colormap is longer than one color/column, graded_colormap creates
% grades for each existing color/column.

mfactor = flip(linspace(0,1,Ngrades+1));
mfactor = mfactor(1:end-1);
newcolr = [];
for i = 1:size(colr,1)
    newcolr = [newcolr ; repmat(colr(i,:),Ngrades,1).*mfactor'];
end