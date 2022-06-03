function Fcontent = rmdotfiles(Fcontent)

ind = [];
for i = 1:length(Fcontent)
    if strcmp(Fcontent(i).name,'.') || strcmp(Fcontent(i).name,'..')
        ind = [ind i];
    end
end

Fcontent(ind) = [];
