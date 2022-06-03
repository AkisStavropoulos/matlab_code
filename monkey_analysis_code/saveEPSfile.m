function saveEPSfile(foldername,filename)

cd(foldername);
filename = [filename '.eps'];
print(filename,'-depsc2','-tiff','-r300','-painters');
disp([filename ' saved.']);