function cmap = bluemap(N)
%% blue to white map
cmap = [linspace(0,0.8,N)' linspace(0,0.8,N)' ones(N,1)];
