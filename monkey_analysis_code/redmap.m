function cmap = redmap(N)
%% red to white map
cmap = [ones(N,1) linspace(0,0.8,N)' linspace(0,0.8,N)'];
