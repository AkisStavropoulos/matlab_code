%% Traveled distance for a given control as a function of tau

% set parameters
x = 400;
theta = 80;
T = 8.5;
dt = 1/60;

% set taus and vmax
N = 8;
tau = linspace(0.5,4,N); % [0.5 1 1.5 2 2.5 3 3.5 4];

vmax = (x/T)*( 1 ./ (-1 + 2*(tau./T) .* log((1 + exp(T./tau))/2))) ;
wmax = (theta/T)*( 1 ./ (-1 + 2*(tau./T) .* log((1 + exp(T./tau))/2))) ;

% set target distance
x_tar = 0; y_tar = x;

% generate appropriate control input for a given tau
ii = 6;
trltime =  2.*tau(ii).*acosh( exp( y_tar./(2.*tau(ii).*vmax(ii)) ) );
sw = switchtime(tau(ii), trltime);

% generate trajectories
for i = 1:N

[err,x_sub,y_sub{i},v,w,v_final,w_final] = gen_sim_traj(x_tar,y_tar,tau(i),vmax(i),wmax(i),dt,trltime,sw,0);
end

figure;plot(tau,cellfun(@(x) x(end),y_sub)); hline(y_tar,'k--'); vline(tau(ii),'r-');
xlabel('time constant [s]'); ylabel('traveled distance [cm]'); title('distance traveled for a given control as a function of \tau');