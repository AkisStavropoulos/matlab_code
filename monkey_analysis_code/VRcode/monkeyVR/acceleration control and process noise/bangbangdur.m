function traveltime = bangbangdur(x,tau,vmax)
%% calculate duration of bang-bang trajectory

traveltime = 2.*tau.*acosh( exp( x./(2.*tau.*vmax) ) );
