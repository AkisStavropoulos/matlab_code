function [tar, sub, b] = TargetsVsResponses(trials,plt)

%% Extract target and response positions
startindx = arrayfun(@(x) find(~isnan(x.continuous.xmp),1), trials);
tar.x = arrayfun(@(x,i) x.prs.xfp - x.continuous.xmp(i), trials, startindx);
tar.y = arrayfun(@(x,i) x.prs.yfp - x.continuous.ymp(i), trials, startindx);
sub.x = arrayfun(@(x,i) x.continuous.xmp(end) - x.continuous.xmp(i), trials, startindx);
sub.y = arrayfun(@(x,i) x.continuous.ymp(end) - x.continuous.ymp(i), trials, startindx);

tar.r = sqrt(tar.x.^2 + tar.y.^2);
tar.th = atan2d(tar.x,tar.y);
sub.r = sqrt(sub.x.^2 + sub.y.^2);
sub.th = atan2d(sub.x,sub.y);

b.r = regress(sub.r(:),tar.r(:));
b.th = regress(sub.th(:),tar.th(:));

%% Plot

if plt
    
figure; 
subplot(1,2,1); hold on;
plot(tar.r,sub.r,'.'); plot(0:500,0:500,'k--'); axis equal; axis([0 500 0 500]); 
xlabel('target distance [cm]'); ylabel('response [cm]');

subplot(1,2,2); hold on;
plot(tar.th,sub.th,'.'); plot(-50:50,-50:50,'k--'); vline(0,'k');hline(0,'k'); axis equal; axis([-50 50 -50 50]); 
xlabel('target angle [deg]'); ylabel('response [deg]');

end