function [t_saccade,sac_mag,peakvel] = DetectSaccades(zle,yle,zre,yre,prs,dt,plt3)

%% detect saccade times - based on Larsson et al., 2013
% dt = prs.dt;
plt = 0;
plt1 = 0;
plt2 = 0;

% Define filter for smoothing of velocity magnitude
sig = 5*prs.filtwidth; %filter width
sz = 10*prs.filtsize; %filter size
t2 = linspace(-sz/2, sz/2, sz);
h = exp(-t2.^2/(2*sig^2));
h = h/sum(h); % normalise filter to ensure area under the graph of the data is not altered

% % take derivative of eye position = eye velocity
% dze = diff(0.5*(ch.zle + ch.zre));  dye = diff(0.5*(ch.yle + ch.yre));
% de = sqrt(dze.^2 + dye.^2); % speed of eye movement
% da = sqrt(diff(dze).^2 + diff(dye).^2); % acceleration of eye movement
% de_smooth = medfilt1(de,30)/dt;     da_smooth = medfilt1(da,30)/dt;
% de_smooth = conv(de,h,'same')/dt;   da_smooth = conv(da,h,'same')/dt;

% define filter for derivative of position and velocity x,y
n = round((1/dt)/250)*2; 
b = 1/60;
g = b*[ones(1,n) 0 -ones(1,n)];

% approximate saccadic intervals
dze = conv(0.5*(zle + zre),g,'same')/dt; dye = conv(0.5*(yle + yre),g,'same')/dt;
aze = conv(dze,g,'same');                   aye = conv(dye,g,'same');
de = sqrt(dze.^2 + dye.^2); % speed of eye movement
da = sqrt(aze.^2 + aye.^2); % acceleration of eye movement
de_smooth = de; % conv(de,h,'same');           
da_smooth = conv(da,h,'same');

sd_aze = nanstd(aze);                       sd_aye = nanstd(aye);           sd_da = nanstd(da);

lamda = prs.lambda;
Iz = abs(aze) > lamda*sd_aze;               Iy = abs(aye) > lamda*sd_aye;   Izy = abs(da) > lamda*sd_da;
I = Iz + Iy > 0;
min_sac_dur = prs.saccade_duration(1);
max_sac_dur = prs.saccade_duration(2);
sac_refr_period = prs.sac_refr_period; % seconds
wn = .12; %seconds, time window
max_back_iter = 200; % No. of max allowable backwards iterations to find saccade onset

% Cut-off last part of trial
trl_cutoff = ceil(prs.trl_cutoff/dt);
if any(I(end-trl_cutoff:end))
    I(end-trl_cutoff:end) = 0;
    indcut = find(diff(I(1:end-trl_cutoff)) == 1,1,'last');
    I(indcut-1:end) = 0;
end

% another smoothed-derivatives technique (Savitzky-Golay filtering, Nystrom & Holmqvist, 2010)
if 0
order = 2;
frames = 81;
[b,g] = sgolay(order,frames);

eyepos = [(0.5*(ch.zle + ch.zre)).^2   (0.5*(ch.yle + ch.yre)).^2];
        de_sg = [];  da_sg = [];
p=1;    de_sg(:,1) = conv(eyepos(:,1), factorial(p)/(-dt)^p * g(:,p+1), 'same');
        de_sg(:,2) = conv(eyepos(:,2), factorial(p)/(-dt)^p * g(:,p+1), 'same');
p=2;    da_sg(:,1) = conv(eyepos(:,1), factorial(p)/(-dt)^p * g(:,p+1), 'same');
        da_sg(:,2) = conv(eyepos(:,2), factorial(p)/(-dt)^p * g(:,p+1), 'same');
de_sg = sqrt(sum(de_sg.^2,2));          da_sg = sqrt(sum(da_sg.^2,2));
end
%
% check saccade refractory period
reset = 1;  n = 1;
while reset
    a = 1;
    if I(n) == 0 
        while I(n+a) == 0
            a = a+1;
            if n+a == length(I);    reset = 0;  break;  end
        end
        if a*dt < sac_refr_period
            I(n:n+a) = 1;
        end
    end
    n = n+a;
    if n == length(I);  reset = 0;  end
end

% check saccade minimum duration
reset = 1;  n = 1;
while reset
    a = 1;
    if I(n) == 1
        while I(n+a) == 1
            a = a+1;
            if n+a == length(I);    reset = 0;  break;  end
        end
        if a*dt < min_sac_dur
            I(n:n+a) = 0;
        end
    end
    n = n+a;
    if n == length(I);  reset = 0;  end
end

% check saccade maximum duration
reset = 1;  n = 1;
while reset
    a = 1;
    if I(n) == 1
        while I(n+a) == 1
            a = a+1;
            if n+a == length(I);    reset = 0;  break;  end
        end
        if a*dt > max_sac_dur
            I(n:n+a) = 0;
            if n+a == length(I);    reset = 0;  break;  end
        end
    end
    n = n+a;
    if n == length(I);  reset = 0;  end
end

if plt2
    figure;plot(X);grid on;ylim([-50 50]);hold on;plot(I*20,'k')
    vline(t.beg/dt,'k');vline(t.end/dt,'r');    
end
figwidth = 2*10^4;

% Cut off first part of trial
I(1:wn/dt+1) = 0;

% saccade onset and offset detection
temp_on = find([0 ; diff(I)] == 1);  temp_off = find([0; diff(I)] == -1);

if ~isempty(temp_on)
    if temp_on(1) > temp_off(1)
        temp_off = temp_off(2:end);
    end
    temp_on = temp_on(1:length(temp_off));
    t_1 = .006; % seconds
    t_2 = .008; % seconds
    delta = 50; % degrees
    beta = 30; % degrees
    for n = 1:length(temp_on)
        [peakvel(n),peakind(n)] = max(de_smooth(temp_on(n):temp_off(n)));
        peakind(n) = peakind(n) + temp_on(n) - 1;
        gamma = (1/3)*(atan2d(dze(peakind(n)-1),dye(peakind(n)-1))...
            + atan2d(dze(peakind(n)),dye(peakind(n)))...
            + atan2d(dze(peakind(n)+1),dye(peakind(n)+1)));
        
        % find saccade onset
        reset = 1;  i = 1;  count1 = 0; count2 = 0; count3 = 0; a = []; SacOnset1 = []; SacOnset2 = []; SacOnset3 = [];
        while reset
            %
            a(i) = atan2d(dze(peakind(n)-i),dye(peakind(n)-i));
            SacOnset1(i) = abs(a(i) - gamma) > delta;
            %
            if i>1
                eps = a(i) - a(i-1);
                SacOnset2(i) = abs(eps) > beta;
            else
                SacOnset2(i) = 0;
            end
            %
            pre_de = nanmean(abs(de(peakind(n)-i-wn/dt:peakind(n)-i)));
            t_de = de(peakind(n)-i);
            SacOnset3(i)  = abs(t_de) < pre_de;
            %
            if SacOnset1(i) == 1 % crit 1: deviation from main direction
                count1 = count1+1;
                if count1 >= t_1/dt
                    sac_on(n) = peakind(n)-i + count1;  if plt; vline(sac_on(n),'b');   xlim([sac_on(n)-figwidth sac_on(n)+figwidth]);  end
                    reset = 0;
                end
            else
                count1 = 0;
            end
            if SacOnset2(i) == 1 % crit 2: inconsistent sample-to-sample direction
                count2 = count2+1;
                if count2 >= t_2/dt && de_smooth(peakind(n)-i) < .2*peakvel(n)
                    sac_on(n) = peakind(n)-i;   if plt; vline(sac_on(n),'b');   xlim([sac_on(n)-figwidth sac_on(n)+figwidth]);  end
                    reset = 0;
                end
            else
                count2 = 0;
            end
            if SacOnset3(i) == 1 % crit 3: smaller displacement at edges of saccade
                count3 = count3+1;
                if count3 >= t_2/dt && de_smooth(peakind(n)-i) < .2*peakvel(n)
                    sac_on(n) = peakind(n)-i;   if plt; vline(sac_on(n),'b');   xlim([sac_on(n)-figwidth sac_on(n)+figwidth]);  end
                    reset = 0;
                end
            else
                count3 = 0;
            end
            
            i = i+1;
            if (peakind(n)-i-wn/dt) < 1; reset = 0; sac_on(n) = nan; end
        end
        % find saccade offset
        reset = 1;  i = 1;   count1 = 0;    count2 = 0; count3 = 0; a = []; SacOffset1 = []; SacOffset2 = []; SacOffset3 = [];
        while reset
            %
            a(i) = atan2d(dze(peakind(n)+i),dye(peakind(n)+i));
            SacOffset1(i) = abs(a(i) - gamma) > delta;
            %
            if i>1
                eps = a(i) - a(i-1);
                SacOffset2(i) = abs(eps) > beta;
            else
                SacOffset2(i) = 0;
            end
            %
            post_de = nanmean(abs(de(peakind(n)+i:peakind(n)+i+wn/dt)));
            t_de = de(peakind(n)+i);
            SacOffset3(i)  = abs(t_de) < pre_de;
            %
            if SacOffset1(i) == 1 % crit 1: deviation from main direction
                count1 = count1+1;
                if count1 >= t_1/dt
                    sac_off(n) = peakind(n)+i - count1; if plt; vline(sac_off(n),'m');   xlim([sac_off(n)-figwidth sac_off(n)+figwidth]);   end
                    reset = 0;
                end
            else
                count1 = 0;
            end
            if SacOffset2(i) == 1 % crit 2: inconsistent sample-to-sample direction
                count2 = count2+1;
                if count2 >= t_2/dt && de_smooth(peakind(n)+i) < .2*peakvel(n)
                    sac_off(n) = peakind(n)+i; if plt; vline(sac_off(n),'m');   xlim([sac_off(n)-figwidth sac_off(n)+figwidth]);    end
                    reset = 0;
                end
            else
                count2 = 0;
            end
            if SacOffset3(i) == 1 % crit 3: smaller displacement at edges of saccade
                count3 = count3+1;
                if count3 >= t_2/dt && de_smooth(peakind(n)+i) < .2*peakvel(n)
                    sac_off(n) = peakind(n)+i; if plt;  vline(sac_off(n),'m');   xlim([sac_off(n)-figwidth sac_off(n)+figwidth]);   end
                    reset = 0;
                end
            else
                count3 = 0;
            end
            
            i = i+1;
            if (peakind(n)+i+wn/dt) > length(de); reset = 0; sac_off(n) = peakind(n)+i; end
        end
    end
    temp_on = temp_on; temp_off = temp_off;
    TEMP = temp_on(:) - sac_on(:);
    maxindx = find(TEMP > max_back_iter);
    falseindx = find((sac_off(:)-sac_on(:))*dt > max_sac_dur);
    nanindx = find(isnan(sac_on));
    rmindx = unique([falseindx(:) ; maxindx(:) ; nanindx(:)]);
    sac_on(rmindx) = [];  sac_off(rmindx) = [];
    temp_on(rmindx) = []; temp_off(rmindx) = [];
    
    % check that the eye has changed position (valid saccade)
    ON = [sac_on(:) temp_on(:)];    OFF = [sac_off(:) temp_off(:)];
    
    minmag = prs.min_sac_mag; % minimum saccade magnitude
    overlap_pct = prs.overlap_pct;
    Nsamples = round(prs.peri_sac_dur/dt);
    skipsamp_pre = 0;     skipsamp_post = 0;    sac_crit = [];  SacIndx = [];   sac_mag = [];
    for n = 1:size(ON,1)
        for k = 1:size(ON,2)
            % pre-saccade (w/ fixes for first saccades)
            if ~isempty(SacIndx); N = SacIndx(end); else; N = n-1; end
            if N > 0
                if (ON(n,k)-Nsamples) - OFF(N,k) < - 20
                    skipsamp_pre = abs((ON(n,k)-Nsamples) - OFF(N,k));
                else; skipsamp_pre  = 0;
                end
            else
                if ON(n,k)-Nsamples > 0
                    skipsamp_pre  = 0;
                else; skipsamp_pre  = -(ON(n,k)-Nsamples) + 1;
                end
            end
            
            LEpos_pre = [zle(ON(n,k)-Nsamples+skipsamp_pre:ON(n,k)) yle(ON(n,k)-Nsamples+skipsamp_pre:ON(n,k))];
            LEpos_pre = [nan(skipsamp_pre,2) ; LEpos_pre];
            LEpre = [];
            for zz = 1:size(LEpos_pre,2)
                x = [1:size(LEpos_pre,1)];  XX = [ones(length(x),1) x(:)];
                y = LEpos_pre(:,zz);
                [c]=regress(y,XX);
                LEpre(:,zz) = LEpos_pre(:,zz) - XX*c + XX(end,:)*c;
            end
            REpos_pre = [zre(ON(n,k)-Nsamples+skipsamp_pre:ON(n,k)) yre(ON(n,k)-Nsamples+skipsamp_pre:ON(n,k))];
            REpos_pre = [nan(skipsamp_pre,2) ; REpos_pre];
            REpre  = [];
            for zz = 1:size(REpos_pre,2)
                x = [1:size(REpos_pre,1)];  XX = [ones(length(x),1) x(:)];
                y = REpos_pre(:,zz);
                [c]=regress(y,XX);
                REpre(:,zz) = REpos_pre(:,zz) - XX*c + XX(end,:)*c;
            end
            BEpos_pre = [0.5*(LEpre(:,1) + REpre(:,1))  0.5*(LEpre(:,2) + REpre(:,2))];
            % sanity check
            if plt1
                figure;plot(REpos_pre(:,zz));hold on;hold on;plot(XX*c,'r');plot(REpre(:,zz),'b')
            end
            %
            % post-saccade (w/ fixes for last saccades)
            if n < size(ON,1)
                if (OFF(n,k)+Nsamples) - ON(n+1,k) > 20
                    skipsamp_post = (OFF(n,k)+Nsamples) - ON(n+1,k);
                else; skipsamp_post = 0; end
            else
                skipsamp_post = 0;
            end
            
            LEpos_post = [zle(OFF(n,k):OFF(n,k)+Nsamples-skipsamp_post) yle(OFF(n,k):OFF(n,k)+Nsamples-skipsamp_post)];
            LEpos_post = [LEpos_post ; nan(skipsamp_post,2)];
            LEpost = [];
            for zz = 1:size(LEpos_post,2)
                x = [1:size(LEpos_post,1)];  XX = [ones(length(x),1) x(:)];
                y = LEpos_post(:,zz);
                [c]=regress(y,XX);
                LEpost(:,zz) = LEpos_post(:,zz) - XX*c + XX(1,:)*c;
            end
            REpos_post = [zre(OFF(n,k)+skipsamp_post:OFF(n,k)+Nsamples) yre(OFF(n,k)+skipsamp_post:OFF(n,k)+Nsamples)];
            REpos_post = [REpos_post ; nan(skipsamp_post,2)];
            REpost  = [];
            for zz = 1:size(REpos_post,2)
                x = [1:size(REpos_post,1)];  XX = [ones(length(x),1) x(:)];
                y = REpos_post(:,zz);
                [c]=regress(y,XX);
                REpost(:,zz) = REpos_post(:,zz) - XX*c + XX(1,:)*c;
            end
            BEpos_post = [0.5*(LEpost(:,1) + REpost(:,1))  0.5*(LEpost(:,2) + REpost(:,2))];
            BEpost(n,k) = sqrt(sum(nanmean(BEpos_post,1).^2));
            % sanity check
            if plt1
                hold on;plot(REpos_post(:,zz));hold on;hold on;plot(XX*c,'r');plot(REpost(:,zz),'k')
            end
            % saccade magnitude
            LEsac_mag = sqrt((nanmean(LEpos_pre(:,1))-nanmean(LEpos_post(:,1))).^2 + (nanmean(LEpos_pre(:,2))-nanmean(LEpos_post(:,2))).^2);
            REsac_mag = sqrt((nanmean(REpos_pre(:,1))-nanmean(REpos_post(:,1))).^2 + (nanmean(REpos_pre(:,2))-nanmean(REpos_post(:,2))).^2);
            BEsac_mag = sqrt((nanmean(BEpos_pre(:,1))-nanmean(BEpos_post(:,1))).^2 + (nanmean(BEpos_pre(:,2))-nanmean(BEpos_post(:,2))).^2);
            sac_mag(n,k) = BEsac_mag;
            % Pre- and Post-saccade position variability
            pre_sd_z = prctile(BEpos_pre(:,1),overlap_pct);
            pre_sd_y = prctile(BEpos_pre(:,2),overlap_pct);
            post_sd_z = prctile(BEpos_post(:,1),overlap_pct);
            post_sd_y = prctile(BEpos_post(:,2),overlap_pct);
            % test overlap (requires at least one non-overlapping dimension)
            ind_z = any(sum( post_sd_z <= [pre_sd_z; flip(pre_sd_z)] ).*sum( post_sd_z >= [pre_sd_z; flip(pre_sd_z)] ));
            ind_y = any(sum( post_sd_y <= [pre_sd_y; flip(pre_sd_y)] ).*sum( post_sd_y >= [pre_sd_y; flip(pre_sd_y)] ));
            overlap(n) = and(ind_z,ind_y);
        end
        if sac_mag(n,1) > minmag  && ~overlap(n) % abs(diff(BEpost(n,:))) < minmag
            sac_crit = [sac_crit sac_on(n)];
            SacIndx = [SacIndx n];
            if plt2; vline(sac_crit(end),'g');   xlim([sac_on(n)-figwidth sac_on(n)+figwidth]);   end
        end
    end
    if ~isempty(sac_crit)
        t_saccade = sac_crit(:)*dt;
        sac_on = sac_on(SacIndx);       sac_off = sac_off(SacIndx);
        peakvel = peakvel(SacIndx);     sac_mag = sac_mag(SacIndx,1);
    else
        t_saccade = [];
        sac_on = [];       sac_off = [];
        peakvel = [];      sac_mag = [];
    end
else
    t_saccade = [];
    sac_on = [];       sac_off = [];
    peakvel = [];      sac_mag = [];
end
%% Test saccade detection criteria
if plt3
figure;plot(sac_mag,peakvel,'.');   xlabel('saccade magnitude');    ylabel('Peak Velocity');    xlim([0 60]);  %ylim([0 120]);
XX = [ones(length(sac_mag),1) sac_mag(:)];
y = peakvel(:);
c = regress(y,XX);
res_err = peakvel(:) - XX*c;
hold on;    plot(XX(:,end),XX*c,'r'); title(['slope = ' num2str(c(end),3) ', \sigma_{res} = ' num2str(mean(abs(res_err)))]);

OutlierIndx = abs(res_err) > 15;
end
