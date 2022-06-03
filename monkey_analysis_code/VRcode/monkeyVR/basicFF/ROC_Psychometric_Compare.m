clear
clc
fclose('all');

%% 
set1= s(1,1);
set2= s(1,2);
% set3= s(1,3);
% set4= s(1,4);

MasterCartSet1= set1.master_cart
MasterCartSet2= set2.master_cart
% MasterCartSet3= set3.master_cart
% MasterCartSet4= set4.master_cart

set1Filter= set1.OnOff
set2Filter= set2.OnOff
% set3Filter= set3.OnOff
% set4Filter= set4.OnOff

MasterCartSet1= MasterCartSet1(set1Filter==0.4,:)
MasterCartSet2= MasterCartSet2(set2Filter==0.4,:)
% MasterCartSet3= MasterCartSet3(set3Filter~=0.4,:)
% MasterCartSet4= MasterCartSet4(set4Filter~=0.4,:)

MasterCartMain=  vertcat(MasterCartSet1,MasterCartSet2);
%MasterCartMain= MasterCartSet1;

% get rid of trials when target was on

% Get rid of trials that monkey skipped
% if moved less than 1 meter
traveledDist=[];
DistFromTarget=[];

for i=1:length(MasterCartMain)
    traveledDist(i) = sqrt(dot((MasterCartMain(i,3:4)-MasterCartMain(i,1:2)),(MasterCartMain(i,3:4)-MasterCartMain(i,1:2))));
    DistFromTarget(i)= sqrt(dot((MasterCartMain(i,5:6)-MasterCartMain(i,1:2)),(MasterCartMain(i,5:6)-MasterCartMain(i,1:2))));
end


MasterCartMaincleaned= MasterCartMain((traveledDist./DistFromTarget)>.4 &(traveledDist./DistFromTarget)<1.6 ,:);


X_monk= MasterCartMaincleaned(:,3:4);
X_fly= MasterCartMaincleaned(:,5:6);

maxrewardwin=4;
npermutations=10;
[rewardwin, pCorrect, pCorrect_shuffled_mu] = ComputeROCFirefly(X_fly,X_monk,maxrewardwin,npermutations)

figure;
hold on;
plot(rewardwin,pCorrect,'k','linewidth',2);
plot(rewardwin,pCorrect_shuffled_mu,'Color',[.5 .5 .5],'linewidth',2);
xlabel('Hypothetical reward window (m)'); ylabel('Fraction of rewarded trials');
legend({'true','shuffled'});


figure;
hold on;
plot(pCorrect_shuffled_mu,pCorrect,'k','linewidth',2); 
plot(0:1,0:1,'--k');
xlabel('Shuffled accuracy'); ylabel('Actual accuracy');

%% Now Do this for all the datasets that you have
% You can find screen data in example file Z:\Data\Monkey2_newzdrive\Jimmy\U-probe\Mar 10 2021\neural data\Pre-processing X E

for i=1:length(behv_stats.pos_rel.x_stop)   
    MonkPosX(i) = trials_behv(i).continuous.xmp(end,end)/100;
    MonkPosY(i) = trials_behv(i).continuous.ymp(end,end)/100;
    TargPosX(i) = trials_behv(i).continuous.xfp(end,end)/100;
    TargPosY(i) = trials_behv(i).continuous.yfp(end,end)/100;
    OldFilter(i)= trials_behv(i).prs.fly_duration;
end

MonkPosX= MonkPosX(:,OldFilter==0.3);
MonkPosY= MonkPosY(:,OldFilter==0.3);
TargPosX= TargPosX(:,OldFilter==0.3);
TargPosY= TargPosY(:,OldFilter==0.3);

% Calculate traveled distance and distance from target at the end of the
% trial
for i=1:length(MonkPosX)
    traveledDistOld(i) = sqrt(dot([MonkPosX(i),MonkPosY(i)],[MonkPosX(i),MonkPosY(i)]));
    DistFromTargetOld(i)= sqrt(dot([TargPosX(i),TargPosY(i)],[TargPosX(i),TargPosY(i)]));
end

MonkPosX = MonkPosX(:,(traveledDistOld./DistFromTargetOld)>.4 &(traveledDistOld./DistFromTargetOld)<1.6);
MonkPosY = MonkPosY(:,(traveledDistOld./DistFromTargetOld)>.4 &(traveledDistOld./DistFromTargetOld)<1.6);
TargPosX = TargPosX(:,(traveledDistOld./DistFromTargetOld)>.4 &(traveledDistOld./DistFromTargetOld)<1.6);
TargPosY = TargPosY(:,(traveledDistOld./DistFromTargetOld)>.4 &(traveledDistOld./DistFromTargetOld)<1.6);

X_flyOld= [TargPosX', TargPosY' ];
X_monkOld= [MonkPosX', MonkPosY' ];

maxrewardwin=4;
npermutations=10;
[rewardwinOld, pCorrectOld, pCorrect_shuffled_muOld] = ComputeROCFirefly(X_flyOld,X_monkOld,maxrewardwin,npermutations)


figure(1)
hold on;
plot(rewardwinOld,pCorrectOld,'r','linewidth',2);
plot(rewardwinOld,pCorrect_shuffled_muOld,'--r','linewidth',2);
xlabel('Hypothetical reward window (m)'); ylabel('Fraction of rewarded trials');
legend({'trueHeadset','shuffledHeadset','trueScreen','shuffledScreen'});


figure(2)
hold on;
plot(pCorrect_shuffled_muOld,pCorrectOld,'r','linewidth',2); 
plot(0:1,0:1,'--k');
xlabel('Shuffled accuracy'); ylabel('Actual accuracy');
legend({'Headset','ChanceLevel','Screen'});


