%% Vergence Analysis

% this is a code written by Nastaran Arfaei at angelaki Lab NYU to test
% whether monkeys fuse the image properly when using the refitted VR
% headset. 

%% step1
% housekeeping
close all; clear all; clc;

% Path stuff
addpath('Z:\Projects\Monkey_Unity_VR\data');
addpath('Z:\Projects\Monkey_Unity_VR\code');
addpath('Z:\Projects\Monkey_Unity_VR\code\cbrewer');

% color stuff
reds = cbrewer('seq', 'Reds', 10);

% Directory stuff
human_directory = 'Z:\Projects\Monkey_Unity_VR\data';

%% THINGS TO CHANGE

newtable = [calibrationdata160;calibrationdata1370;calibrationdata4928;calibrationdata5325;calibrationdata581;calibrationdata7677;calibrationdata936];
%newtable=[calibrationdata6030];
monkey_name = 'Jimmy';
date_name = '20210622';
name='calibration_data_1101.txt';
OurData= newtable;
%% clean up the first 0.5 second
idx = [];

for  ii = 2:length(OurData.Mode)
    if OurData.Mode(ii)==2 && OurData.Mode(ii-1)==3;
        idx=[idx, ii];
    end
end

for  ii = length(idx):-1:1
    OurData(idx(ii):idx(ii)+80,:)=[];
end
%OurData(0:80,:)=[];

OurData(OurData.Status=='none',:)=[];

%%

Mode=table2array(OurData(:, 6));
Locations= table2array(OurData(:, 7:8));
Index=table2array(OurData(:, 5));
confidence=table2array(OurData(:, 4));

LeftX=[];
RightX=[];
BioncularX=[];
LeftY=[];
RightY=[];
BioncularY=[];

LeftXSTE=[];
RightXSTE=[];
BioncularXSTE=[];
LeftYSTE=[];
RightYSTE=[];
BioncularYSTE=[];

for ii = 0:10
    thisLeftX= Locations(Index==ii & Mode==1 & confidence>0.9 ,1);
    thisLeftY= Locations(Index==ii & Mode==1 & confidence>0.9,2);
    thisRightX=Locations(Index==ii & Mode==0 & confidence>0.9,1);
    thisRightY=Locations(Index==ii & Mode==0 & confidence>0.9,2);
    thisBothX= Locations(Index==ii & Mode==2 & confidence>0.9,1);
    thisBothY= Locations(Index==ii & Mode==2 & confidence>0.9,2);
    
    LeftX(ii+1)= nanmean(thisLeftX);
    RightX(ii+1)= nanmean(thisRightX);
    BioncularX(ii+1)= nanmean(thisBothX);
    LeftY(ii+1)= nanmean(thisLeftY);
    RightY(ii+1)= nanmean(thisRightY);
    BioncularY(ii+1)= nanmean(thisBothY);
    
    %record standard errors too to see if significant
    LeftXSTE(ii+1)= nanstd(thisLeftX)/sqrt(42);
    RightXSTE(ii+1)= nanstd(thisRightX)/sqrt(42);
    BioncularXSTE(ii+1)= nanstd(thisBothX)/sqrt(42);
    LeftYSTE(ii+1)=  nanstd(thisLeftY)/sqrt(42);
    RightYSTE(ii+1)= nanstd(thisRightY)/sqrt(42);
    BioncularYSTE(ii+1)= nanstd(thisBothY)/sqrt(42);
    
end
%%

figure(1)
y = [LeftX; RightX];
error= [LeftXSTE; RightXSTE];
error=error';
b = bar(y', 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(y');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',y',error,'k','linestyle','none');
hold off


figure(2)
y = [LeftY; RightY];
error= [LeftYSTE; RightYSTE];
error=error';
b = bar(y', 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(y');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',y',error,'k','linestyle','none');
hold off



figure(1)
y = [LeftX; RightX;BioncularX];
error= [LeftXSTE; RightXSTE;BioncularXSTE];
error=error';
b = bar(y', 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(y');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',y',error,'k','linestyle','none');
hold off
%style



figure(2)
y = [LeftY; RightY;BioncularY];
error= [LeftYSTE; RightYSTE;BioncularYSTE];
error=error';
b = bar(y', 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(y');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',y',error,'k','linestyle','none');
hold off
%style
%%
%added 7/8 - Josh

syms k_1 k_2 p_1 p_2 real
syms r x y

distortionX = subs(x * (1 + k_1 * r^2 + k_2 * r^4) + 2 * p_1 * x * y + p_2 * (r^2 + 2 * x^2), r, sqrt(x^2 + y^2));
distortionY = subs(y * (1 + k_1 * r^2 + k_2 * r^4) + 2 * p_2 * x * y + p_1 * (r^2 + 2 * y^2), r, sqrt(x^2 + y^2));

syms X Y Positive

eq1 = X == distortionX;
eq2 = Y == distortionY;

parameters = [k_1 k_2 p_1 p_2];
parameterValues = [0.2 0 0 0];
eq1 = expand(subs(eq1, parameters, parameterValues));
eq2 = expand(subs(eq2, parameters, parameterValues));

Result1 = solve([eq1, eq2], [x, y], 'MaxDegree', 3, 'Real', true, 'ReturnConditions', true);

%%
this=table2array(OurData(:, 7:9));
figure(3)
scatter(LeftX,LeftY);
hold on
scatter(RightX,RightY);
scatter(BioncularX,BioncularY);

figure(3)
scatter(table2array(locations(:,2)),table2array(locations(:,3)),'k','.');
grid on

polLocations = atan2d(Locations(:,1),Locations(:,2));
figure(4)
plot(table2array(OurData(50000:100000,2)),Locations(50000:100000,1))
hold on
plot(table2array(OurData(50000:100000,2)),Locations(50000:100000,2))

%%

OurData=calibrationdata312;

Mode=table2array(OurData(:, 6));
Locations= table2array(OurData(:, 7:8));
Index=table2array(OurData(:, 5));
confidence=table2array(OurData(:, 4));

LeftX1=[];
RightX1=[];
BioncularX1=[];
LeftY1=[];
RightY1=[];
BioncularY1=[];

LeftXSTE1=[];
RightXSTE1=[];
BioncularXSTE1=[];
LeftYSTE1=[];
RightYSTE1=[];
BioncularYSTE1=[];

for ii = 1:10
    thisLeftX1= Locations(Index==ii & Mode==1 & confidence>0.9,1);
    thisLeftY1= Locations(Index==ii & Mode==1 & confidence>0.9,2);
    thisRightX1=Locations(Index==ii & Mode==0 & confidence>0.9,1);
    thisRightY1=Locations(Index==ii & Mode==0 & confidence>0.9,2);
    thisBothX1= Locations(Index==ii & Mode==2 & confidence>0.9,1);
    thisBothY1= Locations(Index==ii & Mode==2 & confidence>0.9,2);
    
    LeftX1(ii)= nanmean(thisLeftX1);
    RightX1(ii)= nanmean(thisRightX1);
    BioncularX1(ii)= nanmean(thisBothX1);
    LeftY1(ii)= nanmean(thisLeftY1);
    RightY1(ii)= nanmean(thisRightY1);
    BioncularY1(ii)= nanmean(thisBothY1);
    
    %record standard errors too to see if significant
    LeftXSTE1(ii)= nanstd(thisLeftX1)/sqrt(length(thisLeftX1~=NaN));
    RightXSTE1(ii)= nanstd(thisRightX1)/sqrt(length(thisRightX1~=NaN));
    BioncularXSTE1(ii)= nanstd(thisBothX1)/sqrt(length(thisBothX1~=NaN));
    LeftYSTE1(ii)=  nanstd(thisLeftY1)/sqrt(length(thisLeftY1~=NaN));
    RightYSTE1(ii)= nanstd(thisRightY1)/sqrt(length(thisRightY1~=NaN));
    BioncularYSTE1(ii)= nanstd(thisBothY1)/sqrt(length(thisBothY1~=NaN));
    
end

figure(1)
y = [(LeftX+LeftX1)/2; (RightX+RightX1)/2];
bar(y')


% figure(2)
% y = [LeftX; RightX;BioncularX];
% bar(y')


figure(2)
y = [(LeftY+LeftY1)/2; (RightY+RightY1)/2];
bar(y')

figure(3)
scatter(this(Mode==1 ,1),this(Mode==1 ,2),this(Mode==1 ,3));
hold on
scatter3(this(Mode==0 ,1),this(Mode==0 ,2),this(Mode==0 ,3),'.');
scatter3(this(Mode==2 ,1),this(Mode==2 ,2),this(Mode==2 ,3));










((LeftY+LeftY1)/2)-((RightY+RightY1)/2)






BioncularYSTEPlot= nanstd(LeftY(1),LeftY1(1))/sqrt(2);





this=table2array(OurData(:, 7:9));
figure(3)
scatter3(this(Mode==1 ,1),this(Mode==1 ,2),this(Mode==1 ,3));
hold on
scatter3(this(Mode==0 ,1),this(Mode==0 ,2),this(Mode==0 ,3));



NormRight= table2array(OurData(:, 17:19));

figure(5)
scatter3(NormRight(Mode==1 ,1),NormRight(Mode==1 ,2),NormRight(Mode==1 ,3),'.');
hold on
scatter3(NormRight(Mode==0 ,1),NormRight(Mode==0 ,2),NormRight(Mode==0 ,3),'.');









