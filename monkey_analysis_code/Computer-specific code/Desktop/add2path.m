%% folder for path

% choose code to load
experiment = input(['Which experiment do you want to analyze? \n' '1 Motion Cueing, 2 Savin-Angelaki, 3 monkey VR, 4 humanVR: \n']);
% Add Data folder
addpath(genpath('C:\Users\ges6\Documents\MATLAB\Data'));

% Add Code
addpath('G:\My Drive\MATLAB\Code');

if experiment==1 % choose which code to load
    addpath(genpath('G:\My Drive\MATLAB\Code\Single firefly with motion cuing'));
    rmpath(genpath('G:\My Drive\MATLAB\Code\Savin-Angelaki'));
    
elseif experiment==2
    addpath(genpath('G:\My Drive\MATLAB\Code\Savin-Angelaki'));
    rmpath(genpath('G:\My Drive\MATLAB\Code\Single firefly with motion cuing'));

elseif experiment==3
    addpath('G:\My Drive\MATLAB\Code\VRcode');
    addpath(genpath('G:\My Drive\MATLAB\Code\VRcode\monkeyVR'));
    addpath(genpath('G:\My Drive\MATLAB\Code\Savin-Angelaki'));
    rmpath(genpath('G:\My Drive\MATLAB\Code\Single firefly with motion cuing'));
%     rmpath(genpath('G:\My Drive\MATLAB\Code\Savin-Angelaki'));
elseif experiment==4
    addpath('G:\My Drive\MATLAB\Code\VRcode');
    addpath(genpath('G:\My Drive\MATLAB\Code\VRcode\humanVR'));
    addpath(genpath('G:\My Drive\MATLAB\Code\Savin-Angelaki'));
    rmpath(genpath('G:\My Drive\MATLAB\Code\Single firefly with motion cuing'));
%     rmpath(genpath('G:\My Drive\MATLAB\Code\Savin-Angelaki'));
end

addpath(genpath('G:\My Drive\MATLAB\Code\Computer-specific code\Desktop'));

addpath(genpath('G:\My Drive\MATLAB\Code\Sync folders code'));

addpath(genpath('G:\My Drive\MATLAB\Code\bads-master'));

% Set default Tick Direction
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');

% Set default figure Font Name
set(0,'defaultAxesFontName', 'MyriadPro-Regular')
set(0,'defaultTextFontName', 'MyriadPro-Regular')

clear experiment