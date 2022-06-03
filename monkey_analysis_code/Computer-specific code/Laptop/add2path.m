%% folder for path

% choose which experiment code to load
experiment = input(['Which experiment do you want to analyze? \n' '1 for Motion Cueing, 0 for Perturbations: \n']);
% Add Data folder
addpath(genpath('C:\Users\soopa\Documents\MATLAB\Subjects'));

% Add Code
addpath('G:\My Drive\MATLAB\Code');
addpath('G:\My Drive\MATLAB\Code\regress');

if experiment % choose which code to load

addpath(genpath('G:\My Drive\MATLAB\Code\Single firefly with motion cuing'));
rmpath(genpath('G:\My Drive\MATLAB\Code\Single firefly with perturbations'));
else
addpath(genpath('G:\My Drive\MATLAB\Code\Single firefly with perturbations'));
rmpath(genpath('G:\My Drive\MATLAB\Code\Single firefly with motion cuing'));
end

addpath(genpath('G:\My Drive\MATLAB\Code\Computer-specific code\Laptop'));

addpath(genpath('G:\My Drive\MATLAB\Code\Sync folders code'));

addpath(genpath('G:\My Drive\MATLAB\Code\bads-master'));

% Set default Tick Direction
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');

clear experiment