%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Jean-Paul Noel, NYU Nov 2020, Extracting coherence information for
% traveling wave
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Housekeeping
close all; clear all; clc; 

% saving location; 
save_folder = 'C:\Users\lab\Documents\Savin-Angelaki\Coherence'; 

% Sessions to be extracted
animal = [53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, ...
    51, 51, 51, 51, 51, 51,...
    44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44] ; 

session_num = [110, 108, 93, 123, 130, 126, 122, 124, 116, 114, 113, 92, 120, 136, 134, 133, 132, 128, 115, 111, ...
120, 43, 42, 41, 40, 38,...
221, 220, 219, 218, 217, 216, 215, 213, 212, 211, 209, 208, 207, 206, 205, 203, 202, 200, 190, 189];  

for i = 1:size(animal, 2)
    % I just think it may be cleaner to have each session separately, so I
    % re-start this every time.
    experiment = experiment('monkey-firefly'); 
    
    try
    % Do all the hard work
    experiment.AddSessions_JP(animal(i), session_num(i), {'behv','units', 'lfps','pop'});
    
    % Just keep the stuff I want (for space)
    session = experiment.sessions(end);
    
    % Save
    filename = fullfile(save_folder, [num2str(animal(i)), 's', num2str(session_num(i)), '.mat']); 
    save(filename, 'session', '-v7.3');
    catch
    end
    
    clearvars -except animal session_num i save_folder
    
end
