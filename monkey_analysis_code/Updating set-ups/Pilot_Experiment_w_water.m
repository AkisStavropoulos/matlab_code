%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Jean-Paul Noel - NYU - February 2019
%
% Pilot Experiment to test a bunch of different visual stimuli and see if
% animals run in different direction (in theory because of perceived
% self-motion). The idea for stimuli comes from Dichgans & Brandt (1978)
% Figure 9, and the idea for looking at running comes from anectodal
% evidence from Mohler et al. (2005)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% housekeeping
close all; clear; clc; sca;
delete(instrfindall)
delete(instrfindall)

%% for Debugging
debug = 0;
if debug % one of these has to be set to 1 for debugging to do anything
    Gabor = 1;
    VerticalLines = 0;
    Numdots = 0;
end


%% Acquire session information
if ~debug
    try
        %enter subjects information
        repeat = 1;
        while (repeat)
            inputPrompt = {'Mouse ID:'};
            defaultAnswers = {''};
            answer = inputdlg(inputPrompt,'Run Information',1, defaultAnswers);
            [subjString] = deal(answer{:});
            if isempty(subjString) || ~ischar(subjString)
                errorMsg = errordlg('Fields require strings and all fields are required','Input Error');
                repeat = 1;
                uiwait(errorMsg);
            else
                repeat = 0;
            end
        end
        
        sex = input('Subject gender (m or f): ','s');
        data.sex=sex;
        
        % USEFUL TOOL FOR COUNTING DAYS
        %d1 = datetime('08-nov-2018'); d2 = datetime('11-feb-2019');
        %between (d1, d2), and assume that 1 month = 30 days
        
        age = input('Subject age in days: ','s');
        data.age = str2double(age);
        
        data.date = date;
        
        % create output file for subject
        progFile = [pwd,'\progress\',subjString, '_', date, '_water.mat'];
        if exist(progFile) == 2 %if progress file exists
            errorMsg = errordlg('This file already exists, try again from','Input Error');
            uiwait(errorMsg);
            return
        else %create settings if no prog file
            trials = []; %create empty Trls array
        end
    catch
        return
    end
end

%% parameters for Experiment
% making sure shuffling is actually shuffling
rng('default')

% Gabor moving - different sizes (6) 3 in one direction, and 3 in the other
% direction...
% Vertical lines moving - different velocities (6)
% Dots moving - different number of dots (6)
% Lots of dots moving - different levels of coherence (6)

conditions = 1:24; % these will be a total of 20 different conditions (we will also do a long black period at the beginning and end)
nrep = 15;
conditions = repmat(conditions, 1, nrep);
conditions = Shuffle(conditions);



trials(:, 1) = Shuffle(conditions);
velocities_stimuli = cell(size(trials, 1), 4);  % initializing a cell array 'velocities' where we are going to store the velocities for x1, x2, y1, and y2
velocities_no_stimuli = cell(size(trials, 1), 4);  % initializing a cell array 'velocities' where we are going to store the velocities for x1, x2, y1, and y2
% I decided to story velocities
% instead of position, as it's
% easier to integrate numerically
% than derive.
pre = cell(1, 4);
post = cell(1, 4);


% We will have to do 5 seconds of stimuli and and 5 seconds off, to get the
% experiment to be 1 hour
TimePerPresentation = 5; % 5 seconds on, 5 seconds off
baseline_time = 60; % 1 minute of baseline before the experiment

%% Parameters for Water Reward
water_threshold_pre_post = 20e3;
water_threshold_trial = 10e3;
water_time = '40';


%% Starting Arduino for ball position data
[mr, ~, ~, ~, ~] = InitializeArduino_JP_PILOT( RigParameters );

%% Starting Arduino for Water
comPort = 'COM6';
[s, flag] = INIT_Serial(comPort);


%% PsychtoolBox starting
try %at this point we are going to do everything within a catch loop to try to minimize annoying errors with the psychtoolbox
    % Skips errors on my macbook, for some reason I too, get errors if I don't
    % include this step. I don't use Psychtoolbox to run my experiments, so I
    % haven't digged into the reason for this..
    
    %% Common accross all types of visual presentation
    
    Screen('Preference', 'SkipSyncTests', 1);
    
    % Here we call some default settings for setting up Psychtoolbox
    PsychDefaultSetup(2);
    
    % Screen Number - this is to make sure anything is presented on the
    % projector
    screenNumber = 1; %max(Screen('Screens'));
    
    % Define black, white and grey
    white = WhiteIndex(screenNumber);
    grey = white / 2;
    black = BlackIndex(screenNumber);
    
    % hidding cursor
    if ~debug
        HideCursor;
    end
    
    % Open an on screen window
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey);
    
    
    tStart = GetSecs;
    timedout = false;
    
    %% Starting water stuff for pre
    master_water_x1 = [0];
    master_water_x2 = [0];
    master_water_y1 = [0];
    master_water_y2 = [0];
    counter_water = 2;
    distance = 0;
    
    
    %%
    while ~timedout % while loop for pre
        [ x1, y1, x2, y2 ] = ReadArduino_JP_PILOT( mr );
        pre{1} = [pre{1}, x1];
        pre{2} = [pre{2}, y1];
        pre{3} = [pre{3}, x2];
        pre{4} = [pre{4}, y2];
        
        if distance < water_threshold_pre_post % If the animal wasn't run past a certain distance
            master_water_x1(counter_water) = x1 + master_water_x1(counter_water-1); % we are adding velocities here, which gives us distance
            master_water_y1(counter_water) = y1 + master_water_y1(counter_water-1);
            master_water_x2(counter_water) = x2 + master_water_x2(counter_water-1);
            master_water_y2(counter_water) = y2 + master_water_y2(counter_water-1);
            
            x = master_water_x1(counter_water) + master_water_x2(counter_water); % this is just because of how the sensor works
            y = master_water_x1(counter_water) - master_water_x2(counter_water);
            
            distance = sqrt( (x-0)^2 + (y-0)^2); % this is to keep track not only of forward distance, but both components
            
            counter_water = counter_water + 1;
        else
            counter_water = 2; % re-starting everything for the water
            
            master_water_x1 = [0];
            master_water_x2 = [0];
            master_water_y1 = [0];
            master_water_y2 = [0];
            
            distance = 0;
            
            fprintf(s, ['1,' water_time]); %actually giving the water
            
            
        end
        
        tCurrent = GetSecs;
        if tCurrent - tStart > baseline_time
            timedout = true;
        end
        
        
    end
    
    
    % Get the size of the on screen window
    [screenXpixels, screenYpixels] = Screen('WindowSize', window);
    
    % Query the frame duration
    ifi = Screen('GetFlipInterval', window);
    
    % Get the centre coordinate of the window
    [xCenter, yCenter] = RectCenter(windowRect);
    
    % Maximum priority level
    topPriorityLevel = MaxPriority(window);
    
    % Set up alpha-blending for smooth (anti-aliased) lines
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
    %% Gabor Parameters and building
    
    % Obvious Parameters
    orientation = 90;
    contrast = 1.0;
    aspectRatio = 1.0;
    
    % number of Gabors
    nGabors = 1;
    
    % Position the Gabor
    xPos = xCenter;
    yPos = yCenter + 200;
    
    % Drift speed for the 2D global motion
    degPerSec = 360 * 4;
    degPerFrame =  degPerSec * ifi;
    
    %% Vertical Lines Parameters and Building
    
    % Grating size in pixels
    gratingSizePix = 2600;
    %gratingSizePix = 1300;
    
    % Grating frequency in cycles / pixel
    freqCyclesPerPix = 0.005;
    
    % Define Half-Size of the grating image.
    texsize = gratingSizePix / 2;
    
    % First we compute pixels per cycle rounded to the nearest pixel
    pixPerCycle = ceil(1 / freqCyclesPerPix);
    
    % Frequency in Radians
    freqRad = freqCyclesPerPix * 2 * pi;
    
    % This is the visible size of the grating
    visibleSize = 2 * texsize + 1;
    
    
    % Define our grating. Note it is only 1 pixel high. PTB will make it a full
    % grating upon drawing
    x = meshgrid(-texsize:texsize + pixPerCycle, 1);
    grating = grey * cos(freqRad*x) + grey;
    
    % Make a two layer mask filled with the background colour
    mask = ones(1, numel(x), 2) * white;
    
    % Contrast for our contrast modulation mask: 0 = mask has no effect, 1 = mask
    % will at its strongest part be completely opaque frameCounter.e. 0 and 100% contrast
    % respectively
    %contrast = 0.8;
    contrast = 1;
    
    % Place the grating in the 'alpha' channel of the mask
    mask(:, :, 2)= grating .* contrast;
    
    % Make our grating mask texture
    gratingMaskTex = Screen('MakeTexture', window, mask);
    
    % Make a black and white noise mask half the size of our grating. This will
    % be scaled upon drawing to make a "chunky" noise texture which our grating
    % will mask. Note the round function in here. For this demo we are simply
    % rounding the size to the nearest pixel, leaving PTB to do some scaling
    %   noise = rand(round(visibleSize / 2)) .* white; %noise
    noise = zeros(round(visibleSize/2)); % black
    %  noise = ones(round(visibleSize/2)); % while
    
    % Make our noise texture
    noiseTex = Screen('MakeTexture', window, noise);
    
    % Make a destination rectangle for our textures and center this on the
    % screen
    dstRect = [0 0 visibleSize visibleSize-1700];
    %dstRect = [0 0 visibleSize visibleSize];
    dstRect = CenterRect(dstRect, windowRect);
    dstRect(2) = dstRect(2) + 400;
    
    % We set PTB to wait one frame before re-drawing
    waitframes = 1;
    
    % Calculate the wait duration
    waitDuration = waitframes * ifi;
    
    % Recompute pixPerCycle, this time without the ceil() operation from above.
    % Otherwise we will get wrong drift speed due to rounding errors
    pixPerCycle = 1 / freqCyclesPerPix;
    
    %% Dot parameters and bulding
    monRefresh = 1/ifi; % frames per second
    
    mon_horizontal_cm  	= 5;
    %mon_horizontal_cm  = 38;
    %view_dist_cm 		= 50;
    view_dist_cm 		= 50;
    
    % diameter/length of side of aperture
    %apD = 50; % CHANGE THIS NUMBER TO MAKE THE AREA THE DOTS OCCUPY LARGER OR SMALLER
    apD = 50;
    
    % Everything is initially in coordinates of visual degrees, convert to pixels
    % (pix/screen) * (screen/rad) * rad/deg
    ppd = pi * windowRect(3) / atan(mon_horizontal_cm/view_dist_cm/2)  / 360;
    
    d_ppd = floor(apD/10 * ppd);
    
    
    %% starting the water stuff for experimental trials
    master_water_x1 = [0];
    master_water_x2 = [0];
    master_water_y1 = [0];
    master_water_y2 = [0];
    counter_water = 2;
    distance = 0;
    
    for t = 1:size(trials, 1)
        
        fprintf ('Trial number %d out of 360\n\n\n', t); % so you know on what trial you are.
        
        condition = conditions(t); % the particular conditions for this trial
        
        if ismember(condition, 1:6) % meaning it's part of the Gabor dimension experiment
            
            % Dimensions
            if condition == 1
                gaborDimPix = 300; % THIS IS GOING TO BE A CRITICAL DIMENSION FOR SIZE
                degPerFrameGabors = 30;
            elseif condition ==  2
                gaborDimPix = 600;
                degPerFrameGabors = 30;
            elseif condition == 3
                gaborDimPix = 900;
                degPerFrameGabors = 30;
            elseif condition == 4
                gaborDimPix = 300;
                degPerFrameGabors = -30;
            elseif condition == 5
                gaborDimPix = 600;
                degPerFrameGabors = -30;
            elseif condition == 6
                gaborDimPix = 900;
                degPerFrameGabors = -30;
            end
            
            % Sigma of Gaussian
            sigma = gaborDimPix / 6;
            
            % Spatial Frequency (Cycles Per Pixel)
            % One Cycle = Grey-Black-Grey-White-Grey i.e. One Black and One White Lobe
            numCycles = 10;
            freq = numCycles / gaborDimPix;
            
            % Build a procedural gabor texture
            gabortex = CreateProceduralGabor(window, gaborDimPix, gaborDimPix,...
                [], [0.5 0.5 0.5 0.0], 1, 0.5);
            
            % Make the destination rectangles for all the Gabors in the array
            baseRect = [0 0 gaborDimPix gaborDimPix];
            %allRects = nan(4, 1);
            allRects = CenterRectOnPointd(baseRect, xPos, yPos);
            
            % Randomise the Gabor orientations and determine the drift speeds of each gabor.
            % This is given by multiplying the global motion speed by the cosine
            % difference between the global motion direction and the global motion.
            % Here the global motion direction is 0. So it is just the cosine of the
            % angle we use. We re-orientate the array when drawing
            gaborAngles = 90;
            
            % Randomise the phase of the Gabors and make a properties matrix. We could
            % if we want have each Gabor with different properties in all dimensions.
            % Not just orientation and drift rate as we are doing here.
            % This is the power of using procedural textures
            phaseLine = rand .* 360;
            propertiesMat = repmat([NaN, freq, sigma, contrast, aspectRatio, 0, 0, 0],...
                nGabors, 1);
            propertiesMat(:, 1) = phaseLine';
            
        elseif ismember(condition, 7:12) % it's part of the velocity of vertical lines experiment
            
            if condition == 7
                cyclesPerSecond = 1;
            elseif condition == 8
                cyclesPerSecond = 5;
            elseif condition == 9
                cyclesPerSecond = 10;
            elseif condition == 10
                cyclesPerSecond = -1;
            elseif condition == 11
                cyclesPerSecond = -5;
            elseif condition == 12
                cyclesPerSecond = -10;
            end
            
            % Translate requested speed of the grating (in cycles per second) into
            % a shift value in "pixels per frame"
            shiftPerFrame = cyclesPerSecond * pixPerCycle * waitDuration;
            
            % Set the frame counter to zero, we need this to 'drift' our grating
            frameCounter = 0;
            
        elseif ismember(condition, 13:18) % is part of the number of dots experimemt
            
            % Dot stuff
            %coh = 0.512;
            coh = 1;
            speed = 3;
            dotSize = 20;
            
            maxDotsPerFrame = 2000; % By trial and error and depends on graphics card, % THIS WILL ALSO BE A KEY PARAMETER FOR THE EXPERIMENT
            
            % Number of dots per video frame = 16.7 dots per sq.deg/sec,
            %dotspersq = 16.7; %per sq.deg/sec, % THIS WILL ALSO BE A KEY PARAMETER FOR THE EXPERIMENT
            if condition == 13
                dotspersq = 10; %per sq.deg/sec, % THIS WILL ALSO BE A KEY PARAMETER FOR THE EXPERIMENT
                direction = 0; % 0 is rightward, 180 is leftward
                this_x = [];
            elseif condition == 14
                dotspersq = 100; %per sq.deg/sec, % THIS WILL ALSO BE A KEY PARAMETER FOR THE EXPERIMENT
                direction = 0; % 0 is rightward, 180 is leftward
                this_x = [];
            elseif condition == 15
                dotspersq = 200; %per sq.deg/sec, % THIS WILL ALSO BE A KEY PARAMETER FOR THE EXPERIMENT
                direction = 0; % 0 is rightward, 180 is leftward
                this_x = [];
            elseif condition == 16
                dotspersq = 10; %per sq.deg/sec, % THIS WILL ALSO BE A KEY PARAMETER FOR THE EXPERIMENT
                direction = 180; % 0 is rightward, 180 is leftward
                this_x = [];
            elseif condition == 17
                dotspersq = 100; %per sq.deg/sec, % THIS WILL ALSO BE A KEY PARAMETER FOR THE EXPERIMENT
                direction = 180; % 0 is rightward, 180 is leftward
                this_x = [];
            elseif condition == 18
                dotspersq = 200; %per sq.deg/sec, % THIS WILL ALSO BE A KEY PARAMETER FOR THE EXPERIMENT
                direction = 180; % 0 is rightward, 180 is leftward
                this_x = [];
            end
            
            % When rounding up, do not exceed the number of dots that can be plotted in
            % a video frame.
            ndots = min(maxDotsPerFrame, ceil(dotspersq * apD .* apD * 0.01 / monRefresh));
            % dxdy is an N x 2 matrix that gives jumpsize in units on 0..1
            %   deg/sec * Ap-unit/deg * sec/jump = unit/jump
            dxdy = repmat(speed * 10/apD * (3/monRefresh) ...
                * [cos(pi*direction/180.0) -sin(pi*direction/180.0)], ndots,1);
            
            % ARRAYS, INDICES for loop
            ss = rand(ndots*3, 2); % array of dot positions raw [xposition, yposition]
            
            % Divide dots into three sets
            Ls = cumsum(ones(ndots,3)) + repmat([0 ndots ndots*2], ndots, 1);
            loopi = 1; % Loops through the three sets of dots
            
        else % is part of the coherence of dots experiment
            
            % Dot stuff
            if condition == 19
                coh = .5;
                direction = 0;
                this_x = [];
            elseif condition == 20
                coh = .75;
                direction = 0;
                this_x = [];
            elseif condition == 21
                coh = 1;
                direction = 0;
                this_x = [];
            elseif condition == 22
                coh = .5;
                direction = 180;
                this_x = [];
            elseif condition == 23
                coh = .75;
                direction = 180;
                this_x = [];
            elseif condition == 24
                coh = 1;
                direction = 180;
                this_x = [];
            end
            
            speed = 3;
            dotSize = 20;
            
            maxDotsPerFrame = 2000; % By trial and error and depends on graphics card, % THIS WILL ALSO BE A KEY PARAMETER FOR THE EXPERIMENT
            
            % Number of dots per video frame = 16.7 dots per sq.deg/sec,
            %dotspersq = 16.7; %per sq.deg/sec, % THIS WILL ALSO BE A KEY PARAMETER FOR THE EXPERIMENT
            dotspersq = 200; %per sq.deg/sec, % THIS WILL ALSO BE A KEY PARAMETER FOR THE EXPERIMENT
            % When rounding up, do not exceed the number of dots that can be plotted in
            % a video frame.
            ndots = min(maxDotsPerFrame, ceil(dotspersq * apD .* apD * 0.01 / monRefresh));
            
            % dxdy is an N x 2 matrix that gives jumpsize in units on 0..1
            %   deg/sec * Ap-unit/deg * sec/jump = unit/jump
            dxdy = repmat(speed * 10/apD * (3/monRefresh) ...
                * [cos(pi*direction/180.0) -sin(pi*direction/180.0)], ndots,1);
            
            % ARRAYS, INDICES for loop
            ss = rand(ndots*3, 2); % array of dot positions raw [xposition, yposition]
            
            % Divide dots into three sets
            Ls = cumsum(ones(ndots,3)) + repmat([0 ndots ndots*2], ndots, 1);
            loopi = 1; % Loops through the three sets of dots
            
            
        end
        
        %% Finally; we start presenting the stimuli and recording position sensor data
        
        % Perform initial flip to gray background and sync us to the retrace:
        vbl = Screen('Flip', window);
        
        % Animation loop
        
        if ismember(condition, 1:6)
            
            tStart = GetSecs;
            % repeat until a valid key is pressed or we time out
            timedout = false;
            while ~timedout
                
                % Set the right blend function for drawing the gabors
                Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');
                
                % Batch draw all of the Gabors to screen
                Screen('DrawTextures', window, gabortex, [], allRects, gaborAngles - 90,...
                    [], [], [], [], kPsychDontDoRotation, propertiesMat');
                
                % Change the blend function to draw an antialiased fixation point
                % in the centre of the array
                Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
                
                % Flip our drawing to the screen
                vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                
                % Increment the phase of our Gabors
                phaseLine = phaseLine + degPerFrameGabors;
                propertiesMat(:, 1) = phaseLine';
                
                tCurrent = GetSecs;
                if tCurrent - tStart > TimePerPresentation
                    timedout = true;
                end
                
                [ x1, y1, x2, y2 ] = ReadArduino_JP_PILOT( mr );
                velocities_stimuli{t, 1} = [velocities_stimuli{t, 1}, x1];
                velocities_stimuli{t, 2} = [velocities_stimuli{t, 2}, y1];
                velocities_stimuli{t, 3} = [velocities_stimuli{t, 3}, x2];
                velocities_stimuli{t, 4} = [velocities_stimuli{t, 4}, y2];
                
                %% Start of water if statement (probably should be made into a function
                if distance < water_threshold_trial % if the animal hasn't run enough
                    
                    master_water_x1(counter_water) = x1 + master_water_x1(counter_water-1); % this is the sum of velocities, which gives you position
                    master_water_y1(counter_water) = y1 + master_water_y1(counter_water-1);
                    master_water_x2(counter_water) = x2 + master_water_x2(counter_water-1);
                    master_water_y2(counter_water) = y2 + master_water_y2(counter_water-1);
                    
                    x = master_water_x1(counter_water) + master_water_x2(counter_water); % this is because of how the sensors work
                    y = master_water_x1(counter_water) - master_water_x2(counter_water);
                    
                    distance = sqrt( (x-0)^2 + (y-0)^2);
                    
                    counter_water = counter_water + 1;
                else % the mouse ran enough, so water is given
                    counter_water = 2;
                    
                    master_water_x1 = [0];
                    master_water_x2 = [0];
                    master_water_y1 = [0];
                    master_water_y2 = [0];
                    
                    distance = 0;
                    
                    fprintf(s, ['1,' water_time]);
                end
                %% End of water if statement
                
            end % end of visual presentation time
            tStart = GetSecs;
            timedout = false;
            while ~timedout
                
                % Set the right blend function for drawing the gabors
                Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');
                
                % draw gray background
                Screen('FillRect', window, [grey, grey, grey], windowRect);
                
                % Change the blend function to draw an antialiased fixation point
                % in the centre of the array
                Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
                
                % Flip our drawing to the screen
                vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                
                tCurrent = GetSecs;
                if tCurrent - tStart > TimePerPresentation
                    timedout = true;
                end
                
                [ x1, y1, x2, y2 ] = ReadArduino_JP_PILOT( mr );
                velocities_no_stimuli{t, 1} = [velocities_no_stimuli{t, 1}, x1];
                velocities_no_stimuli{t, 2} = [velocities_no_stimuli{t, 2}, y1];
                velocities_no_stimuli{t, 3} = [velocities_no_stimuli{t, 3}, x2];
                velocities_no_stimuli{t, 4} = [velocities_no_stimuli{t, 4}, y2];
                
                %% Start of water if statement (probably should be made into a function
                if distance < water_threshold_trial % if the animal hasn't run enough
                    
                    master_water_x1(counter_water) = x1 + master_water_x1(counter_water-1); % this is the sum of velocities, which gives you position
                    master_water_y1(counter_water) = y1 + master_water_y1(counter_water-1);
                    master_water_x2(counter_water) = x2 + master_water_x2(counter_water-1);
                    master_water_y2(counter_water) = y2 + master_water_y2(counter_water-1);
                    
                    x = master_water_x1(counter_water) + master_water_x2(counter_water); % this is because of how the sensors work
                    y = master_water_x1(counter_water) - master_water_x2(counter_water);
                    
                    distance = sqrt( (x-0)^2 + (y-0)^2);
                    
                    counter_water = counter_water + 1;
                else % the mouse ran enough, so water is given
                    counter_water = 2;
                    
                    master_water_x1 = [0];
                    master_water_x2 = [0];
                    master_water_y1 = [0];
                    master_water_y2 = [0];
                    
                    distance = 0;
                    
                    fprintf(s, ['1,' water_time]);
                end
                %% End of water if statement
                
                
                %%
            end
            
        elseif ismember(condition, 7:12)
            
            frameCounter = 0;
            
            tStart = GetSecs;
            timedout = false;
            while ~timedout
                
                % Calculate the xoffset for our window through which to sample our
                % grating
                xoffset = mod(frameCounter * shiftPerFrame, pixPerCycle);
                
                % Now increment the frame counter for the next loop
                frameCounter = frameCounter + 1;
                
                % Define our source rectangle for grating sampling
                srcRect = [xoffset 0 xoffset+visibleSize visibleSize];
                
                % Draw noise texture to the screen
                Screen('DrawTexture', window, noiseTex, [], dstRect, []);
                
                % Draw grating mask
                Screen('DrawTexture', window, gratingMaskTex, srcRect, dstRect, []);
                
                % Flip to the screen on the next vertical retrace
                vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                
                tCurrent = GetSecs;
                if tCurrent - tStart > TimePerPresentation
                    timedout = true;
                end
                
                [ x1, y1, x2, y2 ] = ReadArduino_JP_PILOT( mr );
                velocities_stimuli{t, 1} = [velocities_stimuli{t, 1}, x1];
                velocities_stimuli{t, 2} = [velocities_stimuli{t, 2}, y1];
                velocities_stimuli{t, 3} = [velocities_stimuli{t, 3}, x2];
                velocities_stimuli{t, 4} = [velocities_stimuli{t, 4}, y2];
                
                %% Start of water if statement (probably should be made into a function
                if distance < water_threshold_trial % if the animal hasn't run enough
                    
                    master_water_x1(counter_water) = x1 + master_water_x1(counter_water-1); % this is the sum of velocities, which gives you position
                    master_water_y1(counter_water) = y1 + master_water_y1(counter_water-1);
                    master_water_x2(counter_water) = x2 + master_water_x2(counter_water-1);
                    master_water_y2(counter_water) = y2 + master_water_y2(counter_water-1);
                    
                    x = master_water_x1(counter_water) + master_water_x2(counter_water); % this is because of how the sensors work
                    y = master_water_x1(counter_water) - master_water_x2(counter_water);
                    
                    distance = sqrt( (x-0)^2 + (y-0)^2);
                    
                    counter_water = counter_water + 1;
                else % the mouse ran enough, so water is given
                    counter_water = 2;
                    
                    master_water_x1 = [0];
                    master_water_x2 = [0];
                    master_water_y1 = [0];
                    master_water_y2 = [0];
                    
                    distance = 0;
                    
                    fprintf(s, ['1,' water_time]);
                end
                %% End of water if statement
                
            end
            tStart = GetSecs;
            timedout = false;
            while ~timedout
                
                % Set the right blend function for drawing the gabors
                Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');
                
                % draw gray background
                Screen('FillRect', window, [grey, grey, grey], windowRect);
                
                % Change the blend function to draw an antialiased fixation point
                % in the centre of the array
                Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
                
                % Flip our drawing to the screen
                vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                
                tCurrent = GetSecs;
                if tCurrent - tStart > TimePerPresentation
                    timedout = true;
                end
                
                [ x1, y1, x2, y2 ] = ReadArduino_JP_PILOT( mr );
                velocities_no_stimuli{t, 1} = [velocities_no_stimuli{t, 1}, x1];
                velocities_no_stimuli{t, 2} = [velocities_no_stimuli{t, 2}, y1];
                velocities_no_stimuli{t, 3} = [velocities_no_stimuli{t, 3}, x2];
                velocities_no_stimuli{t, 4} = [velocities_no_stimuli{t, 4}, y2];
                
                
                %% Start of water if statement (probably should be made into a function
                if distance < water_threshold_trial % if the animal hasn't run enough
                    
                    master_water_x1(counter_water) = x1 + master_water_x1(counter_water-1); % this is the sum of velocities, which gives you position
                    master_water_y1(counter_water) = y1 + master_water_y1(counter_water-1);
                    master_water_x2(counter_water) = x2 + master_water_x2(counter_water-1);
                    master_water_y2(counter_water) = y2 + master_water_y2(counter_water-1);
                    
                    x = master_water_x1(counter_water) + master_water_x2(counter_water); % this is because of how the sensors work
                    y = master_water_x1(counter_water) - master_water_x2(counter_water);
                    
                    distance = sqrt( (x-0)^2 + (y-0)^2);
                    
                    counter_water = counter_water + 1;
                else % the mouse ran enough, so water is given
                    counter_water = 2;
                    
                    master_water_x1 = [0];
                    master_water_x2 = [0];
                    master_water_y1 = [0];
                    master_water_y2 = [0];
                    
                    distance = 0;
                    
                    fprintf(s, ['1,' water_time]);
                end
                %% End of water if statement
            end
            
        elseif ismember(condition, 13:24)
            
            frames = 0;
            
            tStart = GetSecs;
            timedout = false;
            while ~timedout
                
                % Set the right blend function for drawing the gabors
                Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');
                
                % Get ss & xs from the big matrices. xs and ss are matrices that have
                % stuff for dots from the last 2 positions + current.
                % Ls picks out the previous set (1:5, 6:10, or 11:15)
                Lthis  = Ls(:,loopi); % Lthis picks out the loop from 3 times ago, which
                % is what is then moved in the current loop
                this_s = ss(Lthis,:);  % this is a matrix of random #s - starting positions
                
                % 1 group of dots are shown in the first frame, a second group are shown
                % in the second frame, a third group shown in the third frame. Then in
                % the next frame, some percentage of the dots from the first frame are
                % replotted according to the speed/direction and coherence, the next
                % frame the same is done for the second group, etc.
                
                % Update the loop pointer
                loopi = loopi+1;
                
                if loopi == 4
                    loopi = 1;
                end
                
                % Compute new locations
                % L are the dots that will be moved
                L = rand(ndots,1) < coh;
                this_s(L,:) = this_s(L,:) + dxdy(L,:);	% Offset the selected dots
                
                if sum(~L) > 0  % if not 100% coherence
                    this_s(~L,:) = rand(sum(~L),2);	% get new random locations for the rest
                end
                
                % Wrap around - check to see if any positions are greater than one or
                % less than zero which is out of the aperture, and then replace with a
                % dot along one of the edges opposite from direction of motion.
                
                N = sum((this_s > 1 | this_s < 0)')' ~= 0;
                if sum(N) > 0
                    xdir = sin(pi*direction/180.0);
                    ydir = cos(pi*direction/180.0);
                    % Flip a weighted coin to see which edge to put the replaced dots
                    if rand < abs(xdir)/(abs(xdir) + abs(ydir))
                        this_s(find(N==1),:) = [rand(sum(N),1) (xdir > 0)*ones(sum(N),1)];
                    else
                        this_s(find(N==1),:) = [(ydir < 0)*ones(sum(N),1) rand(sum(N),1)];
                    end
                end
                
                % Convert to stuff we can actually plot
                this_x(:,1:2) = floor(d_ppd(1) * this_s); % pix/ApUnit
                % This assumes that zero is at the top left, but we want it to be in the
                % center, so shift the dots up and left, which just means adding half of
                % the aperture size to both the x and y direction.
                dot_show = (this_x(:,1:2) - d_ppd/2)';
                
                % After all computations, flip
                % Screen('Flip', window,0,dontclear);
                % vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                % Now do next drawing commands
                Screen('DrawDots', window, dot_show, dotSize, [255 255 255], [xCenter, yCenter]);
                
                % Presentation
                %Screen('DrawingFinished',window,dontclear);
                vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                
                frames = frames + 1;
                
                % Update the arrays so xor works next time
                xs(Lthis, :) = this_x;
                ss(Lthis, :) = this_s;
                
                tCurrent = GetSecs;
                if tCurrent - tStart > TimePerPresentation
                    timedout = true;
                end
                
                [ x1, y1, x2, y2 ] = ReadArduino_JP_PILOT( mr );
                velocities_stimuli{t, 1} = [velocities_stimuli{t, 1}, x1];
                velocities_stimuli{t, 2} = [velocities_stimuli{t, 2}, y1];
                velocities_stimuli{t, 3} = [velocities_stimuli{t, 3}, x2];
                velocities_stimuli{t, 4} = [velocities_stimuli{t, 4}, y2];
                
                
                %% Start of water if statement (probably should be made into a function
                if distance < water_threshold_trial % if the animal hasn't run enough
                    
                    master_water_x1(counter_water) = x1 + master_water_x1(counter_water-1); % this is the sum of velocities, which gives you position
                    master_water_y1(counter_water) = y1 + master_water_y1(counter_water-1);
                    master_water_x2(counter_water) = x2 + master_water_x2(counter_water-1);
                    master_water_y2(counter_water) = y2 + master_water_y2(counter_water-1);
                    
                    x = master_water_x1(counter_water) + master_water_x2(counter_water); % this is because of how the sensors work
                    y = master_water_x1(counter_water) - master_water_x2(counter_water);
                    
                    distance = sqrt( (x-0)^2 + (y-0)^2);
                    
                    counter_water = counter_water + 1;
                else % the mouse ran enough, so water is given
                    counter_water = 2;
                    
                    master_water_x1 = [0];
                    master_water_x2 = [0];
                    master_water_y1 = [0];
                    master_water_y2 = [0];
                    
                    distance = 0;
                    
                    fprintf(s, ['1,' water_time]);
                end
                %% End of water if statement
            end
            
            tStart = GetSecs;
            timedout = false;
            while ~timedout
                
                % Set the right blend function for drawing the gabors
                Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');
                
                % draw gray background
                Screen('FillRect', window, [grey, grey, grey], windowRect);
                
                % Change the blend function to draw an antialiased fixation point
                % in the centre of the array
                Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
                
                % Flip our drawing to the screen
                vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                
                tCurrent = GetSecs;
                if tCurrent - tStart > TimePerPresentation
                    timedout = true;
                end
                
                [ x1, y1, x2, y2 ] = ReadArduino_JP_PILOT( mr );
                velocities_no_stimuli{t, 1} = [velocities_no_stimuli{t, 1}, x1];
                velocities_no_stimuli{t, 2} = [velocities_no_stimuli{t, 2}, y1];
                velocities_no_stimuli{t, 3} = [velocities_no_stimuli{t, 3}, x2];
                velocities_no_stimuli{t, 4} = [velocities_no_stimuli{t, 4}, y2];
                
                
                %% Start of water if statement (probably should be made into a function
                if distance < water_threshold_trial % if the animal hasn't run enough
                    
                    master_water_x1(counter_water) = x1 + master_water_x1(counter_water-1); % this is the sum of velocities, which gives you position
                    master_water_y1(counter_water) = y1 + master_water_y1(counter_water-1);
                    master_water_x2(counter_water) = x2 + master_water_x2(counter_water-1);
                    master_water_y2(counter_water) = y2 + master_water_y2(counter_water-1);
                    
                    x = master_water_x1(counter_water) + master_water_x2(counter_water); % this is because of how the sensors work
                    y = master_water_x1(counter_water) - master_water_x2(counter_water);
                    
                    distance = sqrt( (x-0)^2 + (y-0)^2);
                    
                    counter_water = counter_water + 1;
                else % the mouse ran enough, so water is given
                    counter_water = 2;
                    
                    master_water_x1 = [0];
                    master_water_x2 = [0];
                    master_water_y1 = [0];
                    master_water_y2 = [0];
                    
                    distance = 0;
                    
                    fprintf(s, ['1,' water_time]);
                end
                %% End of water if statement
            end
            
            
        end % end of all different types of condition it can be
        
    end % end of all trials
    
    %% Starting water stuff for post
    master_water_x1 = [0];
    master_water_x2 = [0];
    master_water_y1 = [0];
    master_water_y2 = [0];
    counter_water = 2;
    distance = 0;
    
    tStart = GetSecs;
    timedout = false;
    while ~timedout
        [ x1, y1, x2, y2 ] = ReadArduino_JP_PILOT( mr );
        post{1} = [post{1}, x1];
        post{2} = [post{2}, y1];
        post{3} = [post{3}, x2];
        post{4} = [post{4}, y2];
        
        
        if distance < water_threshold % If the animal wasn't run past a certain distance
            master_water_x1(counter_water) = x1 + master_water_x1(counter_water-1); % we are adding velocities here, which gives us distance
            master_water_y1(counter_water) = y1 + master_water_y1(counter_water-1);
            master_water_x2(counter_water) = x2 + master_water_x2(counter_water-1);
            master_water_y2(counter_water) = y2 + master_water_y2(counter_water-1);
            
            x = master_water_x1(counter_water) + master_water_x2(counter_water); % this is just because of how the sensor works
            y = master_water_x1(counter_water) - master_water_x2(counter_water);
            
            distance = sqrt( (x-0)^2 + (y-0)^2); % this is to keep track not only of forward distance, but both components
            
            counter_water = counter_water + 1;
        else
            counter_water = 2; % re-starting everything for the water
            
            master_water_x1 = [0];
            master_water_x2 = [0];
            master_water_y1 = [0];
            master_water_y2 = [0];
            
            distance = 0;
            
            fprintf(s, ['1,' water_time]); %actually giving the water
            
            
        end
        
        tCurrent = GetSecs;
        if tCurrent - tStart > baseline_time
            timedout = true;
        end
    end
    
    
    sca;
    
    %% save all the data
    save(progFile, 'velocities_stimuli', 'velocities_no_stimuli', 'conditions', 'pre', 'post', 'data')
    fclose(instrfindall)
    
    disp('EVERYTHING WENT WELL!!!')
    
catch
    
    sca;
    commandwindow;
    % re-setting everything to default
    fclose('all');
    Screen('CloseAll');
    ShowCursor;
    ListenChar;
    Priority(0);
    fprintf('Error!!!\n');
    fclose(instrfindall)
end
