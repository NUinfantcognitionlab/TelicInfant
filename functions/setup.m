function [calculationsMap, colorsMap, screenInfoMap] = setup()
    %The following lines set up the Psychtoolbox environment.
    Screen('Preference', 'SkipSyncTests', 1);
    %The previous line sets the experiment to run screen sychronization tests,
    %to make timing accurate. If the experiment won't run, change the 0 to a 1
    %to skip those tests, at the risk of innacurate timing. (Make sure you try
    %changing the screen resolution to another scale and back! That might fix
    %things without skipping sync tests.)
    close all;
    sca
    PsychDefaultSetup(2);
    screens = Screen('Screens');
    % in the event that you need to run on a separate monitor,
    % screenNumber would need to change
    screenNumber = max(screens);
    rng('shuffle');
    KbName('UnifyKeyNames');

    % cond = input('Condition r (right alternating) or l (left alternating): ', 's');
    % alternatingSide = condcheck(cond);
    alternatingSide = 'left';
    %Define Colors
    white = WhiteIndex(screenNumber);
    black = BlackIndex(screenNumber);
    grey = white/2;
    rgbgrey = [.5 .5 .5];
    %%%Screen Stuff
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
    %opens a window in the most external screen and colors it)
    Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    %Anti-aliasing or something? It's from a tutorial
    ifi = Screen('GetFlipInterval', window);
    %Drawing intervals; used to change the screen to animate the image
    [screenXpixels, screenYpixels] = Screen('WindowSize', window);
    %The size of the screen window in pixels
    [xCenter, yCenter] = RectCenter(windowRect);
    %The center of the screen window
    % Each screen is the same size, calculated here
    stimXpixels = (screenXpixels/30)*10;
    % The x axis center is different for each screen (calculated here), but the
    % y center (above, line 40) is the same
    leftxCenter = (screenXpixels/30)*5;
    rightxCenter = (screenXpixels/30)*25;
    % ListenChar(0);
    HideCursor(screenNumber);

    %all times in seconds
    LoopTime = .75;
    framesPerLoop = round(LoopTime / ifi) + 1;
    framesPerObjectLoop = 80;
    minSpace = 10;
    %the minimum possible number of frames between steps
    breakTime = .25;
    %The number of seconds for each pause
    displayTime = .5;
    %number of seconds for which to display images
    blankscreenTime = .3;
    %Length of space between loops presentation
    % textsize = 40;
    % textspace = 1.5;
    xGridSpaces = 4;
    yGridSpaces = 5;
    scale = screenYpixels / 37;%previously 15
    %Matlab's strings are stupid, so I have quotes and quotes with spaces in
    %variables here
    quote = '''';
    squote = ' ''';
    vbl = Screen('Flip', window);
    Priority(MaxPriority(window));
    % change this to shrink the grey background
    baseRect = [0 0 stimXpixels screenYpixels];

    % the alpha value and function here are part of making the transparant image background
    [imagename, ~, alpha] = imread('star.png');
    imagename(:,:,4) = alpha(:,:);

    % Get the size of the image
    [s1, s2, ~] = size(imagename);

    % Here we check if the image is too big to fit on the screen and abort if
    % it is. See ImageRescaleDemo to see how to rescale an image.
    if s1 > screenYpixels || s2 > screenYpixels
        disp('ERROR! Image is too big to fit on the screen');
        sca;
        return;
    end

    % Make the image into a texture
    imageTexture = Screen('MakeTexture', window, imagename);

    % read parameters from file
    parametersKeyList = csvread('TelicInfantParameters.csv',1,0);;
    % parametersKeyList = parametersKeyList(randperm(size(parametersKeyList,1)),:);

    % Parameters interpretation:
    % 1: 4 ovals, object scale 1
    % 2: 4 ovals, object scale 2
    % 3: 8 ovals, object scale 1
    % 4: 8 ovals, object scale .5
    paramValues = [4, 1;
                   4, 2;
                   8, 1;
                   8,.5];

    % I'm making a bunch of containers that match variable names to values.
    % This way, I can just pass the container to a function and have all the myriad
    % info I need to write to the screen, without writing and re-writing all
    % the variables in order to the argument list.
    calculationsMap = containers.Map({'framesPerObjectLoop', 'framesPerLoop', 'minSpace', 'breakTime', ...
    'displayTime', 'blankscreenTime', 'scale', 'xGridSpaces', 'yGridSpaces', 'paramValues', 'parametersKeyList'}, {framesPerObjectLoop, framesPerLoop, minSpace, ...
    breakTime, displayTime, blankscreenTime, scale, xGridSpaces, yGridSpaces, paramValues, parametersKeyList});
    colorsMap = containers.Map({'screenBlack', 'screenWhite', 'screenGrey', ...
    'rgbgrey'}, {black, white, grey, rgbgrey});
    screenInfoMap = containers.Map({'window', 'vbl', 'ifi', 'baseRect', ...
    'screenXpixels', 'screenYpixels', 'stimXpixels', 'xCenter', 'yCenter', 'leftxCenter', 'rightxCenter', 'screenNumber', 'imageTexture'}, {window, vbl, ifi, baseRect, screenXpixels, (screenYpixels/3)*3, stimXpixels, xCenter, yCenter, leftxCenter, rightxCenter, screenNumber, imageTexture});
end