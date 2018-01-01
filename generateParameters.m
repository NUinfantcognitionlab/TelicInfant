function [] = generateParameters()
    [calculationsMap, colorsMap, screenInfoMap] = runSetup();

    parametersKeyList = calculationsMap('parametersKeyList');
    translatedParameters = calculationsMap('translatedParameters');

    for p = 1:size(parametersKeyList, 1)
        constantParams = translatedParameters(parametersKeyList(p,1),:);
        alternatingParams = translatedParameters(parametersKeyList(p,2),:);
        % generate constant stimuli
        % check constant stimuli for overlap/off-screen
        % generate alternating stimuli
        % check alternating stimuli for overlap/off-screen
        % save to csv as four rows: constant xpoints, constant ypoints, alternating xpoints, alternating ypoints
    end
    

    sca;
    Priority(0);
end








%Generate a set of x,y points for a single ellipse
function [xpoints, ypoints] = plotEllipse(numberOfFrames, ellipseScale);
    xpoints = [];
    ypoints = [];
    majorAxis = 2*ellipseScale;
    minorAxis = 1*ellipseScale;
    centerX = 0;
    centerY = 0;
    theta = linspace(0,2*pi,numberOfFrames);

    %Start with the basic, unrotated ellipse
    x = (majorAxis/2) * sin(theta) + centerX;
    y = (minorAxis/2) * cos(theta) + centerY;


    %It doesn't start from the right part of the ellipse, so I'm gonna
    %shuffle it around so it does. (this is important I promise)
    %It also adds in some extra frames to smooth the transition between
    %ellipses
    start = round((numberOfFrames)/4);
    x2 = [x(start:numberOfFrames) x(2:start)];
    y2 = [y(start:numberOfFrames) y(2:start)];
    %Finally, accumulate the points in full points arrays for easy graphing
    %and drawing
    xpoints = [xpoints x2];
    ypoints = [ypoints y2];
end

% Rotates one ellipse (with points xpoints and ypoints) a random number of degrees
function [rotatedxpoints rotatedypoints] = rotateEllipse(xpoints, ypoints)
    totalpoints = length(xpoints);


    %In this process, I wind up copying things because I might back up to a
    %different point, and I don't want my calculations to mess with each other.
    %(like, if I change a point, I want the calculations for future points to
    %be calculated from the static previous graph, and not from any changes I
    %just made.

    %So, I have a couple variables that are just copies of the point sets. It's
    %important, I promise.

    %rotate
    rotatedxpoints = xpoints;
    rotatedypoints = ypoints;
    f = randi(360);

    for m = 1:totalpoints-1
        rotatedxpoints(m) = xpoints(m)*cos(f) - ypoints(m)*sin(f);
        rotatedypoints(m) = ypoints(m)*cos(f) + xpoints(m)*sin(f);
    end
end

function [final_xpoints, final_ypoints] = rotateSection(xpoints, ypoints)
    nx = xpoints;
    ny = ypoints;
    xwidth = max(xpoints);
    yheight = max(ypoints);
    totalpoints = length(xpoints);

    petalnum = 0;

    %In this process, I wind up copying things because I might back up to a
    %different point, and I don't want my calculations to mess with each other.
    %(like, if I change a point, I want the calculations for future points to
    %be calculated from the static previous graph, and not from any changes I
    %just made.

    %So, I have a couple variables that are just copies of the point sets. It's
    %important, I promise.

    %Move to origin, since the pieces are a bit uneven
    for m = 1:totalpoints-1
        nx(m) = xpoints(m) - xwidth*.25;
        ny(m) = ypoints(m) - yheight*.25;
    end

    %rotate
    copy_nx = nx;
    copy_ny = ny;
    f = randi(360);

    for m = 1:totalpoints-1
        final_xpoints(m) = nx(m)*cos(f) - ny(m)*sin(f);
        final_ypoints(m) = ny(m)*cos(f) + nx(m)*sin(f);
    end

end

%generates a grid of points to plot the objects on
function [gridCoordinates] = generateGrid(xspaces, yspaces)
    gridCoordinates = [];
    for x = 1:xspaces
        for y = 1:yspaces
            gridCoordinates = [gridCoordinates; [(1+((x-1)*2))/(2*xspaces), (1+((y-1)*2))/(2*yspaces)]];
        end
    end
end

% TODO: remember why I wrote this (may not be necessary with format changes)
function [gridCoordinates] = generateEventsGrid(xspaces, yspaces)
    borderShrink = .05;
    xstart = 1/(2*xspaces) + borderShrink;
    xend = ((2*xspaces)-1)/(2*xspaces) - borderShrink;
    xelements = linspace(xstart, xend, xspaces);
    xcolumn = [];
    for i = xelements
        xcolumn = [xcolumn ; repmat(i, yspaces, 1)];
    end

    ystart = 1/(2*yspaces) + borderShrink;
    yend = ((2*yspaces)-1)/(2*yspaces) - borderShrink;
    yelements = linspace(ystart, yend, yspaces);
    ycolumn = repmat(yelements', xspaces, 1);

    gridCoordinates = [xcolumn ycolumn];
end

% transposes one ellipse to the appropriate position on the grid
function [xpoints ypoints] = transposeEllipse(screenInfoMap, xpoints, ypoints, xsideoffset, gridPosition)
    xpoints = xpoints + xsideoffset + (gridPosition(1)*screenInfoMap('stimXpixels'));
    ypoints = ypoints + (gridPosition(2)*screenInfoMap('screenYpixels'));
end



function [parameters] = readParameters()
  parameters = csvread('TelicInfantParameters.csv',1,0);
  % disp(parameters(1,:))
end

function [imageTexture] = generateImgTexture(imagePath, screenYpixels, window)
    % the alpha value and function here are part of making the transparant image background
    [imagename, ~, alpha] = imread(imagePath);
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
end

%Sets up variables for the experiment
function [calculationsMap, colorsMap, screenInfoMap] = runSetup()
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
    imageTexture = generateImgTexture('star.png', screenYpixels, window);

    % read parameters from file
    parametersKeyList = readParameters();
    % parametersKeyList = parametersKeyList(randperm(size(parametersKeyList,1)),:);

    % Parameters interpretation:
    % 1: 4 ovals, object scale 1
    % 2: 4 ovals, object scale 2
    % 3: 8 ovals, object scale 1
    % 4: 8 ovals, object scale .5
    translatedParameters = [4, 1;
                   4, 2;
                   8, 1;
                   8,.5];

    % I'm making a bunch of containers that match variable names to values.
    % This way, I can just pass the container to a function and have all the myriad
    % info I need to write to the screen, without writing and re-writing all
    % the variables in order to the argument list.
    calculationsMap = containers.Map({'framesPerObjectLoop', 'framesPerLoop', 'minSpace', 'breakTime', ...
    'displayTime', 'blankscreenTime', 'scale', 'xGridSpaces', 'yGridSpaces', 'translatedParameters', 'parametersKeyList'}, {framesPerObjectLoop, framesPerLoop, minSpace, ...
    breakTime, displayTime, blankscreenTime, scale, xGridSpaces, yGridSpaces, translatedParameters, parametersKeyList});
    colorsMap = containers.Map({'screenBlack', 'screenWhite', 'screenGrey', ...
    'rgbgrey'}, {black, white, grey, rgbgrey});
    screenInfoMap = containers.Map({'window', 'vbl', 'ifi', 'baseRect', ...
    'screenXpixels', 'screenYpixels', 'stimXpixels', 'xCenter', 'yCenter', 'leftxCenter', 'rightxCenter', 'screenNumber', 'imageTexture'}, {window, vbl, ifi, baseRect, screenXpixels, (screenYpixels/3)*3, stimXpixels, xCenter, yCenter, leftxCenter, rightxCenter, screenNumber, imageTexture});
end