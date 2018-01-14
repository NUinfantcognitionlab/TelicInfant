function [] = generateParameters()
    [calculationsMap, colorsMap, screenInfoMap] = runSetup();

    parametersKeyList = calculationsMap('parametersKeyList');
    translatedParameters = calculationsMap('translatedParameters');

    constNat = zeros(150, 300);
    constUnnat = zeros(150, 300);
    alternateNat = zeros(150, 300);
    alternateUnnat = zeros(150, 300);

    constNatBreaks = zeros(75, 8);
    constUnnatBreaks = zeros(75, 8);
    alternateNatBreaks = zeros(75, 8);
    alternateUnnatBreaks = zeros(75, 8);

    for p = 1:size(parametersKeyList, 1)
        constantParams = translatedParameters(parametersKeyList(p,1),:);
        alternatingParams = translatedParameters(parametersKeyList(p,2),:);

        scale = calculationsMap('scale');

        gridCoordinates = generateGrid(calculationsMap('xGridSpaces'), calculationsMap('yGridSpaces'));
        twoscalegridCoordinates = generateGrid(2,4);


        %CONST NAT
        generationInvalid = true;
        %while there isn't a valid generated set of points
        while generationInvalid
            gridPositions = [1:(calculationsMap('xGridSpaces')*calculationsMap('yGridSpaces'))];
            gridPositions = gridPositions(randperm(length(gridPositions)));

            twoscalegridPositions = [1:(2*4)];
            twoscalegridPositions = twoscalegridPositions(randperm(length(twoscalegridPositions)));
            %generate new set of points
            [xpoints ypoints breakList] = generateNaturalCoordinateSet(calculationsMap, screenInfoMap, ...
            constantParams(1), constantParams(2), 0, scale, gridCoordinates, gridPositions, twoscalegridCoordinates, twoscalegridPositions, calculationsMap('framesPerLoop'));            
            
            generationInvalid = false;
            previousBreak = 1;
            nextBreak = breakList(1);
            breakIndex = 2;
            
            %for every point in xpoints
            for q = 1:size(xpoints, 2)
                if xpoints(q) <40 || xpoints(q)>screenInfoMap('stimXpixels')-40 || ypoints(q)<40 || ypoints(q)>screenInfoMap('screenYpixels')-40
                    disp('off screen!!')
                    generationInvalid = true;
                    break;
                end

                if q == nextBreak    
                    if breakIndex <= numel(breakList)
                        previousBreak = nextBreak;
                        nextBreak = breakList(breakIndex);
                        breakIndex = breakIndex + 1;
                        continue;
                    end
                end
                % if a point is the next break, shift the sections to examine (the +/- 1s are to exclude checking breaks with each other)
                pointRange = [1:previousBreak-1,nextBreak:size(xpoints, 2)-1];
                if previousBreak == 1
                    pointRange = [nextBreak:size(xpoints, 2)-1];
                end
                for r = pointRange
                    if pdist([xpoints(q),ypoints(q);xpoints(r),ypoints(r)]) < 3;
                        disp('overlap!!')
                        generationInvalid = true;
                        break;
                    end
                end

                if generationInvalid
                    break;
                end
            end
            
        end
        % plot(xpoints, ypoints);
        constNat(p*2, 1:numel(xpoints)) = xpoints;
        constNat(p*2+1, 1:numel(ypoints)) = ypoints;
        % constNatBreaks = [constNatBreaks;breakList];
        constNatBreaks(p, 1:numel(breakList)) = breakList;
        % break;

        %%%%%%%%%%%%%%%%%%CONST UNNAT
        generationInvalid = true;
        while generationInvalid
            gridPositions = [1:(calculationsMap('xGridSpaces')*calculationsMap('yGridSpaces'))];
            gridPositions = gridPositions(randperm(length(gridPositions)));

            twoscalegridPositions = [1:(2*4)];
            twoscalegridPositions = twoscalegridPositions(randperm(length(twoscalegridPositions)));
            %generate new set of points
            [xpoints ypoints breakList] = generateRandomCoordinateSet(calculationsMap, screenInfoMap, ...
            constantParams(1), constantParams(2), 0, scale, gridCoordinates, gridPositions, twoscalegridCoordinates, twoscalegridPositions, calculationsMap('framesPerLoop'));
             
            generationInvalid = false;
            previousBreak = 1;
            nextBreak = breakList(1);
            breakIndex = 2;
            
            %for every point in xpoints
            for q = 1:size(xpoints, 2)
                if xpoints(q) <40 || xpoints(q)>screenInfoMap('stimXpixels')-40 || ypoints(q)<40 || ypoints(q)>screenInfoMap('screenYpixels')-40
                    disp('off screen!!')
                    generationInvalid = true;
                    break;
                end

                if q == nextBreak    
                    if breakIndex <= numel(breakList)
                        previousBreak = nextBreak;
                        nextBreak = breakList(breakIndex);
                        breakIndex = breakIndex + 1;
                        continue;
                    end
                end
                % if a point is the next break, shift the sections to examine (the +/- 1s are to exclude checking breaks with each other)
                pointRange = [1:previousBreak-1,nextBreak:size(xpoints, 2)-1];
                if previousBreak == 1
                    pointRange = [nextBreak:size(xpoints, 2)-1];
                end
                for r = pointRange
                    if pdist([xpoints(q),ypoints(q);xpoints(r),ypoints(r)]) < 3;
                        disp('overlap!!')
                        generationInvalid = true;
                        break;
                    end
                end

                if generationInvalid
                    break;
                end
            end
            
        end
        % plot(xpoints, ypoints);
        constUnnat(p*2, 1:numel(xpoints)) = xpoints;
        constUnnat(p*2+1, 1:numel(ypoints)) = ypoints;
        constUnnatBreaks(p, 1:numel(breakList)) = breakList;
        % break;



        %%%%%%%%%%%%%%%%%%ALTERNATE UNNAT
        generationInvalid = true;
        while generationInvalid
            gridPositions = [1:(calculationsMap('xGridSpaces')*calculationsMap('yGridSpaces'))];
            gridPositions = gridPositions(randperm(length(gridPositions)));

            twoscalegridPositions = [1:(2*4)];
            twoscalegridPositions = twoscalegridPositions(randperm(length(twoscalegridPositions)));
            [xpoints ypoints breakList] = generateNaturalCoordinateSet(calculationsMap, screenInfoMap, ...
            alternatingParams(1), alternatingParams(2), 0, scale, gridCoordinates, gridPositions, twoscalegridCoordinates, twoscalegridPositions, calculationsMap('framesPerLoop'));
            
            generationInvalid = false;
            previousBreak = 1;
            nextBreak = breakList(1);
            breakIndex = 2;
            
            %for every point in xpoints
            for q = 1:size(xpoints, 2)
                if xpoints(q) <40 || xpoints(q)>screenInfoMap('stimXpixels')-40 || ypoints(q)<40 || ypoints(q)>screenInfoMap('screenYpixels')-40
                    disp('off screen!!')
                    generationInvalid = true;
                    break;
                end

                if q == nextBreak    
                    if breakIndex <= numel(breakList)
                        previousBreak = nextBreak;
                        nextBreak = breakList(breakIndex);
                        breakIndex = breakIndex + 1;
                        continue;
                    end
                end
                % if a point is the next break, shift the sections to examine (the +/- 1s are to exclude checking breaks with each other)
                pointRange = [1:previousBreak-1,nextBreak:size(xpoints, 2)-1];
                if previousBreak == 1
                    pointRange = [nextBreak:size(xpoints, 2)-1];
                end
                for r = pointRange
                    if pdist([xpoints(q),ypoints(q);xpoints(r),ypoints(r)]) < 3;
                        disp('overlap!!')
                        generationInvalid = true;
                        break;
                    end
                end

                if generationInvalid
                    break;
                end
            end
            
        end

        alternateNat(p*2, 1:numel(xpoints)) = xpoints(1,:);
        alternateNat(p*2+1, 1:numel(ypoints)) = ypoints(1,:);
        alternateNatBreaks(p, 1:numel(breakList)) = breakList;
        % break;



        %%%%%%%%%%%%%%%%%%ALTERNATE UNNAT
        generationInvalid = true;
        while generationInvalid
            gridPositions = [1:(calculationsMap('xGridSpaces')*calculationsMap('yGridSpaces'))];
            gridPositions = gridPositions(randperm(length(gridPositions)));

            twoscalegridPositions = [1:(2*4)];
            twoscalegridPositions = twoscalegridPositions(randperm(length(twoscalegridPositions)));
            [xpoints ypoints breakList] = generateRandomCoordinateSet(calculationsMap, screenInfoMap, ...
            alternatingParams(1), alternatingParams(2), 0, scale, gridCoordinates, gridPositions, twoscalegridCoordinates, twoscalegridPositions, calculationsMap('framesPerLoop'));
            
            
            generationInvalid = false;
            previousBreak = 1;
            nextBreak = breakList(1);
            breakIndex = 2;
            
            %for every point in xpoints
            for q = 1:size(xpoints, 2)
                if xpoints(q) <40 || xpoints(q)>screenInfoMap('stimXpixels')-40 || ypoints(q)<40 || ypoints(q)>screenInfoMap('screenYpixels')-40
                    disp('off screen!!')
                    generationInvalid = true;
                    break;
                end

                if q == nextBreak    
                    if breakIndex <= numel(breakList)
                        previousBreak = nextBreak;
                        nextBreak = breakList(breakIndex);
                        breakIndex = breakIndex + 1;
                        continue;
                    end
                end
                % if a point is the next break, shift the sections to examine (the +/- 1s are to exclude checking breaks with each other)
                pointRange = [1:previousBreak-1,nextBreak:size(xpoints, 2)-1];
                if previousBreak == 1
                    pointRange = [nextBreak:size(xpoints, 2)-1];
                end
                for r = pointRange
                    if pdist([xpoints(q),ypoints(q);xpoints(r),ypoints(r)]) < 3;
                        disp('overlap!!')
                        generationInvalid = true;
                        break;
                    end
                end

                if generationInvalid
                    break;
                end
            end
            
        end
        alternateUnnat(p*2, 1:numel(xpoints)) = xpoints;
        alternateUnnat(p*2+1, 1:numel(ypoints)) = ypoints;
        alternateUnnatBreaks(p, 1:numel(breakList)) = breakList;
        % break;
        
    end

    csvwrite('ConstNat.csv',constNat(2:size(constNat, 1),:));
    csvwrite('ConstUnnat.csv',constUnnat(2:size(constUnnat, 1),:));
    csvwrite('AlternateNat.csv',alternateNat(2:size(alternateNat, 1),:));
    csvwrite('AlternateUnnat.csv',alternateUnnat(2:size(alternateUnnat, 1),:));

    csvwrite('ConstNatBreaks.csv',constNatBreaks);
    csvwrite('ConstUnnatBreaks.csv',constUnnatBreaks);
    csvwrite('AlternateNatBreaks.csv',alternateNatBreaks);
    csvwrite('AlternateUnnatBreaks.csv',alternateUnnatBreaks);
    

    sca;
    Priority(0);
end




% Generates the x-y coordinate pairs for the ellipses on the grid, but does not draw them.
function [xpoints ypoints breakList] = generateNaturalCoordinateSet(calculationsMap, screenInfoMap, ...
  numberOfLoops, ellipseScale, xsideoffset, scale, gridCoordinates, gridPositions, twoscalegridCoordinates, twoscalegridPositions, framesPerLoop)
    xpoints = [];
    ypoints = [];
    breakList = [];
    if numberOfLoops == 8
        framesPerLoop = round(calculationsMap('framesPer8LoopAnimation')/numberOfLoops); 
    else
        framesPerLoop = round(calculationsMap('framesPer4LoopAnimation')/numberOfLoops);
    end
    for i = 1:numberOfLoops
        [xEllipse, yEllipse] = plotEllipse(framesPerLoop, ellipseScale);
        [xEllipse, yEllipse] = rotateEllipse(xEllipse, yEllipse);
        xEllipse = xEllipse .* scale;
        yEllipse = yEllipse .* scale;

        % if the scale is big, use the smaller grid. otherwise, use the normal grid
        if(ellipseScale==2)
            [xEllipse, yEllipse] = transposeEllipse(screenInfoMap, xEllipse, yEllipse, twoscalegridCoordinates(twoscalegridPositions(i),:));
        else
            [xEllipse, yEllipse] = transposeEllipse(screenInfoMap, xEllipse, yEllipse, gridCoordinates(gridPositions(i),:));
        end

        xpoints = [xpoints xEllipse];
        ypoints = [ypoints yEllipse];
        breakList = [breakList numel(xpoints)];
    end
end

function [processedxPoints processedyPoints breakList] = generateRandomCoordinateSet(calculationsMap, screenInfoMap, ...
  numberOfLoops, ellipseScale, xsideoffset, scale, gridCoordinates, gridPositions, twoscalegridCoordinates, twoscalegridPositions, framesPerLoop)
    if numberOfLoops == 8
        framesPerLoop = round(calculationsMap('framesPer8LoopAnimation')/numberOfLoops); 
    else
        framesPerLoop = round(calculationsMap('framesPer4LoopAnimation')/numberOfLoops);
    end
    
    [xpoints, ypoints] = getEllipseSetPoints(numberOfLoops, framesPerLoop, ellipseScale);
    % sort the breaks by ascending to deal with each chunk of points separately
    breakList = sort(generateBreakList('random', numel(xpoints), numberOfLoops, calculationsMap('minSpace')), 'ascend');
    % placeholder for breaking up the points and index for progressing through grid positions
    tempi = 1;
    gridIndex = 1;
    processedxPoints = [];
    processedyPoints = [];
    for i = breakList(1:end)
        % get the points from the previous index to the next
        currentxPoints = xpoints(tempi:i);
        currentyPoints = ypoints(tempi:i);
        % rotate and transpose
        [currentxPoints, currentyPoints] = rotateSection(currentxPoints, currentyPoints);
        currentxPoints = currentxPoints .* scale;
        currentyPoints = currentyPoints .* scale;

        % if the scale is big, use the smaller grid. otherwise, use the normal grid
        if(ellipseScale==2)
            [currentxPoints, currentyPoints] = transposeEllipse(screenInfoMap, currentxPoints, currentyPoints, twoscalegridCoordinates(twoscalegridPositions(gridIndex),:));
        else
            [currentxPoints, currentyPoints] = transposeEllipse(screenInfoMap, currentxPoints, currentyPoints, gridCoordinates(gridPositions(gridIndex),:));
        end

        processedxPoints = [processedxPoints currentxPoints];
        processedyPoints = [processedyPoints currentyPoints];

        % reset iterators
        tempi = i;
        gridIndex = gridIndex + 1;
    end
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
function [xpoints ypoints] = transposeEllipse(screenInfoMap, xpoints, ypoints, gridPosition)
    xpoints = xpoints + (gridPosition(1)*screenInfoMap('stimXpixels'));
    ypoints = ypoints + (gridPosition(2)*screenInfoMap('screenYpixels'));
end

function [parameters] = readParameters()
  parameters = csvread('TelicInfantParameters.csv',1,0);
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

    minSpace = 10;
    %the minimum possible number of frames between steps
    breakTime = .25;

    %will ideally be used to generate loopTime and framesPerLoop on the fly
    totalTimePerAnimation = 4.5;
    breakFrames = round(breakTime / ifi);
    framesPer4LoopAnimation = round(totalTimePerAnimation / ifi) - breakFrames*4;
    framesPer8LoopAnimation = round(totalTimePerAnimation / ifi) - breakFrames*8;

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
    calculationsMap = containers.Map({'totalTimePerAnimation','framesPer4LoopAnimation', 'framesPer8LoopAnimation', 'framesPerLoop', 'minSpace', 'breakTime', ...
    'displayTime', 'blankscreenTime', 'scale', 'xGridSpaces', 'yGridSpaces', 'translatedParameters', 'parametersKeyList'}, {totalTimePerAnimation, framesPer4LoopAnimation, framesPer8LoopAnimation, framesPerLoop, minSpace, ...
    breakTime, displayTime, blankscreenTime, scale, xGridSpaces, yGridSpaces, translatedParameters, parametersKeyList});
    colorsMap = containers.Map({'screenBlack', 'screenWhite', 'screenGrey', ...
    'rgbgrey'}, {black, white, grey, rgbgrey});
    screenInfoMap = containers.Map({'window', 'vbl', 'ifi', 'baseRect', ...
    'screenXpixels', 'screenYpixels', 'stimXpixels', 'xCenter', 'yCenter', 'leftxCenter', 'rightxCenter', 'screenNumber', 'imageTexture'}, {window, vbl, ifi, baseRect, screenXpixels, (screenYpixels/3)*3, stimXpixels, xCenter, yCenter, leftxCenter, rightxCenter, screenNumber, imageTexture});
    sca;
end

function [xpoints, ypoints] = getEllipseSetPoints(numberOfLoops, framesPerLoop, ellipseScale)
    %OK, so, the ellipses weren't lining up at the origin very well, so
    %smoothframes designates a few frames to smooth this out. It uses fewer
    %frames for the ellipse, and instead spends a few frames going from the
    %end of the ellipse to the origin.
    smoothframes = 0;
    doublesmooth = smoothframes*2;
    xpoints = [];
    ypoints = [];
    majorAxis = 2*ellipseScale;
    minorAxis = 1*ellipseScale;
    centerX = 0;
    centerY = 0;
    theta = linspace(0,2*pi,framesPerLoop-smoothframes);
    %The orientation starts at 0, and ends at 360-360/numberOfLoops
    %This is to it doesn't make a complete circle, which would have two
    %overlapping ellipses.
    orientation = linspace(0,360-round(360/numberOfLoops),numberOfLoops);
    for i = 1:numberOfLoops
        %orientation calculated from above
        loopOri=orientation(i)*pi/180;

        %Start with the basic, unrotated ellipse
        initx = (majorAxis/2) * sin(theta) + centerX;
        inity = (minorAxis/2) * cos(theta) + centerY;

        %Then rotate it
        x = (initx-centerX)*cos(loopOri) - (inity-centerY)*sin(loopOri) + centerX;
        y = (initx-centerX)*sin(loopOri) + (inity-centerY)*cos(loopOri) + centerY;
        %then push it out based on the rotation
        for m = 1:numel(x)
            x2(m) = x(m) + (x(round(numel(x)*.75)) *1);
            y2(m) = y(m) + (y(round(numel(y)*.75)) *1);
        end

        %It doesn't start from the right part of the ellipse, so I'm gonna
        %shuffle it around so it does. (this is important I promise)
        %It also adds in some extra frames to smooth the transition between
        %ellipses
        start = round((framesPerLoop-smoothframes)/4);
        x3 = [x2(start:framesPerLoop-smoothframes) x2(2:start) linspace(x2(start),0,smoothframes)];
        y3 = [y2(start:framesPerLoop-smoothframes) y2(2:start) linspace(y2(start),0,smoothframes)];
        %Finally, accumulate the points in full points arrays for easy graphing
        %and drawing
        xpoints = [xpoints x3];
        ypoints = [ypoints y3];
    end
end

% generates a list of break points to separate ellipses/pieces naturally or randomly
function [Breaks] = generateBreakList(breakType, totalpoints, loops, minSpace)
    if strcmp(breakType, 'equal')
        %Breaks = 1 : totalpoints/loops : totalpoints;
        Breaks = linspace(totalpoints/loops, totalpoints, loops);
        Breaks = arrayfun(@(x) round(x),Breaks);

    elseif strcmp(breakType, 'random')
        %tbh I found this on stackoverflow and have no idea how it works
        %http://stackoverflow.com/questions/31971344/generating-random-sequence-with-minimum-distance-between-elements-matlab/31977095#31977095
        if loops >1
            numberOfBreaks = loops - 1;
            %The -10 accounts for some distance away from the last point,
            %which I add on separately.
            E = (totalpoints-10)-(numberOfBreaks-1)*minSpace;

            ro = rand(numberOfBreaks+1,1);
            rn = E*ro(1:numberOfBreaks)/sum(ro);

            s = minSpace*ones(numberOfBreaks,1)+rn;

            Breaks=cumsum(s)-1;

            Breaks = reshape(Breaks, 1, length(Breaks));
            Breaks = arrayfun(@(x) round(x),Breaks);

            %I'm adding one break on at the end, otherwise I'll end up with
            %more "pieces" than in the equal condition.
            Breaks = [Breaks totalpoints];
        else
            Breaks = [totalpoints];
        end
        

    else
        Breaks = [];
    end
end