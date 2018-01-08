function [] = TelicInfant()
    % Set up initial conditions via the command line
    condition = input('Condition e (event) or o (object): ', 's');
    displayType = displayTypeInputcheck(condition);
    condition = input('Condition r (right alternating) or l (left alternating): ', 's');
    alternatingSide = alternatingInputcheck(condition);
    condition = input('Condition nat (natural) or unnat (unnatural): ', 's');
    breakType = breakTypeInputcheck(condition);

    %runs a bunch of Psychtoolbox setup and variable calculation, and returns some hashmaps storing that information for later retrieval
    % calculationsMap : framesPerObjectLoop, framesPerLoop, minSpace, breakTime, displayTime, blankscreenTime, scale, xGridSpaces, yGridSpaces, translatedParameters, parametersKeyList
    % colorsMap : screenBlack, screenWhite, screenGrey, rgbgrey
    % screenInfoMap : window, vbl, ifi, baseRect, screenXpixels, screenYpixels, stimXpixels, xCenter, yCenter, leftxCenter, rightxCenter, screenNumber, imageTexture
    [calculationsMap, colorsMap, screenInfoMap] = runSetup();
    timePerAnimation = 4.5;
    if strcmp(breakType, 'equal')
        fileTypeString = 'Nat';
    else
        fileTypeString = 'Unnat';
    end
    alternateParametersList = csvread(['Alternate', fileTypeString, '.csv'], 0, 0);
    alternateBreakList = csvread(['Alternate', fileTypeString, 'Breaks.csv'], 0, 0);
    constParametersList = csvread(['Const', fileTypeString, '.csv'], 0, 0);
    constBreakList = csvread(['Const', fileTypeString, 'Breaks.csv'], 0, 0);

    

    % NOTE: The projector mirrors the view, so 'left' here is used to indicate the right side of a non-projected screen.
    % alternatingSide = 'right';
    % breakType = 'random';

    for t = 1:1
        attentionScreen(screenInfoMap, colorsMap);
        if strcmp(displayType, 'object')
            runObjectTrial(calculationsMap, screenInfoMap, colorsMap, alternatingSide, breakType, alternateParametersList, alternateBreakList, constParametersList, constBreakList);
        else
            runEventsFromPoints(calculationsMap, screenInfoMap, colorsMap, alternatingSide, breakType, alternateParametersList, alternateBreakList, constParametersList, constBreakList);
        end
        if strcmp(alternatingSide, 'left')
            alternatingSide = 'right';
        else
            alternatingSide = 'left';
        end
    end
    endingScreen(screenInfoMap, colorsMap);

    sca
    Priority(0);
end


%%%%%%%%%
% FUNCTIONS FOR OBJECTS

% Manages running the Object condition for an amount of time. Uses showObjectStimuli to display stims
function [] = runObjectTrial(calculationsMap, screenInfoMap, colorsMap, alternatingSide, breakType,  alternateParametersList, alternateBreakList, constParametersList, constBreakList)
    parametersKeyList = calculationsMap('parametersKeyList');
    p = 1;
    translatedParameters = calculationsMap('translatedParameters');
    % this boolean is calculated up here to make sure the conditional during stim presentation is as fast as possible
    leftAlternating = strcmp(alternatingSide, 'left');
    finalTime = datenum(clock + [0, 0, 0, 0, 0, 10]);
    
    if strcmp(breakType, 'equal')
        generationFunction = @drawObjectsFromPoints;
    else
        generationFunction = @drawObjectsFromPoints;
    end

    while datenum(clock) < finalTime
        % read parameters from the list based on the current trial number
        % constantParams contains both x and ypoints
        constantParams = [constParametersList(p,:); constParametersList(p+1,:)];
        alternatingParams = [alternateParametersList(p,:); alternateParametersList(p+1,:)];
        constantBreakList = constBreakList(((p+1)/2),:)
        alternatingBreakList = alternateBreakList(((p+1)/2),:)
        if leftAlternating
            % because of the projector mirroring, putting the alternation on the right of the screen will be on the left of the projector screen
            showObjectStimuli(calculationsMap, screenInfoMap, colorsMap, breakType, constantParams(1,:), constantParams(2,:), constantBreakList, alternatingParams(1,:), alternatingParams(2,:), alternatingBreakList);
        else
            showObjectStimuli(calculationsMap, screenInfoMap, colorsMap, breakType, alternatingParams(1,:), alternatingParams(2,:), alternatingBreakList, constantParams(1,:), constantParams(2,:), constantBreakList);
        end
        p = p+2;
        if p > numel(constParametersList)/2
            p = 0;
        end
    end
end

% takes parameters for two sides of stimuli and draws both to the screen for the amount of time specified for the display
% Uses generateObjectSet to generate the points for the object and draw them to the screen
function [] = showObjectStimuli(calculationsMap, screenInfoMap, colorsMap, breakType, leftxpoints, leftypoints, leftBreaks, rightxpoints, rightypoints, rightBreaks)
    window = screenInfoMap('window');
    black = colorsMap('screenBlack');
    % plot(leftxpoints, leftypoints)

    disp(leftxpoints(:, 1:5))
    disp(rightxpoints(:, 1:5))
    
    drawObjectsFromPoints(calculationsMap, screenInfoMap, colorsMap, breakType, ...
        'left', leftxpoints, leftypoints, leftBreaks);
    drawObjectsFromPoints(calculationsMap, screenInfoMap, colorsMap, breakType, ...
        'right', rightxpoints, rightypoints, rightBreaks);
    screenInfoMap('vbl') = Screen('Flip', window, screenInfoMap('vbl') + calculationsMap('blankscreenTime'));
    drawBlankScreen(screenInfoMap, colorsMap);
    screenInfoMap('vbl') = Screen('Flip', window, screenInfoMap('vbl')+calculationsMap('displayTime'));
end


function [] = drawObjectsFromPoints(calculationsMap, screenInfoMap, colorsMap, breakType, ...
  screenside, xpoints, ypoints, breakList)
    scale = calculationsMap('scale');
    framesPerLoop = calculationsMap('framesPerObjectLoop');
    window = screenInfoMap('window');
    black = colorsMap('screenBlack');
    xGridSpaces = calculationsMap('xGridSpaces');
    yGridSpaces = calculationsMap('yGridSpaces');

    yCenter = screenInfoMap('yCenter');
    if strcmp(screenside, 'left')
        xCenter = screenInfoMap('leftxCenter');
    else
        xCenter = screenInfoMap('rightxCenter');
        xpoints = xpoints + screenInfoMap('stimXpixels')*2;
    end

    totalpoints = numel(xpoints);
    bgRect = CenterRectOnPointd(screenInfoMap('baseRect'), xCenter, yCenter);
    Screen('FillRect', window, colorsMap('rgbgrey'), bgRect);
    savepoint = 1;
    for p = 1:totalpoints - 2
        if ~any(p == breakList) && ~any(p+1 == breakList)
            Screen('DrawLine', window, black, xpoints(p), ypoints(p), ...
                xpoints(p+1), ypoints(p+1), 5);
        elseif strcmp(breakType, 'equal')
            Screen('DrawLine', window, black, xpoints(p), ypoints(p), ...
                xpoints(savepoint), ypoints(savepoint), 5);
            savepoint = p+1;
        end
    end
end

% Draws sets of object stimuli to the screen, but DOESN'T flip the screen to make them visible.
% The side of the screen, and the size of the loops are specified
% plotEllipse, rotateEllipse, generateGrid, and transposeEllipse are used to generate, rotate, and place ellipses
function [] = drawNaturalObjectSet(calculationsMap, screenInfoMap, colorsMap, ...
  numberOfLoops, ellipseScale, screenside, lineColor)
    scale = calculationsMap('scale');
    framesPerLoop = calculationsMap('framesPerObjectLoop');
    window = screenInfoMap('window');
    black = colorsMap('screenBlack');
    xGridSpaces = calculationsMap('xGridSpaces');
    yGridSpaces = calculationsMap('yGridSpaces');

    yCenter = screenInfoMap('yCenter');
    if strcmp(screenside, 'left')
        xCenter = screenInfoMap('leftxCenter');
        xsideoffset = 0;
    else
        xCenter = screenInfoMap('rightxCenter');
        xsideoffset = screenInfoMap('stimXpixels')*2;
    end

    %create a set of potential poitions on the grid, from 1 to however many spots there are
    gridPositions = [1:(xGridSpaces*yGridSpaces)];
    gridPositions = gridPositions(randperm(length(gridPositions)));

    twoscalegridPositions = [1:(2*4)];
    twoscalegridPositions = twoscalegridPositions(randperm(length(twoscalegridPositions)));


    gridCoordinates = generateGrid(xGridSpaces, yGridSpaces);
    twoscalegridCoordinates = generateGrid(2,4);

    [xpoints ypoints breakList] = generateNaturalCoordinateSet(calculationsMap, screenInfoMap, ...
    numberOfLoops, ellipseScale, xsideoffset, scale, gridCoordinates, gridPositions, twoscalegridCoordinates, twoscalegridPositions, calculationsMap('framesPerLoop'));

    totalpoints = numel(xpoints);

    bgRect = CenterRectOnPointd(screenInfoMap('baseRect'), xCenter, yCenter);
    Screen('FillRect', window, colorsMap('rgbgrey'), bgRect);
    savepoint = 1;

    for p = 1:totalpoints - 2
        if ~any(p == breakList) && ~any(p+1 == breakList)
            Screen('DrawLine', window, lineColor, xpoints(p), ypoints(p), ...
                xpoints(p+1), ypoints(p+1), 5);
        else
            Screen('DrawLine', window, lineColor, xpoints(p), ypoints(p), ...
                xpoints(savepoint), ypoints(savepoint), 5);
            savepoint = p+1;
        end
    end
    Screen('DrawLine', window, lineColor, xpoints(totalpoints-1), ypoints(totalpoints-1), ...
                xpoints(savepoint), ypoints(savepoint), 5);
end

% Generates the x-y coordinate pairs for the ellipses on the grid, but does not draw them.
function [xpoints ypoints breakList] = generateNaturalCoordinateSet(calculationsMap, screenInfoMap, ...
  numberOfLoops, ellipseScale, xsideoffset, scale, gridCoordinates, gridPositions, twoscalegridCoordinates, twoscalegridPositions, framesPerLoop)
    xpoints = [];
    ypoints = [];
    breakList = [];
    for i = 1:numberOfLoops
        [xEllipse, yEllipse] = plotEllipse(framesPerLoop, ellipseScale);
        [xEllipse, yEllipse] = rotateEllipse(xEllipse, yEllipse);
        xEllipse = xEllipse .* scale;
        yEllipse = yEllipse .* scale;

        % if the scale is big, use the smaller grid. otherwise, use the normal grid
        if(ellipseScale==2)
            [xEllipse, yEllipse] = transposeEllipse(screenInfoMap, xEllipse, yEllipse, xsideoffset, twoscalegridCoordinates(twoscalegridPositions(i),:));
        else
            [xEllipse, yEllipse] = transposeEllipse(screenInfoMap, xEllipse, yEllipse, xsideoffset, gridCoordinates(gridPositions(i),:));
        end

        xpoints = [xpoints xEllipse];
        ypoints = [ypoints yEllipse];
        breakList = [breakList numel(xpoints)];
    end
end

function [processedxPoints processedyPoints breakList] = generateRandomCoordinateSet(calculationsMap, screenInfoMap, ...
  numberOfLoops, ellipseScale, xsideoffset, scale, gridCoordinates, gridPositions, twoscalegridCoordinates, twoscalegridPositions, framesPerLoop)
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
            [currentxPoints, currentyPoints] = transposeEllipse(screenInfoMap, currentxPoints, currentyPoints, xsideoffset, twoscalegridCoordinates(twoscalegridPositions(gridIndex),:));
        else
            [currentxPoints, currentyPoints] = transposeEllipse(screenInfoMap, currentxPoints, currentyPoints, xsideoffset, gridCoordinates(gridPositions(gridIndex),:));
        end

        processedxPoints = [processedxPoints currentxPoints];
        processedyPoints = [processedyPoints currentyPoints];

        % reset iterators
        tempi = i;
        gridIndex = gridIndex + 1;
    end
end

% Like the above function, but for randomized breaks instead
function [] = drawRandomObjectSet(calculationsMap, screenInfoMap, colorsMap, ...
  numberOfLoops, ellipseScale, screenside, lineColor)
    scale = calculationsMap('scale');
    framesPerLoop = calculationsMap('framesPerObjectLoop');
    window = screenInfoMap('window');
    black = colorsMap('screenBlack');
    xGridSpaces = calculationsMap('xGridSpaces');
    yGridSpaces = calculationsMap('yGridSpaces');

    xpoints = [];
    ypoints = [];

    yCenter = screenInfoMap('yCenter');
    if strcmp(screenside, 'left')
        xCenter = screenInfoMap('leftxCenter');
        xsideoffset = 0;
    else
        xCenter = screenInfoMap('rightxCenter');
        xsideoffset = screenInfoMap('stimXpixels')*2;
    end

    %create a set of potential poitions on the grid, from 1 to however many spots there are
    gridPositions = [1:(xGridSpaces*yGridSpaces)];
    gridPositions = gridPositions(randperm(length(gridPositions)));

    twoscalegridPositions = [1:(2*4)];
    twoscalegridPositions = twoscalegridPositions(randperm(length(twoscalegridPositions)));


    gridCoordinates = generateGrid(xGridSpaces, yGridSpaces);
    twoscalegridCoordinates = generateGrid(2,4);

    [xpoints ypoints breakList] = generateRandomCoordinateSet(calculationsMap, screenInfoMap, ...
    numberOfLoops, ellipseScale, xsideoffset, scale, gridCoordinates, gridPositions, twoscalegridCoordinates, twoscalegridPositions, calculationsMap('framesPerLoop'));

    % drawing to screen
    totalpoints = numel(xpoints);
    bgRect = CenterRectOnPointd(screenInfoMap('baseRect'), xCenter, yCenter);
    Screen('FillRect', window, colorsMap('rgbgrey'), bgRect);
    savepoint = 1;

    for p = 1:totalpoints - 2
        if ~any(p == breakList) && ~any(p+1 == breakList)
            Screen('DrawLine', window, lineColor, xpoints(p), ypoints(p), ...
                xpoints(p+1), ypoints(p+1), 5);
        end
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
    %The orientation starts at 0, and ends at 360-360/numberOfLoops
    %This is so it doesn't make a complete circle, which would have two
    %overlapping ellipses.
    % orientation = linspace(0,360-round(360/numberOfLoops),numberOfLoops);

    %orientation calculated from above
    % loopOri=0;

    %Start with the basic, unrotated ellipse
    x = (majorAxis/2) * sin(theta) + centerX;
    y = (minorAxis/2) * cos(theta) + centerY;

    %Then rotate it
    % x = (initx-centerX)*cos(loopOri) - (inity-centerY)*sin(loopOri) + centerX;
    % y = (initx-centerX)*sin(loopOri) + (inity-centerY)*cos(loopOri) + centerY;

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

    %Move to origin
    % for m = 1:totalpoints
    %     nx(m) = xpoints(m) - xpoints(halfLoop)/2;
    %     ny(m) = ypoints(m) - ypoints(halfLoop)/2;
    % end

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

    %Move to origin
    for m = 1:totalpoints-1
        % if any(m==Breaks)
        %     petalnum = petalnum+1;
        % end
        % the point becomes the normal point MINUS the outermost part of the section
        % nx(m) = xpoints(m) - xpoints(halfLoop + (numberOfFrames * petalnum))/2;
        % ny(m) = ypoints(m) - ypoints(halfLoop + (numberOfFrames * petalnum))/2;
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


%%%%%%%%%
% FUNCTIONS FOR EVENTS

function [] = runEventsFromPoints(calculationsMap, screenInfoMap, colorsMap, alternatingSide, breakType,  alternateParametersList, alternateBreakList, constParametersList, constBreakList)
    % A bunch of info retrieval up front to avoid dealing with maps during the stimulus presentation
    parametersKeyList = calculationsMap('parametersKeyList');
    translatedParameters = calculationsMap('translatedParameters');
    leftAlternating = strcmp(alternatingSide, 'left');
    %TODO: change the time to 60s
    finalTime = datenum(clock + [0, 0, 0, 0, 0, 25]);
    blankscreenTime = calculationsMap('blankscreenTime');
    window = screenInfoMap('window');
    leftxCenter = screenInfoMap('leftxCenter');
    rightxCenter = screenInfoMap('rightxCenter');
    yCenter = screenInfoMap('yCenter');
    scale = calculationsMap('scale');
    ifi = screenInfoMap('ifi');
    minSpace = calculationsMap('minSpace');
    breakFrames = round(calculationsMap('breakTime') / ifi);
    leftAlternating = strcmp(alternatingSide, 'left');

    % set up parameter cycle
    constantParams = [constParametersList(1,:); constParametersList(2,:)];
    alternatingParams = [alternateParametersList(1,:); alternateParametersList(2,:)];
    constantBreakList = constBreakList(1,:);
    alternatingBreakList = alternateBreakList(1,:);
    alternating_xpoints = alternatingParams(1,:);
    alternating_ypoints = alternatingParams(2,:);
    constant_xpoints = constantParams(1,:);
    constant_ypoints = constantParams(2,:);
    p=3;

    [constant_xpoints, constant_ypoints] = addBreakFrames(constant_xpoints, constant_ypoints, constantBreakList, breakFrames);
    [alternating_xpoints, alternating_ypoints] = addBreakFrames(alternating_xpoints, alternating_ypoints, alternatingBreakList, breakFrames);
    % frame indexing variable f, parameter indexing variable p
    f=1;
    currentAnimationLength = max([numel(constant_xpoints), numel(alternating_ypoints)]);
    while datenum(clock) < finalTime
        % if the current f is through all the frames of the current animation
        if f >= currentAnimationLength
            % if the clock+1 second is before the final time (there's more than one second left to animate)
            if datenum(clock + [0, 0, 0, 0, 0, 1]) < finalTime
                % draw the blank screen and flip (asap flip)
                drawBlankScreen(screenInfoMap, colorsMap);
                screenInfoMap('vbl') = Screen('Flip', window, screenInfoMap('vbl') + 0.5 * ifi);
                %reset animation parameters
                constantParams = [constParametersList(p,:); constParametersList(p+1,:)];
                alternatingParams = [alternateParametersList(p,:); alternateParametersList(p+1,:)];
                constantBreakList = constBreakList(((p+1)/2),:);
                alternatingBreakList = alternateBreakList(((p+1)/2),:);
                alternating_xpoints = alternatingParams(1,:);
                alternating_ypoints = alternatingParams(2,:);
                constant_xpoints = constantParams(1,:);
                constant_ypoints = constantParams(2,:);
                p=p+2;

                [constant_xpoints, constant_ypoints] = addBreakFrames(constant_xpoints, constant_ypoints, constantBreakList, breakFrames);
                [alternating_xpoints, alternating_ypoints] = addBreakFrames(alternating_xpoints, alternating_ypoints, alternatingBreakList, breakFrames);
                % and reset f and the currentAnimationLength
                f=1;
                currentAnimationLength = max([numel(constant_xpoints), numel(alternating_ypoints)]);
                % draw the first point and flip, vbl+blankscreentime
                drawEventFrame(calculationsMap, screenInfoMap, colorsMap, constant_xpoints, constant_ypoints, alternating_xpoints, alternating_ypoints, f)
                screenInfoMap('vbl') = Screen('Flip', window, screenInfoMap('vbl') + blankscreenTime);
                % if the time is too late, do nothing;no flip or anything so the stars stay on screen
                % this way, it lasts the right amount of time without going over
            end
        % otherwisewise, draw the frame for the corresponding f
        else
            if(leftAlternating)
                drawEventFrame(calculationsMap, screenInfoMap, colorsMap, constant_xpoints, constant_ypoints, alternating_xpoints, alternating_ypoints, f)
                % then flip, min time
                screenInfoMap('vbl') = Screen('Flip', window, screenInfoMap('vbl') + 0.5 * ifi);
            else
                drawEventFrame(calculationsMap, screenInfoMap, colorsMap, alternating_xpoints, alternating_ypoints, constant_xpoints, constant_ypoints, f)
                % then flip, min time
                screenInfoMap('vbl') = Screen('Flip', window, screenInfoMap('vbl') + 0.5 * ifi);
            end
        end
        f = f+1;
    end
end

% manages running the events condition for an amount of time.
function [] = runEventsTrial(calculationsMap, screenInfoMap, colorsMap, timePerAnimation, alternatingSide, breakType)
    % A bunch of info retrieval up front to avoid dealing with maps during the stimulus presentation
    parametersKeyList = calculationsMap('parametersKeyList');
    translatedParameters = calculationsMap('translatedParameters');
    leftAlternating = strcmp(alternatingSide, 'left');
    finalTime = datenum(clock + [0, 0, 0, 0, 0, 25]);
    blankscreenTime = calculationsMap('blankscreenTime');
    window = screenInfoMap('window');
    leftxCenter = screenInfoMap('leftxCenter');
    rightxCenter = screenInfoMap('rightxCenter');
    yCenter = screenInfoMap('yCenter');
    scale = calculationsMap('scale');
    ifi = screenInfoMap('ifi');
    minSpace = calculationsMap('minSpace');
    breakFrames = round(calculationsMap('breakTime') / ifi);
    leftAlternating = strcmp(alternatingSide, 'left');

    if strcmp(breakType, 'equal')
        generationFunction = @generateNaturalCoordinateSet;
    else
        generationFunction = @generateRandomCoordinateSet;
    end

    % set up parameter cycle
    constantParams = translatedParameters(parametersKeyList(1,1),:);
    alternatingParams = translatedParameters(parametersKeyList(1,2),:);
    p=2;

    % generate the ellipse sets, and set the current animation Length (the number of frames to run to for this set) to whichever is longer)
    a_framesPerLoop = round((timePerAnimation/constantParams(1)) / ifi) + 1 - breakFrames;
    b_framesPerLoop = round((timePerAnimation/alternatingParams(1)) / ifi) + 1 - breakFrames;
    [constant_xpoints constant_ypoints constantBreakList] = generationFunction(calculationsMap, screenInfoMap, ...
    constantParams(1), constantParams(2), 0, scale, gridCoordinates, gridPositions, twoscalegridCoordinates, twoscalegridPositions, a_framesPerLoop);
    [constant_xpoints, constant_ypoints] = addBreakFrames(constant_xpoints, constant_ypoints, constantBreakList, breakFrames);
    [alternating_xpoints alternating_ypoints alternating_breakList] = generationFunction(calculationsMap, screenInfoMap, ...
    alternatingParams(1), alternatingParams(2), screenInfoMap('stimXpixels')*2, scale, gridCoordinates, gridPositions, twoscalegridCoordinates, twoscalegridPositions, b_framesPerLoop);
    [alternating_xpoints, alternating_ypoints] = addBreakFrames(alternating_xpoints, alternating_ypoints, alternating_breakList, breakFrames);
    % frame indexing variable f, parameter indexing variable p
    f=1;
    currentAnimationLength = max([numel(constant_xpoints), numel(alternating_ypoints)]);
    while datenum(clock) < finalTime
        % if the current f is through all the frames of the current animation
        if f >= currentAnimationLength
            % if the clock+1 second is before the final time
            if datenum(clock + [0, 0, 0, 0, 0, 1]) < finalTime
                % draw the blank screen and flip (asap flip)
                drawBlankScreen(screenInfoMap, colorsMap);
                screenInfoMap('vbl') = Screen('Flip', window, screenInfoMap('vbl') + 0.5 * ifi);

                % retrieve new parameters
                constantParams = translatedParameters(parametersKeyList(p,1),:);
                alternatingParams = translatedParameters(parametersKeyList(p,2),:);
                % generate a new set of animations
                a_framesPerLoop = round((timePerAnimation/constantParams(1)) / ifi) + 1 - breakFrames;
                b_framesPerLoop = round((timePerAnimation/alternatingParams(1)) / ifi) + 1 - breakFrames;
                [constant_xpoints constant_ypoints constant_breakList] = generationFunction(calculationsMap, screenInfoMap, ...
                constantParams(1), constantParams(2), 0, scale, gridCoordinates, gridPositions, twoscalegridCoordinates, twoscalegridPositions, a_framesPerLoop);
                [constant_xpoints, constant_ypoints] = addBreakFrames(constant_xpoints, constant_ypoints, constant_breakList, breakFrames);
                [alternating_xpoints alternating_ypoints alternating_breakList] = generationFunction(calculationsMap, screenInfoMap, ...
                alternatingParams(1), alternatingParams(2), screenInfoMap('stimXpixels')*2, scale, gridCoordinates, gridPositions, twoscalegridCoordinates, twoscalegridPositions, b_framesPerLoop);
                [alternating_xpoints, alternating_ypoints] = addBreakFrames(alternating_xpoints, alternating_ypoints, alternating_breakList, breakFrames);
                % and reset f and the currentAnimationLength, but increment p
                f=1;
                currentAnimationLength = max([numel(constant_xpoints), numel(alternating_ypoints)]);
                p = p+1;
                % draw the first point and flip, vbl+blankscreentime
                drawEventFrame(calculationsMap, screenInfoMap, colorsMap, constant_xpoints, constant_ypoints, alternating_xpoints, alternating_ypoints, f)
                screenInfoMap('vbl') = Screen('Flip', window, screenInfoMap('vbl') + blankscreenTime);
                % if the time is too late, do nothing;no flip or anything so the stars stay on screen
                % this way, it lasts the right amount of time without going over
            end
        % otherwisewise, draw the frame for the corresponding f
        else
            if(leftAlternating)
                drawEventFrame(calculationsMap, screenInfoMap, colorsMap, constant_xpoints, constant_ypoints, alternating_xpoints, alternating_ypoints, f)
                % then flip, min time
                screenInfoMap('vbl') = Screen('Flip', window, screenInfoMap('vbl') + 0.5 * ifi);
            else
                drawEventFrame(calculationsMap, screenInfoMap, colorsMap, alternating_xpoints, alternating_ypoints, constant_xpoints, constant_ypoints, f)
                % then flip, min time
                screenInfoMap('vbl') = Screen('Flip', window, screenInfoMap('vbl') + 0.5 * ifi);
            end
        end
        f = f+1;
    end
end

function [xpoints, ypoints] = addBreakFrames(xpoints, ypoints, breakList, breakFrames)
    breakList = sort(breakList, 'descend');
    for i = breakList
        xRepeat = repelem(xpoints(i-1), breakFrames);
        yRepeat = repelem(ypoints(i-1), breakFrames);
        % repeating a couple frames to avoid drawing any frames connecting one ellipse to another
        xpoints = [xpoints(1:i-1) xRepeat xpoints(i+1:end)];
        ypoints = [ypoints(1:i-1) yRepeat ypoints(i+1:end)];
    end
end

% Display the frameNumber-th frame of each event set. Does NOT actually flip to make it visible on the screen
function [] = drawEventFrame(calculationsMap, screenInfoMap, colorsMap, a_xpoints, a_ypoints, b_xpoints, b_ypoints, frameNumber)
    imageTexture = screenInfoMap('imageTexture');
    window = screenInfoMap('window');
    a_pointslength = numel(a_xpoints);
    b_pointslength = numel(b_xpoints);
    grey = colorsMap('rgbgrey');
    baseRect = screenInfoMap('baseRect');
    yCenter = screenInfoMap('yCenter');
    b_xpoints = b_xpoints + screenInfoMap('stimXpixels')*2;
    % draw gray background
    bgRect = CenterRectOnPointd(baseRect, screenInfoMap('leftxCenter'), yCenter);
    Screen('FillRect', window, grey, bgRect);
    bgRect = CenterRectOnPointd(baseRect, screenInfoMap('rightxCenter'), yCenter);
    Screen('FillRect', window, grey, bgRect);

    % conditional makes sure there's still frames left to display of that event
    if frameNumber <= a_pointslength
        n = frameNumber;
    else
        % if the framenumber is past, just keep drawing the star to the last point
        n = a_pointslength;
    end
    % draw the frame to the screen
    destRect = [a_xpoints(n) - 68/2, ... %left
                a_ypoints(n) - 68/2, ... %top
                a_xpoints(n) + 68/2, ... %right
                a_ypoints(n) + 68/2]; %bottom
    Screen('DrawTexture', window, imageTexture, [], destRect, 0);
    % repeat for other side
    if frameNumber <= b_pointslength
        n = frameNumber;
    else
        n = b_pointslength;
    end
    destRect = [b_xpoints(n) - 68/2, ... %left
                b_ypoints(n) - 68/2, ... %top
                b_xpoints(n) + 68/2, ... %right
                b_ypoints(n) + 68/2]; %bottom
    Screen('DrawTexture', window, imageTexture, [], destRect, 0);
end

%%%%%%%
% FUNCTIONS FOR BOTH

% generates a set of x,y points for a number (numberOfLoops) of ellipses, and arranges them into a "flower"-like shape, like in previous Telic stimuli
% very similar to the plotEllipse function, but creates more than one ellipse, and also rotates and arranges them appropriately
% This function is also used for the events and random-break objects
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

% draw a screen with blank gray background on both sides of the screen
function [] = drawBlankScreen(screenInfoMap, colorsMap)
    window = screenInfoMap('window');
    grey = colorsMap('rgbgrey');
    yCenter = screenInfoMap('yCenter');
    Screen('FillRect', window, colorsMap('screenBlack'));
    leftRect = CenterRectOnPointd(screenInfoMap('baseRect'), screenInfoMap('leftxCenter'), yCenter);
    rightRect = CenterRectOnPointd(screenInfoMap('baseRect'), screenInfoMap('rightxCenter'), yCenter);
    Screen('FillRect', window, grey, leftRect);
    Screen('FillRect', window, grey, rightRect);
end

% display an attention-grabber in the middle of the screen
% Code from this tutorial: http://peterscarfe.com/scaleddotgriddemo.html
function [] = attentionScreen(screenInfoMap, colorsMap)
    window = screenInfoMap('window');
    grey = colorsMap('rgbgrey');
    yCenter = screenInfoMap('yCenter');
    screenYpixels = screenInfoMap('screenYpixels');
    screenXpixels = screenInfoMap('stimXpixels');
    bgRect = CenterRectOnPointd(screenInfoMap('baseRect'), screenInfoMap('xCenter'), yCenter);

    dim = 5;
    [x, y] = meshgrid(-dim:1:dim, -dim:1:dim);

    % Scale this by the distance in pixels we want between each dot
    pixelScale = screenYpixels / (dim * 2 + 2);
    pixelScale = screenYpixels / (dim * 4 + 2);
    x = x .* pixelScale;
    y = y.* pixelScale;

    % Calculate the number of dots
    numDots = numel(x);

    % Make the matrix of positions for the dots into two vectors
    xPosVector = reshape(x, 1, numDots);
    yPosVector = reshape(y, 1, numDots);

    % We can define a center for the dot coordinates to be relaitive to. Here
    % we set the centre to be the centre of the screen
    dotCenter = [screenInfoMap('xCenter') screenInfoMap('yCenter')];

    % Set the color of our dot to be random
    dotColors = rand(3, numDots);

    % Set the maximum size of the dots to 25 pixels
    maxDotSize = 30;

    % Our grid will pulse by taking the absolute value of a sine wave, we
    % therefore set the amplitude of the sine wave to 1 as this will be a
    % multiplicative scale factor ranging between 0 and 1.
    % With 0 the dots will all be on top of one another. With 1 the gird will
    % have its maximum size.
    % These are the parameters for the sine wave
    % See: http://en.wikipedia.org/wiki/Sine_wave
    amplitude = .75;
    frequency = 0.05;
    angFreq = 2 * pi * frequency;
    startPhase = 0;
    time = 0;

    % Sync us and get a time stamp
    screenInfoMap('vbl') = Screen('Flip', window);
    ifi = screenInfoMap('ifi');
    waitframes = 1;

    % Loop the animation until a key is pressed
    while ~KbCheck

        % Scale the grid coordinates
        scaleFactor = abs(amplitude * sin(angFreq * time + startPhase));

        % Scale the dot size. Here we limit the minimum dot size to one pixel
        % by using the max function as PTB won't allow the dot size to scale to
        % zero (sensibly, as you'd be drawing no dots at all)
        thisDotSize = max(4, maxDotSize .* scaleFactor);


        % Screen('FillRect', window, colorsMap('rgbgrey'), bgRect);
        % Draw all of our dots to the screen in a single line of code adding
        % the sine oscilation to the X coordinates of the dots
        Screen('DrawDots', window, [xPosVector; yPosVector] .* scaleFactor,...
            thisDotSize, dotColors, dotCenter, 2);

        % Flip to the screen
        screenInfoMap('vbl')  = Screen('Flip', window, screenInfoMap('vbl') + (waitframes - 0.5) * ifi);

        % Increment the time
        time = time + ifi;

    end

end

% A blank screen at the end of the experiment; press any key to exit it
function [] = endingScreen(screenInfoMap, colorsMap)  
    window = screenInfoMap('window');
    grey = colorsMap('rgbgrey');
    yCenter = screenInfoMap('yCenter');
    screenYpixels = screenInfoMap('screenYpixels');
    screenXpixels = screenInfoMap('stimXpixels');

    screenInfoMap('vbl') = Screen('Flip', window);
    KbWait();
end


%%%%%%%%%%%%%%%%%%%
% GENERAL FUNCTIONS

function [condition] = alternatingInputcheck(condition)
    while ~strcmp(condition, 'r') && ~strcmp(condition, 'l') && ~strcmp(condition, 'right') && ~strcmp(condition, 'left')
        condition = input('Condition must be r or l. Please enter r (right alternating first) or l (left alternating first): ', 's');
    end
    if strcmp(condition, 'r') || strcmp(condition, 'right')
      condition = 'right';
    else
      condition = 'left';
    end
end

function [condition] = breakTypeInputcheck(condition)
    while ~strcmp(condition, 'nat') && ~strcmp(condition, 'unnat') && ~strcmp(condition, 'natural') && ~strcmp(condition, 'unnatural')
        condition = input('Condition must be nat or unnat. Please enter nat (natural, equally-spaced breaks) or unnat (unnatural, randomly-spaced breaks): ', 's');
    end
    if strcmp(condition, 'nat') || strcmp(condition, 'natural')
        condition = 'equal';
    else
        condition = 'random';
    end
end

function [condition] = displayTypeInputcheck(condition)
    while ~strcmp(condition, 'o') && ~strcmp(condition, 'e') && ~strcmp(condition, 'event') && ~strcmp(condition, 'object')
        condition = input('Condition must be e or o. Please enter e (events) or o (objects): ', 's');
    end
    if strcmp(condition, 'e') || strcmp(condition, 'event')
        condition = 'event';
    else
        condition = 'object';
    end
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