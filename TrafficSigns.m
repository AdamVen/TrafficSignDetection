% ELEX 7815 Course Project
% Author: Adam Vengroff
% Date: 11/09/2017
% Purpose: This project identifies traffic signs from images

function TrafficSigns(image, N)
% House Keeping
%tic;
%home;
%close all;
%clear;

% For debugging purposes
%image = 'dog.jpg';
%N = 10;

% Read in image
I1 = imread(image);

% Get image information
dim = size(I1);
xdim = size(I1, 2);
ydim = size(I1, 1);

% Preallocation for speed improvement
R1 = zeros(ydim, xdim, 3, 'uint8'); % Image with only red channel
B1 = zeros(ydim, xdim, 1); % Red edges extracted
T1 = zeros(ydim, xdim, 1); % Edges thickened
S1 = zeros(ydim, xdim, 3, 'uint8'); % Blank image where polygon is drawn
M1 = zeros(ydim, xdim, 3, 'uint8'); % Overlap between S1 and  T1
matchData = zeros(1e6, 4); % Store information on overlap
hsvI = rgb2hsv(I1);

% Variables
N_sides = N; % Number of sides in polygon
t = (1 /(N_sides * 2):1 / N_sides:1)' * 2 * pi; % t vector for polygon
iteration = 1; % number of times images has been checked for polygon fit

% Thresholds for extracting red pixels based on HSI
redHueLower = 10 / 255;
redHueUpper = 170 / 255;
redSatMin = 70 / 255;
redValMin =  50 / 255;

% Thresholds for classifying object as sign
highCertainty = 1.4;
mediumCertainty = 1.2;
for x = 1:1:xdim
    for y = 1:1:ydim
        % identify and extract red pixels 
        if (hsvI(y, x, 1) > redHueUpper || hsvI(y, x, 1) < redHueLower) && (hsvI(y, x, 2) > redSatMin) && (hsvI(y, x, 3) > redValMin)
            R1(y, x, 1) = I1(y, x, 1);
            R1(y, x, 2) = I1(y, x, 2);
            R1(y, x, 3) = I1(y, x, 3);
        end
    end
end

% Apply Edge Detection
E1 = edge(rgb2gray(R1), 'Prewitt', 0.07);

% Thicken Edges
for x = 2:1:xdim
    for y = 2:1:ydim
        if E1(y, x) == 1
           T1(y + 1, x - 1) = 1;
           T1(y + 0, x - 1) = 1;
           T1(y - 1, x - 1) = 1;
           T1(y + 1, x + 0) = 1;
           T1(y - 0, x + 0) = 1;
           T1(y + 1, x + 1) = 1;
           T1(y + 0, x + 1) = 1;
           T1(y - 1, x + 1) = 1;
        end
    end
end
    
% Fill in any gaps in lines and then fill in shapes
s = strel('square', 5);
C1 = imclose(T1, s);
BW2 = imfill(C1,'holes');

% Get largest objects
F1 = bwareafilt(logical(BW2), 2); % use 2 largest objects as safety margin
N1 = bwareafilt(logical(T1), 2); % Image with blank pixels removed

% Get boundaries
[B, L, N, A] = bwboundaries(E1, 'noholes'); % B = Boundaries, N = number of objects, 4 = connectivity

% Error checking - if no red pixels are found exit function
if (max(max(F1))) == 0;
    fprintf('Error with %s - no red pixels detected (threshhold too high?) \n', image)
    return;
end

% Only grab slices of images that contain pixels to improve efficiency
columns = any(N1);
myLoc = 1;

for x = 1:1:xdim
        if columns(x) == 1;
            Slice0(:, myLoc) = N1(:, x);
            Slice0Image(:, myLoc, 1) = I1(:, x, 1);
            Slice0Image(:, myLoc, 2) = I1(:, x, 2);
            Slice0Image(:, myLoc, 3) = I1(:, x, 3);
            myLoc = myLoc + 1;
        end
end

rows = any(Slice0, 2);
myLoc = 1;

for y = 1:1:ydim
        if rows(y) == 1;  
            Slice1(myLoc, :) = Slice0(y, :);
            Slice1Image(myLoc, :, 1) = Slice0Image(y, :, 1);
            Slice1Image(myLoc, :, 2) = Slice0Image(y, :, 2);
            Slice1Image(myLoc, :, 3) = Slice0Image(y, :, 3);
            myLoc = myLoc + 1;
        end
end

startPoint(:, 1) = find(rows, 1);
startPoint(:, 2) = find(columns, 1);

% Reduce resolution to further improve efficiency
multiplier = 128 / size(Slice1, 1);
if multiplier < 1
    Slice1 = imresize(Slice1, multiplier);
end

dim = size(Slice1);
xdim = size(Slice1, 2);
ydim = size(Slice1, 1);

% Cycle through polygon scale factors and position to determine if there's a good match
for scale = (xdim / 5):3:(xdim)
    for xOffset = abs(min(scale*sin(t))):3:xdim - abs(min(scale*sin(t)))
        for yOffset = abs(min(scale*cos(t)))+ 1:3:ydim - max(scale*cos(t))         
            % Get vertices
            polyCoord = [round(scale*sin(t) + xOffset) round(scale*cos(t) + yOffset)];
            polyCoord(N_sides + 1, 1) = polyCoord(1, 1);
            polyCoord(N_sides + 1, 2) = polyCoord(1, 2);
            
            % Create polygon using vertices
            j = 1;
            k = 1;
            
            for i = 1:1:(N_sides * 2) 
                if mod(i, 2) == 1 % odd
                    poly(i) = polyCoord(j, 1);
                    j = j + 1;
                else
                    poly(i) = polyCoord(k, 2);
                    k = k + 1;
                end
            end
            
            
            % Draw Polygon on fresh image 
            S1 = zeros(ydim, xdim, 3, 'uint8');
            S1 = insertShape(S1, 'Polygon', poly, 'Color', 'white');
            %imshow(S1); % uncomment to see polygons in real-time

            % Convert image to logical type
            S2 =  imbinarize(rgb2gray(S1));

            % Get overlap
            Slice1 = imcrop(Slice1,[0 0 xdim ydim]);
            M1 = Slice1 .* S2;

            % Get number of overlapped pixels in images
            commonPixels = sum(M1(:) == 1);
            
            % Save data
            matchData(iteration, :) = [commonPixels scale xOffset yOffset];
            iteration = iteration + 1;     
        end
    end
end

% Calculate best match
[M, I] = max(matchData);
scale = matchData(I(1), 2);
xOffset = matchData(I(1), 3);
yOffset = matchData(I(1), 4);

% Create polygon for Best Fit Shape overlapped on Sliced/Processed Image
% Get vertices
polyCoord = [round(scale*sin(t) + xOffset) round(scale*cos(t) + yOffset)];
polyCoord(N_sides + 1, 1) = polyCoord(1, 1);
polyCoord(N_sides + 1, 2) = polyCoord(1, 2);

% Create polygon using vertices
j = 1;
k = 1;

for i = 1:1:(N_sides * 2) 
    if mod(i, 2) == 1 % odd
        poly(i) = polyCoord(j, 1);
        j = j + 1;
    else
        poly(i) = polyCoord(k, 2);
        k = k + 1;
    end
end

% Create another polygon for Best Fit Shape overlapped on Sliced Raw Image
% Only needed if multiplier < 1
if multiplier < 1
    polyCoord2 = [round((scale*sin(t) + xOffset)/multiplier) round((scale*cos(t) + yOffset)/multiplier)];
    polyCoord2(N_sides + 1, 1) = polyCoord2(1, 1);
    polyCoord2(N_sides + 1, 2) = polyCoord2(1, 2);
else
    polyCoord2 = polyCoord;
end

% Create polygon using vertices
j = 1;
k = 1;

for i = 1:1:(N_sides * 2) 
    if mod(i, 2) == 1 % odd
        poly(i) = polyCoord2(j, 1);
        j = j + 1;
    else
        poly(i) = polyCoord2(k, 2);
        k = k + 1;
    end
end

% Determine if it's a traffic sign
certainty = matchData(I(1), 1) / scale;
if certainty > highCertainty
    str = 'YES';
elseif certainty > mediumCertainty
        str = 'Maybe';
else
        str = 'No';
end

sum(A(:) == 1);

% Plot Results
figure();

subplot(2, 4, 1), imshow(I1);
title('Original');

subplot(2, 4, 2), imshow(R1);
title('Extract Red Pixels');

subplot(2, 4, 3), imshow(E1);
title('Apply Edge Detection');

subplot(2, 4, 4), imshow(T1);
title('Fill in Edges');

subplot(2, 4, 5), imshow(F1);
title('Fill in Objects, Extract Largest Object');

subplot(2, 4, 6), imshow(S1);
hold on;
imshow(Slice1);
plot(polyCoord(:,1), polyCoord(:,2), 'g')
hold off;
title('Slice Image, Overlap Best Fit Shape');

subplot(2, 4, 7)
hold on;
imshow(Slice1Image);
plot(polyCoord2(:,1), polyCoord2(:,2), 'g')
hold off;
title('Best Fit Shape Overlapped Image');

subplot(2, 4, 8)
D1 = zeros(128, 128, 3, 'uint8');
D1 = insertText(D1, [32 32], str, 'FontSize', 33);
imshow(D1);
title('Match Detected?');

time = round(toc);
fprintf('%s - Traffic Sign Detected?: %s Elapsed Time: %d seconds \n', image, str, time)

%toc;
