% ELEX 7815 Course Project
% Author: Adam Vengroff
% Date: 11/09/2017
% Purpose: This project identifies traffic signs from images

% House Keeping
tic;
clear;
home;
close all;

% Read in image
I1 = imread('stopsign5.jpg');
% downsize image to increase efficiency
multiplier = 128 / size(I1, 1);
I1 = imresize(I1, multiplier);
I1 = imcrop(I1,[0 0 128 128]);

% Get image information
dim = size(I1);
xdim = size(I1, 1);
ydim = size(I1, 2);

% Preallocation for speed improvement
R1 = zeros(ydim, xdim, 3, 'uint8'); % Image with only red
B1 = zeros(ydim, xdim, 1); % Red image with edges extracted
T1 = zeros(ydim, xdim, 1); % Image with thickened edges
S1 = zeros(ydim, xdim, 3, 'uint8'); % Image where shape is drawn
M1 = zeros(ydim, xdim, 3, 'uint8'); % Overlap between S1 and  T1
matchData = zeros(1e6, 4); % Store information on overlap

% Variables
N_sides = 8;
t = (1 /(N_sides * 2):1 / N_sides:1)' * 2 * pi;
iteration = 1;

% Extract Red Pixels
for x = 1:1:xdim
    for y = 1:1:ydim
        if I1(x, y, 1) < (I1(x, y, 2) + I1(x, y, 3) + 25) % Offset is to get rid of yellows/oranges
            R1(x, y, 1) = 0;
            R1(x, y, 2) = 0;
            R1(x, y, 3) = 0;
        else
            R1(x, y, 1) = I1(x, y, 1);
            R1(x, y, 2) = I1(x, y, 2);
            R1(x, y, 3) = I1(x, y, 3);
        end
    end
end

% Apply Edge Detection
E1 = edge(rgb2gray(R1), 'Prewitt', 0.1);

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

% Get largest objects
T1 = bwareafilt(logical(T1), 3); % 3 largest objects

% Get boundaries
[B,L,N,A] = bwboundaries(E1, 'noholes'); % B = Boundaries, N = number of objects, 4 = connectivity

% Calculate and plot initial polygon coordinates and link to figure to
% allow for real time image updating
% polyCoord = [round((xdim / 30)*sin(t) + abs(min((xdim / 30)*sin(t))) + 1) round((xdim / 30)*cos(t) + abs(min((xdim / 30)*cos(t)))+ 1)];
% figure();
% fig = plot(polyCoord(:,1), polyCoord(:,2), 'g');
% fig.XDataSource = 'polyCoord(:,1)';
% fig.YDataSource = 'polyCoord(:,2)';
% xlim([0 xdim])
% ylim([0 ydim])

% Cycle through polygon scale factors and position to determine if there's a good match
for scale = (xdim / 30) + 10:3:(xdim / 2.9)
    for xOffset = abs(min(scale*sin(t))) + 1:3:xdim - max(scale*sin(t))
        for yOffset = abs(min(scale*cos(t)))+ 1:3:ydim - max(scale*cos(t))
            % Uncomment to see polygon updates to be seen in real time (warning: slows down process considerably)
            %refreshdata
            %drawnow
            
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
            S1 = zeros(ydim, xdim, 3, 'uint8'); % Image where shape is drawn
            S1 = insertShape(S1, 'Polygon', poly, 'Color', 'white');
            %imshow(S1);
            
            % Convert image to logical type
            S2 =  imbinarize(rgb2gray(S1));

            
            % Get overlap
            T1 = imcrop(T1,[0 0 xdim ydim]);
            M1 = T1 .* S2;
            
            % Get number of pixels in images
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
S1 = zeros(ydim, xdim, 3, 'uint8'); % Image where shape is drawn
S1 = insertShape(S1, 'Polygon', poly, 'Color', 'white');

sum(A(:) == 1);

% Plot Results
figure();

subplot(2, 4, 1), imshow(I1);
title('Original');

subplot(2, 4, 2), imshow(R1);
title('Red Pixels Extracted');

subplot(2, 4, 3), imshow(E1);
title('Edges');

subplot(2, 4, 4), imshow(T1);
title('Thickened Edges, Largest Object');

subplot(2, 4, 5)

% Plot Object Boundaries, taken from MATLAB's bwboundaries() documentation
imshow(T1); hold on;
colors=['b' 'g' 'r' 'c' 'm' 'y'];
for k=1:length(B), 
  boundary = B{k};
  cidx = mod(k,length(colors))+1;
  plot(boundary(:,2), boundary(:,1),...
       colors(cidx),'LineWidth',2);

  %randomize text position for better visibility
  rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
  col = boundary(rndRow,2); row = boundary(rndRow,1);
  h = text(col+1, row-1, num2str(L(row,col)));
  set(h,'Color',colors(cidx),'FontSize',14,'FontWeight','bold');
end

title('Object Boundaries');

subplot(2, 4, 6), imshow(S1);
title('Best Fit Shape');

subplot(2, 4, 7)
hold on;
imshow(I1);
plot(polyCoord(:,1), polyCoord(:,2), 'g')
hold off;
title('Best Fit Shape Overlapped');

toc;
