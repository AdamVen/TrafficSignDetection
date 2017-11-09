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
I1 = imread('stopsign.jpg');

% Get image information
dim = size(I1);
xdim = size(I1, 1);
ydim = size(I1, 2);

% Preallocation for speed improvement
R1 = zeros(ydim, xdim, 3, 'uint8'); % Image with only red
B1 = zeros(ydim, xdim, 1); % Red image with edges extracted
T1 = zeros(ydim, xdim, 1); % Image with thickened edges
S1 = zeros(ydim, xdim, 1); % Image where shape is drawn

% Extract Red Pixels
for x = 1:1:xdim
    for y = 1:1:ydim
        if I1(y, x, 1) < (I1(y, x, 2) + I1(y, x, 3) + 25) % Offset is to get rid of yellows/oranges
            R1(y, x, 1) = 0;
            R1(y, x, 2) = 0;
            R1(y, x, 3) = 0;
        else
            R1(y, x, 1) = I1(y, x, 1);
            R1(y, x, 2) = I1(y, x, 2);
            R1(y, x, 3) = I1(y, x, 3);
        end
    end
end

% Apply Edge Detection
E1 = edge(rgb2gray(R1), 'Prewitt', 0.1);

% Thicken Edges
for x = 1:1:xdim
    for y = 1:1:ydim
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


% Get boundaries
[B,L,N,A] = bwboundaries(T1, 'noholes'); % B = Boundaries, N = number of objects

for k = 1:length(B) % loop throughs object
    for z = 1:1:length(B{k})
        B1(B{k}(z, 1), B{k}(z, 2)) = 1;
    end 
end

for scale = (xdim / 30):3:(xdim / 3)
    N_sides = 8;
    t = (1 /(N_sides * 2):1 / N_sides:1)' * 2 * pi;
    polyCoord = [round(scale*[sin(t)] + xdim / 2) round(scale*[cos(t)]) + ydim / 2]; % get vertices
    
    
end



% plot thickvertices
for x = 1:1:xdim
    for y = 1:1:ydim
        for k = 1:1:length(polyCoord)
            if and(x == polyCoord(k, 1), y == polyCoord(k, 2))
                S1(y + 1, x - 1) = 1;
                S1(y + 0, x - 1) = 1;
                S1(y - 1, x - 1) = 1;
                S1(y + 1, x + 0) = 1;
                S1(y + 0, x + 0) = 1;
                S1(y - 1, x + 0) = 1;
                S1(y + 1, x + 1) = 1;
                S1(y + 0, x + 1) = 1;
                S1(y - 1, x + 1) = 1;
            end
        end
    end
end

% % Draw vertical lines in polygon
% for k = 1:1:length(polyCoord) - 1
%     if polyCoord(k, 1) == polyCoord(k + 1, 1)
%         yCoord = polyCoord(k, 2)
%         yEndCoord = polyCoord(k, 2)
% 
%         for yCoord = polyCoord(k - 1, 2):1:yEndCoord
%             S1(yCoord, polyCoord(k - 1, 2)) = 1;
%         end
%     end
% end

    
% 
% horizontal translation
% if polyCoord(k, 2) == polyCoord(k + 1, 2)
%     // draw horizonnaly
% 
% diagonal translation
% else
%     
% // draw diagonally





figure();
imshow(S1);





% Plot Results
figure();

subplot(2, 4, 1), imshow(I1);
title('Original');

subplot(2, 4, 2), imshow(R1);
title('Red Pixels Extracted');

subplot(2, 4, 3), imshow(E1);
title('Edges');

subplot(2, 4, 4), imshow(T1);
title('Reconstructed Edges');

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

title(' Object Boundaries');

subplot(2, 4, 6), imshow(S1);
title('Shape Comparison');

toc;
