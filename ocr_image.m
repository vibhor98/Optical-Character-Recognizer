
colorImage = imread('airport3.jpg');
I = rgb2gray(colorImage);
[mserRegions, mserConnComp] = detectMSERFeatures(I, ...
    'RegionAreaRange', [200 8000], 'ThresholdDelta', 0.8);

figure
imshow(I)
hold on
plot(mserRegions, 'showPixelList', true, 'showEllipses', false)
title('MSER Regions')
hold off

% Removing non-text region based on Geometric properties.
% Using regionProps to measure MSER Properties
mserStats = regionprops(mserConnComp, 'BoundingBox', 'Eccentricity',...
    'Solidity', 'Extent', 'Euler', 'Image');

% Compute the aspect ratio based on BoundingBox
bbox = vertcat(mserStats.BoundingBox);
w = bbox(:,3);
h = bbox(:,4);
aspectRatio = w./h;

% Threshold the data to determine which regions to remove
filter = aspectRatio' > 3;
filter = filter | [mserStats.Eccentricity] > .995;
filter = filter | [mserStats.Solidity] < .3;
filter = filter | [mserStats.Extent] < 0.2 | [mserStats.Extent] > 0.9;
filter = filter | [mserStats.EulerNumber] < -4;

% Remove regions
mserStats(filter) = [];
mserRegions(filter) = [];

% Show remaining regions
figure
imshow(I)
hold on
plot(mserRegions, 'showPixelList', true, 'showEllipses', false)
title('After removing non-text regions based on Geometric Properties')
hold off

% Removing the non-text region using stroke width variation
% It is calculated using distance transform and binary thinning operation.

% First find binary image of regions and then pad it to avoid boundary
% effects.
regionImage = mserStats(6).Image;
regionImage = padarray(regionImage, [1 1]);

% distance transform
distanceImage = bwdist(~regionImage);
% morphological - thinning opertion
skeletonImage = bwmorph(regionImage, 'thin', inf);

strokeWidthImage = regionImage;
strokeWidthImage(~skeletonImage) = 0;

% Plot regionImage with strokeWidthImage
figure
subplot(1, 2, 1)
imagesc(regionImage)
title('Region Image')

subplot(1, 2, 2)
imagesc(strokeWidthImage)
title('Stroke Width Image')

% Computing the Stroke Width Metric
strokeWidthValues = distanceImage(skeletonImage);
strokeWidthMetric = std(strokeWidthValues) / mean(strokeWidthValues);

% Threshold the Stroke Width Metric to remove non-text regions
strokeWidthThreshold = 0.4;
strokeWidthFilterIdx = strokeWidthMetric > strokeWidthThreshold;

% Loop over every region to detect the non-textual region
for j = 1:numel(mserStats)
    regionImage = mserStats(j).Image;
    regionImage = padarray(regionImage, [1 1], 0);

    distanceImage = bwdist(~regionImage);
    skeletonImage = bwmorph(regionImage, 'thin', inf);

    strokeWidthValues = distanceImage(skeletonImage);
    strokeWidthMetric = std(strokeWidthValues) / mean(strokeWidthValues);
    strokeWidthFilterIdx(j) = strokeWidthMetric > strokeWidthThreshold;
end

% Remove regions based on the Stroke Width Variation
mserRegions(strokeWidthFilterIdx) = [];
mserStats(strokeWidthFilterIdx) = [];

% Show remaining regions
figure
imshow(I)
hold on
plot(mserRegions, 'showPixelList', true, 'showEllipses', false)
title('After removing non-text regions based on Stroke Width Variation')
hold off

% Merging separate letters into single word or line.

% Get bounded boxes for all regions
bbox = vertcat(mserStats.BoundingBox);
% Convert (x, y, width, height) notation to (xmin, ymin, xmax, ymax)
xmin = bbox(:,1);
ymin = bbox(:,2);
xmax = xmin + bbox(:,3) - 1;
ymax = ymin + bbox(:,4) - 1;

% Expand the bounded boxes by small amount
expansionAmount = 0.05;
xmin = (1-expansionAmount) * xmin;
xmax = (1+expansionAmount) * xmax;
ymin = (1-expansionAmount) * ymin;
ymax = (1+expansionAmount) * ymax;

% Clip the boxes to be within the image boundaries
xmin = max(xmin, 1);
ymin = max(ymin, 1);
xmax = min(xmax, size(I, 2));
ymax = min(ymax, size(I, 1));

% Show the expanded bounding boxes
ExpandedBBoxes = [xmin ymin xmax-xmin+1 ymax-ymin+1];
IExpandedBBoxes = insertShape(colorImage, 'Rectangle', ExpandedBBoxes, 'LineWidth', 3);

% Rectangles are expanded and then drawn so that they overlap. The
% overlapping rectangles denote the chars. in one line.
figure
imshow(IExpandedBBoxes)
title('Expanded Bounding Boxes text')

% Compute the overlap ratio to check for inline text to combine them
% forming a word
overlapRatio = bboxOverlapRatio(ExpandedBBoxes, ExpandedBBoxes);

% Set the overlap ratio b/w Bounded box and itself to 0 in order to
% simplify the graph
n = size(overlapRatio, 1);
overlapRatio(1:n+1:n^2) = 0;

% Create a graph
g = graph(overlapRatio);
% Find the connected text regions within the graph
componentIndices = conncomp(g);

% Output of conncomp is all the indices of the connected boxes that make
% up a single word or inline text. Merge overlapping boxes into a single
% box by computing the min and max dimensions of a box.
xmin = accumarray(componentIndices', xmin, [], @min);
ymin = accumarray(componentIndices', ymin, [], @min);
xmax = accumarray(componentIndices', xmax, [], @max);
ymax = accumarray(componentIndices', ymax, [], @max);

% Compose the merged bounded boxes using (x, y, width, height) format.
textBBoxes = [xmin ymin xmax-xmin+1 ymax-ymin+1];

% Remove overlapping bins containing only 1 text
numRegionsInGroup = histcounts(textBBoxes);

% histcounts divides the given array into bins of overlapping group
%textBBoxes(numRegionsInGroup==1, :) = [];
% Remove all boxes with no overlapping

% Show final text detection result
ITextRegion = insertShape(colorImage, 'Rectangle', textBBoxes, 'LineWidth', 3);

figure
imshow(ITextRegion)
title('Detected text')

% Recognize characters within the boxes
ocrText = ocr(I, textBBoxes);
[ocrText.Text]
%class([ocrText.Text])

userPrompt = 'What do you want the computer to say?';
titleBar = 'Text to Speech';
%defaultString = 'Hello World!  MATLAB is an awesome program!';
defaultString = regexprep(char(strcat([ocrText.Text])), '\n+', '');
class(defaultString)
caUserInput = inputdlg(userPrompt, titleBar, 1, {defaultString});
if isempty(caUserInput)
	return;
end % Bail out if they clicked Cancel.
caUserInput = char(caUserInput); % Convert from cell to string.
system( sprintf('say "%s"', caUserInput) )
