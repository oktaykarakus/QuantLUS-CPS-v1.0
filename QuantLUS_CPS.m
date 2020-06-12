%% QUANTLUS - CPS v1.0
% Line Artefact Detection Method based on the Cauchy Proximal Splitting Algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some Important Variables
%       ** params                       : Structure type parameter, which stores all the parameters
%                                         for the cauchy proximal splitting (CPS) algorithm.
%          params.lambda                : CPS algorithm step size value.
%          params.gamma                 : Cauchy penalty scale parameter
%          params.Niter                 : Maximum number of CPS iterations.
%
%       ** detParams                    : Structure type parameter, which stores all the
%               	                     parameters for the detection algorithm.
%          detParams.meanImage          : Mean value of the LUS image.
%          detParams.Fper               : F index value.
%          detParams.numHorLine         : Number of horizontal lines to be detected in addition 
%                                        to the pleural line.
%          detParams.plStart            : Pleural line start angle.
%          detParams.plLineBorder       : Pleural line search range.
%          detParams.subPlLinesBorder   : Other horizontal lines search range.
%          detParams.bLineBorder        : B-lines search range.
%
%       ** theta                        : Radon domain search angle.
%
%       ** wholeRadon                   : Radon image covering the whole theta (H and V lines)
%
%       ** lineImage1, lineImage2       : Matrices storing detected line information for horizontal  
%                                        and B-lines, respectovely. 
%
%       ** linesHor                     : Final detected plural line matrix
%
%       ** linesOtherHor                : Final detected horizontal line matrix
%
%       ** linesVer                     : Final detected plural B-line matrix:
%
%       ** finalImage                   : RBG Image that is the combination of LUS image and
%                                         detected lines, e.g. THE FINAL RESULT. 
%
%       ** simTime                      : Total simulation time.
% 
% This MATLAB Script is also calling the user-defined functions below:
%          functions:
%          1) [radonRecons, imageRecons, centredProbeImage] = CPS_Bline(IM, theta, params);
%          2) [plIdx, subPlIdx, bLineIdx, horInd_i, bInd_i, horInd_j,
%                     bInd_j, F, discardedLines] = Detect_Bline(imageRecons, theta, detParams);
%          
%          Please also see these functions help files for detailed
%          information about input and outout variables corresponding to
%          these functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LICENSE
% 
% This program is free software: you can redistribute it and/or modify it under the terms of the 
% GNU General Public License as published by the Free Software Foundation, either version 3 of 
% the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
% See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program.  If 
% not, see <https://www.gnu.org/licenses/>.
% 
% Copyright (c) Oktay Karakus,PhD 
% University of Bristol, UK
% o.karakus@bristol.ac.uk
% Version 1.0 - May 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCE
%
% [1] Oktay Karakus, Nantheera Anantrasirichai, Amazigh Aguersif, Stein Silva, Adrian Basarab, 
%     and Alin Achim. "Detection of Line Artefacts in Lung Ultrasound Imaging of COVID-19 Patients
%     via Non-Convex Regularization."
%     arXiv preprint arXiv:2005.03080, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clearvars
close all
clc
addpath('.\source functions');
addpath('.\Images');
load('testImage1.mat');              
idx = not(IM < 2);
I = IM./max(IM(:));
[r, c] = size(IM);
detParams.meanImage = mean(I(idx));
detParams.Fper = min(0.5, max(0.25, 1.5*detParams.meanImage));
detParams.numHorLine = 1; 
detParams.plStart = 60;

params.lambda = 5e-5;
params.gamma = 1*(sqrt(params.lambda)/2);
params.Niter = 250;
theta = -60:0.25:90;

tic;
[radonRecons, imageRecons, centredProbeImage] = CPS_Bline(IM, theta, params);

detParams.horLineBorder = round(length(centredProbeImage)/4);
detParams.bLineBorder = round(length(centredProbeImage)/16);

[plIdx, subPlIdx, bLineIdx, horInd_i, bInd_i, horInd_j, bInd_j,...
    validationMetrics, discardedLines] = Detect_Bline(imageRecons, theta, detParams);
simTime = toc; % 
[wholeRadon, rho] = radon(imageRecons, theta);
maxu = max(max(wholeRadon(20:end-20,20:end-20)));
wholeRadon = min(1,max(0,wholeRadon/maxu));
[h, ~] = size(wholeRadon);
iiiddd = find(theta>detParams.plStart, 1);
horInd_j = horInd_j + iiiddd;
thetaWhole = theta;
lineImageBLine = zeros(size(wholeRadon));
lineImageHorLine = zeros(size(wholeRadon));
indexes = [plIdx; subPlIdx; bLineIdx];
dummyLen = length([plIdx; subPlIdx]);
for i = 1:length(indexes)
    if i <= dummyLen
        if not(isempty(horInd_i))
            lineImageHorLine(horInd_i(indexes(i)), horInd_j(indexes(i))) = 1;
        end
    else
        if not(isempty(bInd_i))
            lineImageBLine(bInd_i(indexes(i)), bInd_j(indexes(i))) = 1;
        end
    end
end
lineImageBLine = iradon(lineImageBLine, thetaWhole);
lineImageHorLine = iradon(lineImageHorLine, thetaWhole);
if size(lineImageBLine,1)>length(centredProbeImage)
    lineImageBLine = lineImageBLine(2:end-1,2:end-1);
    lineImageHorLine = lineImageHorLine(2:end-1,2:end-1);
end
lineImageBLine = imresize(lineImageBLine, size(centredProbeImage));
lineImageBLine(1:r, :) = 0;
lineImageHorLine = imresize(lineImageHorLine, size(centredProbeImage));
lineImageHorLine(1:r, :) = 0;
lineImageBLine = max(0, lineImageBLine./max(lineImageBLine(:)));
lineImageHorLine = max(0, lineImageHorLine./max(lineImageHorLine(:)));
lineImageBLine(lineImageBLine>0) = 1;
lineImageHorLine(lineImageHorLine>0) = 1;
tempBLine = lineImageBLine.*centredProbeImage;
tempHorLine = lineImageHorLine.*centredProbeImage;
tempBLine = tempBLine > 0.000001;
tempHorLine = tempHorLine > 0.5;
linesHor =zeros(size(lineImageHorLine));
linesVer =zeros(size(lineImageBLine));
linesHor(imdilate(tempHorLine, strel('disk',4))) = 1;
linesVer(imdilate(tempBLine, strel('disk',2))) = 1;
[a, b] = find(not(linesHor == 0));
T = ones(size(linesVer));
for i = min(b):max(b)
    idx = b==i;
    temp = max(a(idx));
    if temp < round(1.5*r)
        T(1:temp, i) = 0;
    end
end
T(1:min(a), :) = 0;
clear temp a b idx
linesVer = linesVer.*T;

finalImage = centredProbeImage;
finalImage = repmat(finalImage,[1 1 3]);
finalImage(:, :, 2) = finalImage(:, :, 2) + linesVer;
finalImage(:, :, 1) = finalImage(:, :, 1) + linesHor;

figure(1);
set(gcf, 'Position', [100 100 700 500])
imagesc(thetaWhole, rho, wholeRadon); colormap gray; hold on;
h1 = plot(thetaWhole(thetaWhole<=detParams.plStart), ones(size(thetaWhole(thetaWhole<=detParams.plStart)))*detParams.bLineBorder, 'y-', 'LineWidth', 1);
h2 = plot(thetaWhole(thetaWhole<=detParams.plStart), -ones(size(thetaWhole(thetaWhole<=detParams.plStart)))*detParams.bLineBorder, 'y-', 'LineWidth', 1);
h3 = plot(thetaWhole(thetaWhole>detParams.plStart), ones(size(thetaWhole(thetaWhole>detParams.plStart)))*detParams.horLineBorder, 'y-.', 'LineWidth', 1);
h4 = plot(thetaWhole(thetaWhole>detParams.plStart), -ones(size(thetaWhole(thetaWhole>detParams.plStart)))*detParams.horLineBorder, 'y-.', 'LineWidth', 1);
for i = 1:length(indexes)
    if i == length(plIdx)
        if not(isempty(horInd_i))
            h5 = plot(thetaWhole(horInd_j(indexes(i))), rho(horInd_i(indexes(i))), 'rx', 'markersize', 9);
        end
    elseif i <= dummyLen
        if not(isempty(horInd_i))
            h6 = plot(thetaWhole(horInd_j(indexes(i))), rho(horInd_i(indexes(i))), 'r>', 'markersize', 9);
        end
    else
        if not(isempty(bInd_i))
            h7 = plot(thetaWhole(bInd_j(indexes(i))), rho(bInd_i(indexes(i))), 'go', 'markersize', 9);
        end
    end
end
title('Detected Lines in Radon Space')
xlabel('\theta (degrees)')
ylabel('r')
legend('Location','NorthWest', 'FontSize', 10, 'FontWeight', 'normal', 'Color', [0.5 0.5 0.5]);  % this sequence shows just the first 2 legends
set(gca, 'Position', [0.1, 0.1, 0.89, 0.84], 'FontSize', 12, 'FontWeight', 'bold');
legend([h1 h3 h5 h6 h7],{'B-line borders' 'Hor-Line borders' 'Pleural-Line' 'Hor-Lines' 'B-lines'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);imshowpair(IM, finalImage((r+1):end, ((2*r - c)/2 + 1):(end - (2*r - c)/2), :), 'montage');
title(['Detected Lines. (Index values of B-lines are F = [' num2str(validationMetrics) '])']);
axis tight
set(gca, 'Position', [0.0001, 0.0001, 0.99, 0.99])

% Simulation Log
clc
fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf('Line Artefacts Quantification - QuantLUS\n')
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf('Simulation time = %.3f\n', simTime);
fprintf('Total number of detected horizontal lines = %d\n', length([plIdx; subPlIdx]));
fprintf('Total number of detected B-lines = %d\n', (length(bLineIdx) + discardedLines));
fprintf('Total number of validated B-lines = %d\n', length(bLineIdx));
fprintf('Total number of discarded B-lines = %d\n', discardedLines);