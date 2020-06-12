function [plIdx, horLineIdx, bLineIdx, horInd_i, bInd_i, horInd_j, bInd_j, F, discardedLines] = Detect_Bline(imageRecons, theta, detParams)
%
% Function performing the line artefact quantification method in REFERENCE
% (See below) Section IV.
%
% INPUTS
%       ** imageRecons          : CPS Reconstructed LUS image (Probe centred).
%
%       ** theta                : Radon domain search angle.
%
%       ** detParams            : Structure type parameter, which stores all the
%          			  parameters for the detection algorithm. (Please See
%          			  the main function for details)
%
% OUTPUTS
%       ** plIdx                : Vector storing detected plural line's location in Radon space.
%
%       ** horLineIdx           : Vector storing detected horizontal lines' locations  in Radon space.
%
%       ** bLineIdx             : Vector storing detected B-lines' locations in Radon space. 
%
%       ** horInd_i, horInd_j 	: Variables stores all candidate horizontal lines' locations.
%
%       ** bInd_i, bInd_j       : Variables stores all candidate vertical lines' locations.
%
%       ** F                    : Measure index for validation of candidate B-lines.
%
%       ** discardedLines       : Variable storing the number of discarded B-lines.
%
% OTHER IMPORTANT VARIABLES
%       ** thetaAnglesHor 		: Pleural line search angles.
%
%       ** horizontalRadonCoefs	: Radon image calculated within pleural line search range.
%
%       ** thetaAnglesVer 		: B-line search angles.
%
%       ** verticalRadonCoefs 	: Radon image calculated within B-line search range.
%
%       ** lineStats            : Vector which stores some details of the detected line artefacts.
%
% FUNCTIONS
%       ** radonRecons = CauchyProx(Z, gamma, lambda) is the Cauchy proximal operator function.
%
%       ** F = BLinesStats(lineStats, thetaAnglesVer, imageRecons, verticalRadonCoefs, meanImage) 
%               is the function which calculates the measure index F.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
plStart = detParams.plStart;
numHorLine = detParams.numHorLine;
Fper = detParams.Fper;
meanImage = detParams.meanImage;
subPlLinesBorder = detParams.horLineBorder;
plLineBorder = round(subPlLinesBorder/2);
bLineBorder = detParams.bLineBorder;

thetaAnglesHor = theta(theta>plStart);
[horizontalRadonCoefs, verticalRadonAxes] = radon(imageRecons, thetaAnglesHor);
maxVal = max(max(horizontalRadonCoefs(20:end-20,20:end-20)));
horizontalRadonCoefs = min(1,max(0,horizontalRadonCoefs/maxVal));
[h, ~] = size(horizontalRadonCoefs);

Len2 = 1/mean(diff(thetaAnglesHor));
initRadonCoef = imfilter(horizontalRadonCoefs,fspecial('gaussian',7,5));%
initRadonCoef = imdilate(initRadonCoef,strel('disk', 7));

localMaxValues = imregionalmax(initRadonCoef);
localMaxStats = regionprops(localMaxValues,'Centroid');
horInd_i = zeros(length(localMaxStats),1);
horInd_j = zeros(length(localMaxStats),1);
for k = 1:length(localMaxStats)
    horInd_i(k) = round(localMaxStats(k).Centroid(2));
    horInd_j(k) = round(localMaxStats(k).Centroid(1));
end
plInd_i = horInd_i;
plInd_j = horInd_j;
idxI = or(verticalRadonAxes(plInd_i)' < -plLineBorder, verticalRadonAxes(plInd_i)' > plLineBorder)';
idxJ = or(plInd_j < Len2+1, plInd_j > size(horizontalRadonCoefs, 2)-Len2);
plInd_i(or(idxI, idxJ)) = [];
plInd_j(or(idxI, idxJ)) = [];
idxToBeRemoved = horizontalRadonCoefs((plInd_j-1)*h + plInd_i) < 0;
plInd_i(idxToBeRemoved) = [];
plInd_j(idxToBeRemoved) = [];
clear idxI idxJ th toremove

subPlInd_i = horInd_i;
subPlInd_j = horInd_j;
idxI0 = or(verticalRadonAxes(subPlInd_i)' < -subPlLinesBorder, verticalRadonAxes(subPlInd_i)' > subPlLinesBorder)';
idxI00 = and(verticalRadonAxes(subPlInd_i)' > -plLineBorder, verticalRadonAxes(subPlInd_i)' < plLineBorder)';
idxI = or(idxI0,idxI00);
idxJ = or(subPlInd_j < Len2+1, subPlInd_j > size(horizontalRadonCoefs, 2)-Len2);
subPlInd_i(or(idxI, idxJ)) = [];
subPlInd_j(or(idxI, idxJ)) = [];
idxToBeRemoved = horizontalRadonCoefs((subPlInd_j-1)*h + subPlInd_i) < 0;
subPlInd_i(idxToBeRemoved) = [];
subPlInd_j(idxToBeRemoved) = [];

for i = 1:length(plInd_i)
    Matrix(i, 1:3) = [thetaAnglesHor(plInd_j(i)), verticalRadonAxes(plInd_i(i)), horizontalRadonCoefs(plInd_i(i), plInd_j(i))];
end
for i = 1:length(subPlInd_i)
    Matrix0(i, 1:3) = [thetaAnglesHor(subPlInd_j(i)), verticalRadonAxes(subPlInd_i(i)), horizontalRadonCoefs(subPlInd_i(i), subPlInd_j(i))];
end

if not(isempty(plInd_i))
    idx = Matrix(:, 1) > plStart;
    temp = Matrix(idx, 3);
    [plValue, plIdx] = max(temp);
    clear idx temp
    if not(isempty(subPlInd_i))
        idx = Matrix0(:, 1) > plStart;
        temp = Matrix0(idx, 3);
        [horLineValue, horLineIdx] = maxk(temp, numHorLine);
        clear temp idx
        horLineIdx(not(horLineValue > 0.8*plValue)) = [];
        if not(isempty(horLineIdx))
            horLineIdx = horLineIdx + size(Matrix, 1);
        end
    else
        horLineIdx = [];
    end
else
    if not(isempty(subPlInd_i))
        idx = Matrix0(:, 1) > plStart;
        temp = Matrix0(idx, 3);
        [~, plIdx] = max(temp);
        clear idx temp
        horLineIdx = [];
    else
        horLineIdx = [];
        plIdx = [];
    end
end
horInd_i =[plInd_i; subPlInd_i];
horInd_j = [plInd_j; subPlInd_j];


clear localMaxValues initRadonCoef localMaxStats
thetaAnglesVer = theta(theta<=plStart);
verticalRadonCoefs = radon(imageRecons, thetaAnglesVer);
maxVal = max(max(verticalRadonCoefs(20:end-20,20:end-20)));
verticalRadonCoefs = min(1,max(0,verticalRadonCoefs/maxVal));
[h, ~] = size(verticalRadonCoefs);
initRadonCoef = imfilter(verticalRadonCoefs,fspecial('gaussian',7,5));%
initRadonCoef = imdilate(initRadonCoef,strel('disk',7));
localMaxValues = imregionalmax(initRadonCoef);
localMaxStats = regionprops(localMaxValues,'Centroid');
bInd_i = zeros(length(localMaxStats),1);
bInd_j = zeros(length(localMaxStats),1);
for k = 1:length(localMaxStats)
    bInd_i(k) = round(localMaxStats(k).Centroid(2));
    bInd_j(k) = round(localMaxStats(k).Centroid(1));
end
idxI2 = or(verticalRadonAxes(bInd_i)' < -bLineBorder, verticalRadonAxes(bInd_i)' > bLineBorder)';
bInd_i(idxI2) = [];
bInd_j(idxI2) = [];
th = min(max(verticalRadonCoefs((bInd_j-1)*h + bInd_i))*0.25,length(imageRecons)/4);
idxToBeRemoved = verticalRadonCoefs((bInd_j-1)*h + bInd_i) < th;
bInd_i(idxToBeRemoved) = [];
bInd_j(idxToBeRemoved) = [];

if isempty(bInd_i)
    Matrix2 = [-60 0 1];
else
    for i = 1:length(bInd_i)
        Matrix2(i, 1:3) = [thetaAnglesVer(bInd_j(i)), verticalRadonAxes(bInd_i(i)), verticalRadonCoefs(bInd_i(i), bInd_j(i))];        
    end
end

mergeAngle = 8;
mergeAngle2 = 15;
bLineIdx = (1:size(Matrix2, 1))';%find(Matrix2(:, 3) >= min(0.70, max(0.5, plValue)));%min(0.90, max(0.8, plValue)));
if not(isempty(bLineIdx))
    i = 1;
    while i<size(bLineIdx, 1)
        delta = Matrix2(bLineIdx(i + 1), 1) - Matrix2(bLineIdx(i), 1);
        delta2 = abs(Matrix2(bLineIdx(i + 1), 2) - Matrix2(bLineIdx(i), 2));
        % Merging the lines with small angle difference of mergeAngle
        if (delta <= mergeAngle) || ((delta2>20) && ((delta>mergeAngle) && delta<=mergeAngle2))
            if Matrix2(bLineIdx(i + 1), 3) >= Matrix2(bLineIdx(i), 3)
                Matrix2(bLineIdx(i), :) = [];
                bInd_i(bLineIdx(i)) = [];
                bInd_j(bLineIdx(i)) = [];
                bLineIdx = (1:size(Matrix2, 1))';
            else
                Matrix2(bLineIdx(i + 1), :) = [];
                bInd_i(bLineIdx(i+1)) = [];
                bInd_j(bLineIdx(i+1)) = [];
                bLineIdx = (1:size(Matrix2, 1))';
            end
            i = 1;
        else
            i = i + 1;
        end
    end
end

if not(isempty(bLineIdx))
    lineStats(:, 1) = bInd_i(bLineIdx);
    lineStats(:, 2) = bInd_j(bLineIdx);
    F = BLinesStats(lineStats, thetaAnglesVer, imageRecons, verticalRadonCoefs, meanImage);
    idxToBeRemoved = F<(Fper);
    discardedLines = length(find(idxToBeRemoved == 1));
    bLineIdx(idxToBeRemoved) = [];
else
    F = 0;
end