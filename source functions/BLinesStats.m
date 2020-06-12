function F = BLinesStats(lineStats, thetaAnglesVer, imageRecons, verticalRadonCoefs, imageMean)
% Function for calcuating the validation measure index. (For details please
% see REFERENCE - Section IV-C.)
%
% INPUTS
%       ** lineStats            : Vector which stores some details of the detected
%                                 line artefacts.
%
%       ** thetaAnglesVer       : Radon image calculated within B-line search range.
%
%       ** imageRecons          : CPS Reconstructed LUS image (Probe centred).
%
%       ** verticalRadonCoefs 	: Radon image calculated within B-line search range.
%
%       ** imageMean            : Mean value of the LUS.
%                                                 
% OUTPUTS
%       ** F                    : Measure index for validation of candidate B-lines.
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
[h, w] = size(verticalRadonCoefs);
L = size(lineStats, 1);
clear idx
F = zeros(1, L);
for i = 1:L
    indI = lineStats(i, 1);
    indJ = lineStats(i, 2);
    linesmatV = zeros(h, w);
    linesmatV((indJ-1)*h + indI) = 1;
    linesmatV = iradon(linesmatV, thetaAnglesVer);
    if size(linesmatV,1)>length(imageRecons)
        linesmatV = linesmatV(2:end-1,2:end-1);
    end
    lineSize = 2;
    linesmatV = imdilate(linesmatV/max(linesmatV(:)),strel('disk',lineSize));
    linesmatV = imresize(linesmatV, size(imageRecons));
    linesmatV = max(0, linesmatV./(max(linesmatV(:))));
    idx = not(linesmatV == 0);
    half = length(imageRecons)/2;
    idxDown = [zeros(half, 2*half); idx(half+1:end, :)];
    BLine = idxDown;
    BLineIntensity = imageRecons(not(BLine == 0));
    F(i) = mean(BLineIntensity)/imageMean - 1;
end