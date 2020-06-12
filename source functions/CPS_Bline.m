function [radonRecons, imageRecons, centredProbeImage] = CPS_Bline(IM, theta, params)
%
% Function performing the Cauchy proximal splitting (CPS) algorithm for 
% line artefact quantification (For details please see REFERENCE - Section
% III.)
%
% INPUTS
%       ** IM                   : LUS Image.
%
%       ** theta                : Radon domain search angle.
%
%       ** params               : Structure type parameter, which stores all the
%                                 parameters for the cauchy proximal splitting (CPS) algorithm.
%
% OUTPUTS
%       ** radonRecons          : Reconstructed Radon Image.
%       ** imageRecons          : inverse Radon transform of radonRecons.
%       ** centredProbeImage 	: Probe centred image generated from input LUS
%                                 image IM
% FUNCTIONS
%       ** radonRecons = CauchyProx(Z, gamma, lambda) is the Cauchy
%          proximal operator function.
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
%
% For the CPS Algorithm please also refer to:
% [2] O Karakus, P Mayo, and A Achim. "Convergence Guarantees for Non-Convex Optimisation with 
%     Cauchy-Based Penalties" 
%     arXiv preprint arXiv:2003.04798, 2020.
%
% [3] O Karakus, A Achim. (2020): "Cauchy Proximal Splitting (CPS)". 	
%     https://doi.org/10.5523/bris.15y437loa26cr2nx8gnn3l4hzi 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
lambda = params.lambda;
gamma = params.gamma;
Niter = params.Niter;

[r, c] = size(IM);
centredProbeImage = [zeros(r, 2*r); zeros(r, (2*r-c)/2) IM zeros(r, (2*r-c)/2)];
centredProbeImage = centredProbeImage./max(centredProbeImage(:));
C  = @(x) radonT(x, theta);
CT = @(x) radon(x, theta);
CTcentredProbeImage = CT(centredProbeImage);
radonRecons = zeros(size(radon(eye(size(centredProbeImage)), theta)));
delta_x = inf;
iter = 1;
old_X = radonRecons;
while (delta_x(iter) > 1e-3) && (iter < Niter)
    iter = iter + 1;
    Z = radonRecons - lambda*(CT(C(radonRecons)) - CTcentredProbeImage);
    radonRecons = CauchyProx(Z, gamma, lambda);
    delta_x(iter) = max(abs( radonRecons(:) - old_X(:) )) / max(abs(old_X(:)));
    old_X = radonRecons;
    if or((iter > 49 && delta_x(iter) > 1e-1), (iter > 119 && delta_x(iter) > 1e-2))
        break
    end
end
imageRecons = C(radonRecons);
imageRecons = imageRecons/max(imageRecons(:));
