****************************************** QuantLUS - CPS v1.0 **************************************************
*****************************************************************************************************************
This source code package includes the MATLAB source codes for implementing a novel method for line artefacts quantification in lung ultrasound (LUS) images of COVID-19 patients. The method is formulated as a non-convex regularisation problem involving a sparsity-enforcing, Cauchy-based penalty function, and the inverse Radon transform. Moreover, a simple local maxima detection technique in the Radon transform domain, associated with known clinical definitions of line artefacts, is employed. The method accurately identifies both horizontal and vertical line artefacts in LUS images. In order to reduce the number of false and missed detection, the method includes a two-stage validation mechanism, which is performed in both Radon and image domains.

This package includes two folders:

	1) images		: Stores images for line artefacts quantification study. We have only shared an 
				example test image with name: testImage1.mat.
	
	2) source functions	: Stores five source functions:
		2.1) radonT.m
		2.2) BLinesStats.m
		2.3) CauchyProx.m
		2.4) CPS_Bline.m
	        2.5) Detect_Bline.m

and the main MATLAB script:

	1) QuantLUS_CPS.m

*****************************************************************************************************************
LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

Copyright (c) Oktay Karakus <o.karakus@bristol.ac.uk> 
              and 
              Alin Achim <alin.achim@bristol.ac.uk>, 
              09-06-2020, University of Bristol, UK
*****************************************************************************************************************
REFERENCE

[1] O Karakus, N Anantrasirichai, A Aguersif, S Silva, A Basarab, and A Achim. "Detection of Line Artefacts in Lung Ultrasound Imaging 		of COVID-19 Patients via Non-Convex Regularization." arXiv preprint arXiv:2005.03080, 2020.
arXiv link 	: https://arxiv.org/abs/2005.03080

[2] O Karakus, A Achim. (2020): QuantLUS - CPS v1.0 
https://doi.org/10.5523/bris.z47pfkwqivfj2d0qhyq7v3u1i.

*** For the CPS Algorithm please also refer to:

[3] O Karakus, P Mayo, and A Achim. "Convergence Guarantees for Non-Convex Optimisation with Cauchy-Based Penalties" arXiv preprint. 		arXiv preprint arXiv:2003.04798, 2020.
arXiv link 	: https://arxiv.org/abs/2003.04798

[4] O Karakus, A Achim. (2020): "Cauchy Proximal Splitting (CPS)". 	
https://doi.org/10.5523/bris.15y437loa26cr2nx8gnn3l4hzi 

*** For the radonT.m function please also refer to:

[5] O Karakus, and A Achim. (2020): "AssenSAR Wake Detector."
https://doi.org/10.5523/bris.f2q4t5pqlix62sv5ntvq51yjy
*****************************************************************************************************************

