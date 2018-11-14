%% Description
% Generate K-distribution vector
%% Specifiactions
% Other m-files required: none
% MAT-files required: none
%% Authors
% Authors: Ammar Mian, Ph.D., Signal processing at Supelec SONDRA
% Email address: ammar.mian@centralesupelec.fr
%% Copyright
% Copyright 2017 CentraleSupelec
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%% ------------- BEGIN CODE --------------
function x = GenerateKdistributionData(p, N, Sigma, alpha, beta)
% Compute the Gaussian Test of detection for Covariances matrix equality
%             Inputs:
%                 * p = size of vector
%                 * N = number of Observations
%                 * Sigma = The covariance matrix of size p*p 
%                 * alpha = Gamma distribution shape parameter
%                 * beta = Gamma distribution scale parameter
%             Outputs:
%                 *  a matrix of size p*N with each observtion along column dimension


    g = sqrt(gamrnd(alpha,beta,1,N));
    g = ones(p,1)*g;
    x = g.*GenerateGaussianData(p,N,Sigma);

end

%% ------------- END CODE ----------------
