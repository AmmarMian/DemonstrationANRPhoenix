%% Description
% Generate gaussian data
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
function x = GenerateGaussianData(p,N,Sigma)
% A function that generate Gaussian multivaraite data according to a given covaraince matrix.
%             Inputs:
%                 * p = size of vector
%                 * N = number of Observations
%                 * Sigma = The covariance matrix of size p*p 
%             Outputs:
%                 * a matrix of size p*N with each observtion along column dimension

    x = 1/sqrt(2)*(randn(p,N) + 1i*randn(p,N));
    C = chol(Sigma)';
    x = C*x;
end

%% ------------- END CODE ----------------
