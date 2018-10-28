%% Description
% A function to generate multivariate complex normal vectors with both
% covariance and pseudo-covariance parameters as described in:
% Picinbono, B. (1996). Second-order complex random vectors and normal
% distributions. IEEE Transactions on Signal Processing, 44(10), 2637–2640.
%% Specifications
% Other m-files required: None
% MAT-files required: None
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
function [ Z ] = Multivariate_ComplexNormal_Samples( p, N, mu, Gamma, C )
%Multivariate_ComplexNormal_Samples: 
% A function to generate multivariate complex normal vectors with both
% covariance and pseudo-covariance parameters as described in:
% Picinbono, B. (1996). Second-order complex random vectors and normal
% distributions. IEEE Transactions on Signal Processing, 44(10), 2637-2640.
%   Usage: [ Z ] = Multivariate_ComplexNormal_Samples( p, N, Gamma, C )
%          with: * p = size of the vectors
%                * N = number of samples
%                * mu = mean vector of size p*1
%                * Gamma = Covariance matrix of size p*p
%                * C = Pseudo-covariance matrix of size p*p
%        output: * Z = a p*N matrix where each coumn is a sample

    % Construction the 2p*2p matrix to generate real samples
    Gamma_x = 0.5 * real(Gamma + C);
    Gamma_xy = 0.5 * imag(-Gamma + C);
    Gamma_yx = 0.5 * imag(Gamma + C);
    Gamma_y = 0.5 * real(Gamma - C);
    Gamma_2r = [Gamma_x, Gamma_xy; Gamma_yx, Gamma_y];
    
    % Generate real and imaginary part
    mu_2r = [real(mu); imag(mu)].';
    v = mvnrnd(repmat(mu_2r, [N,1]), Gamma_2r).';
    X = v(1:p, :);
    Y = v(p+1:end, :);
    
    % Obtain the complex vector
    Z = X + 1i * Y;

end
%% ------------- END CODE ----------------
