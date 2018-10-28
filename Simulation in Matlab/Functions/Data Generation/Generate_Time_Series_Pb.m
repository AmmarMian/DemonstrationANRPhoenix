%% Description
% A function to generate multivariate time series under conditions
% presented in the paper:
% "New Robust Statistics for Change Detection in Time Series of Multivariate 
%  SAR Images", Transactions on Signal Processing
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
function [ Z ] = Generate_Time_Series_Pb( p, N, T, Sigma, Pb, alpha, beta )
% Generate_Time_Series_Pb:
% A function to generate multivariate time series under conditions
% presented in the paper:
% "New Robust Statistics for Change Detection in Time Series of Multivariate 
%  SAR Images", Transactions on Signal Processing
%   Usage: [ Z ] = Generate_Time_Series_Pb( p, N, Sigma, Pb, alpha, beta )
%          with: * p = size of the vectors
%                * N = number of samples
%                * T = number of dates
%                * Sigma = Covariance matrix of size p*p*T (each for a
%                date)
%                * Pb = either 'MatAndTex', 'Mat', 'Tex' or 'Gaussian'
%                * alpha = scale parameter for Gamma distribution. 1D Array of
%                size T
%                * beta = shape parameter for Gamma distribution. 1D Array of
%                size T
%        output: * Z = a p*N*T matrix where the first dimension is the
%        multivariate vector, the second are the samples for each date and
%        the third corresponds to the different dates


    Z = zeros(p,N,T);
    
    
    switch Pb
        case 'Gaussian' 
            for t=1:T
                Z(:,:,t) = Multivariate_ComplexNormal_Samples(p, N, zeros(p,1), Sigma(:,:,t), zeros(p)); 
            end
        case {'MatAndTex', 'Tex'}
            alpha_t = nan;
            beta_t = nan;
            for t=1:T % We update the texture only if the parameters have changed
                if (alpha(t)~=alpha_t) || (beta(t)~=beta_t)
                    alpha_t = alpha(t);
                    beta_t = beta(t);
                    g = sqrt(gamrnd(alpha_t,beta_t,1,N));
                    g = ones(p,1)*g;
                end
                Z(:,:,t) = g.*Multivariate_ComplexNormal_Samples(p, N, zeros(p,1), Sigma(:,:,t), zeros(p)); 
            end
        case 'Mat'
            for t=1:T % We update the texture every time
                alpha_t = alpha(t);
                beta_t = beta(t);
                g = sqrt(gamrnd(alpha_t,beta_t,1,N));
                g = ones(p,1)*g;
                Z(:,:,t) = g.*Multivariate_ComplexNormal_Samples(p, N, zeros(p,1), Sigma(:,:,t), zeros(p)); 
            end
        otherwise
            warning('Problem of detection not understood. Zeros by default.')
    end

end
%% ------------- END CODE ----------------
