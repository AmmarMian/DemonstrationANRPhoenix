%% Description
% wald statistic as described in :
% "On Multiple Covariance Equality Testing with Application to 
% SAR Change Detection", D.Ciuonzo, V.Carotenuto, A.De Maio
%% Specifiactions
% Other m-files required: SCM.m
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
function twald = CDGaussWald(Data, ~)
% wald statistic as described in :
% "On Multiple Covariance Equality Testing with Application to 
% SAR Change Detection", D.Ciuonzo, V.Carotenuto, A.De Maio
% Syntax:  twald = CDGaussGLRT(Data, Args)
% Inputs:
%    - Data = A 3D array where the first two dimensions corresponds to an array containg
%    each observation at a given date arranged in columns. The thirs
%    specify the date.
%    - Args : unused

% Outputs:
%    twald - the test statistic
    
    [N, K, M] = size(Data);
    vec = @(x) x(:);
    
    % Compute the estimates
    L = 0;
    O = 0;
    Q = 0;
    for m=1:M
       Rm = Data(:,:,m);
       Sigma_m1 = SCM(Rm);
       iSigma_m1 = inv(Sigma_m1);
       if m==1
           Sigma_11 = Sigma_m1;
       else
          L = L + K*trace((eye(N) - Sigma_11*iSigma_m1)^2);
          Q = Q + K*(iSigma_m1 - iSigma_m1*Sigma_11*iSigma_m1);
       end
       O = O + K*kron(iSigma_m1.', iSigma_m1);
    end
    
    twald = real(L - vec(Q)' * (O\vec(Q)));
    
end
%% ------------- END CODE ----------------