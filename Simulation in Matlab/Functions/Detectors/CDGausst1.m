%% Description
% t1 statistic as described in :
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
function t1 = CDGausst1(Data, ~)
% t1 statistic as described in :
% "On Multiple Covariance Equality Testing with Application to 
% SAR Change Detection", D.Ciuonzo, V.Carotenuto, A.De Maio
% Syntax:  t1 = CDGaussGLRT(Data, Args)
% Inputs:
%    - Data = A 3D array where the first two dimensions corresponds to an array containg
%    each observation at a given date arranged in columns. The thirs
%    specify the date.
%    - Args : unused

% Outputs:
%    t1 - the test statistic
    
    [N, K, M] = size(Data);
    Sigma_10 = SCM(reshape(Data, N,M*K));
    iSigma_10 = inv(Sigma_10);
    t1 = 0;
    for m=1:M
        Rm = Data(:,:,m);
        Sigma_m1 = SCM(Rm);
        t1 = t1 + trace( (iSigma_10 * Sigma_m1)^2 )/M;
    end

end
%% ------------- END CODE ----------------