%% Description
% Gaussian GLRT for Change Detection on covariance matrices as described in
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
function ratio = CDGaussGLRT(Data, ~)
% CDGaussGLRT: Compute Gaussian GLRT for CD detection Test
% Syntax:  ratio = CDGaussGLRT(Data, Args)
% Inputs:
%    - Data = A 3D array where the first two dimensions corresponds to an array containg
%    each observation at a given date arranged in columns. The thirs
%    specify the date.
%    - Args : unused

% Outputs:
%    ratio - the test ratio
    
    [p, N, T] = size(Data);
    Sigma_H0 = 0;
    logdenominator = 0;
    for t=1:T
        xkt = Data(:,:,t);
        Sigmat_H1 = SCM(xkt);
        Sigma_H0 = Sigma_H0 + Sigmat_H1/T;
        logdenominator = logdenominator + N*log(abs(det(Sigmat_H1)));
    end

    ratio = real(T*N*log(abs(det(Sigma_H0))) - logdenominator);

end
%% ------------- END CODE ----------------