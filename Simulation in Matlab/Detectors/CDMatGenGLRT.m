%% Description
% SIRV MatGen GLRT for Change Detection on covariance matrices
%% Specifiactions
% Other m-files required: Tyler.m
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
function ratio = CDMatGenGLRT(Data, Args)
% SIRV MatGen GLRT for Change Detection on covariance matrices
% Syntax:  ratio = CDMatGenGLRT(Data, Args)
% Inputs:
%    - Data = A 3D array where the first two dimensions corresponds to an array containg
%    each observation at a given date arranged in columns. The third
%    specify the date.
%    - Args : tol and nitermax for Tyler

% Outputs:
%    ratio - the test ratio
    
    [p, N, T] = size(Data);
    xk0T = reshape(Data, p, N*T);
    Sigma_H0 = Tyler(xk0T, Args);
    iSigma_H0 = inv(Sigma_H0);
    logdenominator = 0;
    lognumerator = T*N*log(abs(det(Sigma_H0)));
    for t=1:T
        xkt = Data(:,:,t);
        Sigmat_H1 = Tyler(xkt, Args);
        lognumerator = lognumerator + sum( p*log(diag( xkt' * iSigma_H0 * xkt)) );
        logdenominator = logdenominator + N*log(abs(det(Sigmat_H1))) + ...
            sum( p*log(diag( xkt' * inv(Sigmat_H1) * xkt) ) );
    end
    
    ratio = real(lognumerator-logdenominator);

end
%% ------------- END CODE ----------------