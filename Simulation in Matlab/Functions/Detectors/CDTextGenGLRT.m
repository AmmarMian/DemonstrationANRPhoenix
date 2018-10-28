%% Description
% SIRV TextGen GLRT for Change Detection on texture
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
function ratio = CDTextGenGLRT(Data, Args)
% SIRV TextGen GLRT for Change Detection on texture
% Syntax:  ratio = CDTextGenGLRT(Data, Args)
% Inputs:
%    - Data = A 3D array where the first two dimensions corresponds to an array containg
%    each observation at a given date arranged in columns. The third
%    specify the date.
%    - Args : tol and nitermax for Tyler

% Outputs:
%    ratio - the test ratio
    
    [p, N, T] = size(Data);
    Sigma_H0 = TylerMC(Data, Args);
    logdenominator = 0;
    lognumerator = 0;
    numTemp = 0;
    denomTemp = 1;
    for t=1:T
        xkt = Data(:,:,t);
        Sigmat_H0 = Sigma_H0(:,:,t);
        iSigmat_H0 = inv(Sigmat_H0);
        Sigmat_H1 = Tyler(xkt,Args);
        iSigmat_H1 = inv(Sigmat_H1);
        lognumerator = lognumerator + N*log(abs(det(Sigmat_H0)));
        logdenominator = logdenominator + N*log(abs(det(Sigmat_H1)));
        numTemp = numTemp + diag( xkt' * iSigmat_H0 * xkt);
        denomTemp = denomTemp .* diag( xkt' * iSigmat_H1 * xkt);
    end
    lognumerator = lognumerator + T*p*sum(log(numTemp));
    logdenominator = logdenominator + T*p*log(T) + p*sum(log(denomTemp));
    
    ratio = real(lognumerator-logdenominator);

end
%% ------------- END CODE ----------------