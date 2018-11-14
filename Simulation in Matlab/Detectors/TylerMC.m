%% Description
% Modified Tyler MatCons fixed point estimator for covariance matrices
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

function [R, niter, err] = TylerMC(Data, Args)
% TylerMC: Compute Tyler fixed point estimate of covariance matrix specific
% to CCG Change detection GLRT 
% Syntax:  [R, niter] = TylerMC(Data, Args)
% Inputs:
%   - Data = A 3D array where the first two dimensions corresponds to an array containg
%    each observation at a given date arranged in columns. The third
%    specify the date.
%    Args - [Tol, Itermax] with
%       Tol = tolerance for Tyler estimation,
%       Itermax = Number of maximum iterations for Tyler estimation.
% Outputs:
%    R - The estimates: size pxpxT, first two dimensions are the matrix and
%    the third the date
%    niter - the number of iterations done
%    err - criterion (for all T matrices)

    % Variables
    tol = Args(1);
    nitermax = Args(2);
    [p, N, T] = size(Data);
    err = inf*ones(T,1);
    % Init
    R = repmat(eye(p),1,1,T);
    niter = 0;
    % Iterations
    while max(err)>tol && niter<nitermax
        a = 0;
        
        % Compute the textures estimate for each pixel (l,c) in W
        for t=1:T 
           xkt = Data(:,:,t);
           Rtemp = R(:,:,t);
           v = chol(Rtemp)' \ xkt;
           v = mean(v.*conj(v)); 
           a = a + v(ones(p,1),:); 
        end
        
        % Compute the estimate for next iteration
        for t=1:T 
           xkt = Data(:,:,t);
           Rtemp = R(:,:,t);
           y = xkt ./ sqrt(a);
           Rnew = y*y'/N; % Numerator
           Rnew = p * Rnew / sum(diag(Rnew)); % Normalize by the trace
           err(t) = norm(Rnew-Rtemp,'fro')/norm(Rtemp,'fro'); % Criterion
           % Update
           R(:,:,t) = Rnew;
        end     
        
        niter = niter+1;
    end

    if niter == nitermax 
        disp('TylerMC: did not converge')
    end

end
%% ------------- END CODE ----------------