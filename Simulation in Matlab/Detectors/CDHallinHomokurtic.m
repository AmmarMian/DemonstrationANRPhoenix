%% Description
% Homokurtic statistic as described in :
% "Optimal Tests for Homogeneity of Covariance, Scale and Shape",
% Marc Hallin and Davy Paindaveine
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
function QNstar = CDHallinHomokurtic(Data, ~)
% Homokurtic statistic as described in :
% "Optimal Tests for Homogeneity of Covariance, Scale and Shape",
% Marc Hallin and Davy Paindaveine
% Syntax:  QNstar = CDHallinHomokurtic(Data, ~)
% Inputs:
%    - Data = A 3D array where the first two dimensions corresponds to an array containg
%    each observation at a given date arranged in columns. The thirs
%    specify the date.
%    - Args : unused

% Outputs:
%    QNstar - the test statistic
    
    [k,ni,m] = size(Data);
%     % From complex to real data
%     Data = cat(1, real(Data), imag(Data));
%     k = 2*k;
    n = m*ni;
    
    % SCM of all the Series
    S = 0;
    for i=1:m
       S = S + 1/m*SCM(Data(:,:,i) - mean(Data(:,:,i),2)); 
    end
    
    % Estimating \hat{k}_k^{(n)}
    temp = 0;
    iS12 = inv(sqrtm(S)); 
    for i=1:m
        mX = mean(Data(:,:,i),2);
        for j=1:ni
            temp = temp + norm(iS12*(Data(:,j,i)-mX))^4;
        end
    end
    kkn = (1/(k*(k+2)*n))*temp - 1;
    
    % Computing statistic
    QNstar = 0;
    for i1=1:m
        for i2=i1+1:m
            Si1 = SCM(Data(:,:,i1) - mean(Data(:,:,i1),2) );
            Si2 = SCM(Data(:,:,i2) - mean(Data(:,:,i2),2) );
            temp = S\(Si1-Si2);
            QNstar_i1i2 = ni/(4*(1+kkn)) * ( trace(temp^2) - (kkn/((k+2)*kkn+2)) * trace(temp)^2 );
            
            QNstar = QNstar +  (2*ni) * QNstar_i1i2;
        end
    end
    QNstar = log(abs((1/n) * QNstar));
end
%% ------------- END CODE ----------------