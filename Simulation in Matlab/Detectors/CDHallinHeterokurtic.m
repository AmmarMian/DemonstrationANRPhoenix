%% Description
% Heterokurtic statistic as described in :
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
function QNcross = CDHallinHeterokurtic(Data, ~)
% Heterokurtic statistic as described in :
% "Optimal Tests for Homogeneity of Covariance, Scale and Shape",
% Marc Hallin and Davy Paindaveine
% Syntax:  QNcross = CDHallinHeterokurtic(Data, ~)
% Inputs:
%    - Data = A 3D array where the first two dimensions corresponds to an array containg
%    each observation at a given date arranged in columns. The thirs
%    specify the date.
%    - Args : unused

% Outputs:
%    QNcross - the test statistic
    

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
    
    % Estimating \hat{E}_{k,i}^{(n)} and \hat{C}_{k,i}^{(n)} 
    Ekin = zeros(m,1);
    Ckin = zeros(m,1);
    Ekn = 0;
    Ckn = 0;
    iS12 = inv(sqrtm(S)); 
    for i=1:m
        mX = mean(Data(:,:,i),2);
        for j=1:ni
            Ekin(i) = Ekin(i) + norm(iS12*(Data(:,j,i)-mX))^4 / ni;
        end
        Ckin(i) = Ekin(i) - k^2;
        Ekn = Ekn + 1/(m*Ekin(i));
        Ckn = Ckn + 1/(m*Ckin(i));
    end
    Ekn = 1/Ekn;
    Ckn = 1/Ckn;

    
    % Computing statistic
    QNcross = 0;
    for i1=1:m
        for i2=i1+1:m
            Si1 = SCM(Data(:,:,i1) - mean(Data(:,:,i1),2) );
            Si2 = SCM(Data(:,:,i2) - mean(Data(:,:,i2),2) );
            temp = S\(Si1-Si2);
            QNcross_i1i2 = (Ckn/(Ckin(i1)*Ckin(i2)))*trace(temp)^2 + ...
                          ((k*(k+2)*Ekn)/(2*Ekin(i1)*Ekin(i2)))*(trace(temp^2) - (trace(temp)^2)/k);
            
            QNcross = QNcross +  (ni^2) * QNcross_i1i2;
        end
    end
    QNcross = log(abs((1/n) * QNcross));
end
%% ------------- END CODE ----------------