%% Description
% Initialisation of detectors structures for Theoretical analysis
%% Specifiactions
% Other m-files required: Detectors function
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
function DetectorsStruct = initTheoChangeDetectors()

    TylerTol = 0.001;
    TylerIterMax = 30;
    
    GaussGLRTStruct = struct;
    GaussGLRTStruct.name = '$\hat{\Lambda}_{\mathrm{G}}$';
    GaussGLRTStruct.detect = @CDGaussGLRT;
    GaussGLRTStruct.Args = '';
    GaussGLRTStruct.scale = 'log';
    
    Gausst1Struct = struct;
    Gausst1Struct.name = '$\hat{\Lambda}_{t_1}$';
    Gausst1Struct.detect = @CDGausst1;
    Gausst1Struct.Args = '';  
    Gausst1Struct.scale = 'linear';
    
    GaussWaldStruct = struct;
    GaussWaldStruct.name = '$\hat{\Lambda}_{\mathrm{wald}}$';
    GaussWaldStruct.detect = @CDGaussWald;
    GaussWaldStruct.Args = '';  
    GaussWaldStruct.scale = 'linear';
    
    
    MatAndTextGLRTStruct = struct;
    MatAndTextGLRTStruct.name = '$\hat{\Lambda}_{\mathrm{MT}}$';
    MatAndTextGLRTStruct.detect = @CDMatAndTextGLRT;
    MatAndTextGLRTStruct.Args = [TylerTol, TylerIterMax];
    MatAndTextGLRTStruct.scale = 'log';
    
    MatGenGLRTStruct = struct;
    MatGenGLRTStruct.name = '$\hat{\Lambda}_{\mathrm{MG}}$';
    MatGenGLRTStruct.detect = @CDMatGenGLRT;
    MatGenGLRTStruct.Args = [TylerTol, TylerIterMax];
    MatGenGLRTStruct.scale = 'log';
    
    MatConstrainedGLRTStruct = struct;
    MatConstrainedGLRTStruct.name = '$\hat{\Lambda}_{\mathrm{MC}}$';
    MatConstrainedGLRTStruct.detect = @CDMatConsGLRT;
    MatConstrainedGLRTStruct.Args = [TylerTol, TylerIterMax];
    MatConstrainedGLRTStruct.scale = 'log';
    
    TextGenGLRTStruct = struct;
    TextGenGLRTStruct.name = '$\hat{\Lambda}_{\mathrm{TG}}$';
    TextGenGLRTStruct.detect = @CDTextGenGLRT;
    TextGenGLRTStruct.Args = [TylerTol, TylerIterMax];
    TextGenGLRTStruct.scale = 'log';
    
    DetectorsStruct = struct;
    DetectorsStruct.list =  [GaussGLRTStruct, Gausst1Struct, GaussWaldStruct, MatAndTextGLRTStruct, MatGenGLRTStruct, TextGenGLRTStruct];
    DetectorsStruct.number = length(DetectorsStruct.list);
    DetectorsStruct.mask = ones(5);


end
%% ------------- END CODE ----------------
