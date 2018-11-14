%% Description
% Read UAVSAR L-band images
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
function [Image, Header] = readUAVSAR(imagePath, imagesize, resizeImage, cropIndexes)
% readUAVSAR - Read UAVSAR L-band images
% Syntax:  [Image, Header] = readUAVSAR(imagePath, imagesize, resizeImage, cropIndexes)
%          WARNING: X = Range and Y = Azimuth
% Inputs:
%    - imagePath = path to image
%    - imagesize = size of image to read (needed absolutely) can be found
%    in annotation file (typically [66664,9426])
%    - resizeImage = If we want to crop images
%    - ceopIndexes = Indexes to crop the image if resizeImage is true :
%     [XLeftUpCorner, YleftUpCorner, XLeftBottomCorner, YRightBottomCorner]
%
% Outputs:
%   - Image = image array
%   - Header = image info

    
    if ~resizeImage
        fileID = fopen(imagePath);
        Image_real = fread(fileID, imagesize, 'float', 4);
        fseek(fileID, 4, 'bof');
        Image_imag = fread(fileID, imagesize, 'float', 4);
        Image = complex(Image_real, Image_imag);
        clear same_real
        clear same_imag
        fclose(fileID);
    else
        numberOfRows = cropIndexes(4) - cropIndexes(2);
        numberOfColumns = cropIndexes(3) - cropIndexes(1);
        Image = zeros(numberOfRows, numberOfColumns);
        parfor row=1:numberOfRows
            fileID = fopen(imagePath);
            fseek(fileID, ( (cropIndexes(2)+row-1)*imagesize(2)+cropIndexes(1))*8, 'bof');
            temp_real = fread(fileID, numberOfColumns, 'float', 4);
            fseek(fileID, ( (cropIndexes(2)+row-1)*imagesize(2)+cropIndexes(1))*8+4, 'bof');
            temp_imag = fread(fileID, numberOfColumns, 'float', 4);
            Image(row,:) = complex(temp_real, temp_imag);
            fclose(fileID);
        end
        clear temp_real
        clear temp_imag
    end
    
        
    Header = struct;
    s = size(Image);
    Header.AzCnt = s(1);
    Header.RgCnt = s(2);
    Header.RgPixelSz = 1.67;
    Header.AzPixelSz = 0.6;
    Header.NominalCenterFreq = 1.26e+9;

    
    % Spatial vectors
    x = linspace(-Header.RgCnt * Header.RgPixelSz / 2, Header.RgCnt * Header.RgPixelSz / 2 - Header.RgPixelSz, Header.RgCnt); % Range
    y = linspace(-Header.AzCnt * Header.AzPixelSz / 2, Header.AzCnt * Header.AzPixelSz / 2 - Header.AzPixelSz, Header.AzCnt); % Azimuth
    c = 3e8; % Speed of light
    kudopcentral = 0;
    kcentral = 2 * Header.NominalCenterFreq / c;
    kx = linspace(kcentral - 1 / (2 * Header.RgPixelSz), kcentral + 1 / (2 * Header.RgPixelSz) - 1 / (Header.RgPixelSz * Header.RgCnt), Header.RgCnt);
    ky = linspace(kudopcentral - 1 / (2 * Header.AzPixelSz), kudopcentral + 1 / (2 * Header.AzPixelSz) - 1 / (Header.AzPixelSz * Header.AzCnt), Header.AzCnt);
    Header.x = x;
    Header.y = y;
    Header.kx = kx;
    Header.ky = ky;
    Header.Nx = length(x);
    Header.Ny = length(y);
end


%% ------------- END CODE ----------------
