%% Description
% Test of CD statistics on real data UAVSAR
%% Specifiactions
% Other m-files required: target creation functions
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
    clear variables
    clc
    close all
    set(0,'defaulttextinterpreter','latex')
    
%% Reading Image
    addpath('..', '../Detectors')
    
    tic;
    ScenePosition = [2891,28891,3491,31251]; % Scene 1
    Nx = ScenePosition(3)-ScenePosition(1);
    Ny = ScenePosition(4)-ScenePosition(2);
    T = 2;
    ImageHyperCube = zeros(Ny, Nx, 3, 2);
    imagePath = '../../../../Data/UAVSAR/SanAnd_26524_09014_007_090423_L090HH_03_BC_s4_1x1.slc';
    [Imagetemp, Header] = readUAVSAR(imagePath, [66664,9426], true, ScenePosition);
    ImageHyperCube(:,:,1,1) = Imagetemp;
    imagePath = '../../../../Data/UAVSAR/SanAnd_26524_09014_007_090423_L090HV_03_BC_s4_1x1.slc';
    [Imagetemp, Header] = readUAVSAR(imagePath, [66664,9426], true, ScenePosition);
    ImageHyperCube(:,:,2,1) = Imagetemp;
    imagePath = '../../../../Data/UAVSAR/SanAnd_26524_09014_007_090423_L090VV_03_BC_s4_1x1.slc';
    [Imagetemp, Header] = readUAVSAR(imagePath, [66664,9426], true, ScenePosition);
    ImageHyperCube(:,:,3) = Imagetemp;

    imagePath = '../../../../Data/UAVSAR/SanAnd_26524_15059_006_150511_L090HH_03_BC_s4_1x1.slc';
    [Imagetemp, Header] = readUAVSAR(imagePath, [66664,9426], true, ScenePosition);
    ImageHyperCube(:,:,1,2) = Imagetemp;
    imagePath = '../../../../Data/UAVSAR/SanAnd_26524_15059_006_150511_L090HV_03_BC_s4_1x1.slc';
    [Imagetemp, Header] = readUAVSAR(imagePath, [66664,9426], true, ScenePosition);
    ImageHyperCube(:,:,2,2) = Imagetemp;
    imagePath = '../../../../Data/UAVSAR/SanAnd_26524_15059_006_150511_L090VV_03_BC_s4_1x1.slc';
    [Imagetemp, Header] = readUAVSAR(imagePath, [66664,9426], true, ScenePosition);
    ImageHyperCube(:,:,3,2) = Imagetemp;
    
    clear Imagetemp
    
    ImageHyperCube = ImageHyperCube(end:-1:1,:,:,:);
    toc
    
%% Detectors

    PolDetectorsStruct = initPolChangeDetectors();
    markers = {'s', 'd', 'o', '*', '+', '^', 'pentagram', 'hexagram'};
    
%% Detection test

     % Non Parrallel
%     Results_Pol = computeCDOnImages(ImageHyperCube, PolDetectorsStruct);
    
    % Parrallel
    nbx = 10;
    nby = 10;
    numberOfWorkers = nbx*nby;
    [My, Mx] = size(PolDetectorsStruct.mask);
    ResParfor = cell(numberOfWorkers,1);
    SubImages = cell(numberOfWorkers,1);
    for index=1:numberOfWorkers
        %2d indexes for image
        subIntervalx = floor(linspace(0,Nx-Mx,nbx+1)) + floor(Mx/2)+1;
        subIntervaly = floor(linspace(0,Ny-My,nby+1)) + floor(My/2)+1;
        [i_y, i_x] = ind2sub([nby, nbx], index);
        indexesY = (subIntervaly(i_y) - floor(My/2):subIntervaly(i_y+1) + floor(My/2));
        indexesX = (subIntervalx(i_x) - floor(Mx/2):subIntervalx(i_x+1) + floor(Mx/2));

        % Transfering the needed data for the sub-image
        SubImages{index} = ImageHyperCube(indexesY, indexesX,:,:);
    end 


    % Compute results
    hbar = parfor_progressbar(numberOfWorkers,sprintf('Computing Pol ...'));  %create the progress bar
    parfor index=1:numberOfWorkers
        ResParfor{index} = computeCDOnImages(SubImages{index}, PolDetectorsStruct);
        hbar.iterate(1);   % update progress by one iteration
    end
    close(hbar);   %close progress bar

    Results_Pol = cell(PolDetectorsStruct.number,1);
    for d=1:PolDetectorsStruct.number
        ResTemp = zeros(Ny-My, Nx-Mx);
        for index=1:numberOfWorkers
            [i_y, i_x] = ind2sub([nby, nbx], index);
            subIntervalx = floor(linspace(0,Nx-Mx,nbx+1));
            subIntervaly = floor(linspace(0,Ny-My,nby+1));
            indexesY = (subIntervaly(i_y)+1:subIntervaly(i_y+1)+1);
            indexesX = (subIntervalx(i_x)+1:subIntervalx(i_x+1)+1);
            ResTemp(indexesY,indexesX,:) = ResParfor{index}{d};
        end
        Results_Pol{d} = ResTemp;
    end

%% ROC

    load('groundTruth_UAVSARScene1')
    groundTruth = groundTruth_UAVSARScene1(floor(My/2)+1:end-floor(My/2), floor(Mx/2)+1:end-floor(Mx/2));
    numberOfPoints = 100;
    figPol = figure;
    fontSize = 12;
    % Polarimetry
    lg = {};
    Res = Results_Pol{1}(end:-1:1,:);
    Nsimu = length(Res(:));
    Pfa=(Nsimu:-1:1)/Nsimu;
    for d=1:PolDetectorsStruct.number
        lg = cat(2, lg, PolDetectorsStruct.list(d).name);
        [Mx, My] = size(PolDetectorsStruct.mask);
        Res = real(Results_Pol{d}(end:-1:1,:));
        Nsimu = length(Res(:));
        lambda = sort(Res(:));
        index = floor(logspace(1,log10(Nsimu),numberOfPoints));
        index = unique(floor(index(:)));
        lambda = lambda(end-index+1);
        Pd = zeros(numberOfPoints,1);
        Pfa = zeros(numberOfPoints,1);
        for i_lambda=1:numberOfPoints 
            Pd(i_lambda) = sum(sum(groundTruth.*(Res > lambda(i_lambda)))) / sum(groundTruth(:));
            Pfa(i_lambda) = sum(sum(not(groundTruth).*(Res > lambda(i_lambda)))) / sum(sum(not(groundTruth)));
        end
        figure(figPol)
        pl = semilogx(Pfa, Pd, 'linestyle', '-', 'marker', markers{d});
        pl.MarkerIndices = floor(linspace(1,numberOfPoints,30));
        hold on
    end
    figure(figPol)
    title('UAVSAR Scene 1')
    legend(lg, 'interpreter', 'latex', 'fontsize', fontSize, 'location', 'northwest')
    xlabel('$\mathrm{P}_{\mathrm{FA}}$', 'interpreter', 'latex')
    ylabel('$\mathrm{P}_{\mathrm{D}}$', 'interpreter', 'latex')
    grid
    grid minor
    textString = {sprintf('$K_1=%d$, $K_2=%d$', Mx, My)};
    annotation('textbox',...
    [0.15, 0.57, 0.1, 0.1], ...
    'Units','characters',...
    'String', textString,...
    'backgroundcolor', 'w', ...
    'FitBoxToText','on',...
    'interpreter', 'latex');
   

%% Plotting
    for t=1:T
        signalHh = ImageHyperCube(:,:,1,t);
        signalHv = ImageHyperCube(:,:,2,t);
        signalVv = ImageHyperCube(:,:,3,t);
        
        %HH-VV
        Rraw = signalHh - signalVv;

        %HV+VH
        Graw = signalHv;

        %HH+VV
        Braw = signalHh + signalVv;

        %Span
        Span = abs(Rraw).^2 + abs(Braw).^2 + abs(Graw).^2;


        R = abs(Rraw);%abs(Rraw./sqrt(Span));
        G = abs(Graw);%abs(Graw./sqrt(Span));
        B = abs(Braw);%abs(Braw./sqrt(Span));

        RGB = cat(3,R,G,B);

        figure
        toShow = RGB;
        imagesc(toShow)
        axis image
        xlabel('x (m)', 'interpreter', 'latex')
        ylabel('y (m)', 'interpreter', 'latex')
%         colorbar
        title(sprintf('Image %d',t), 'interpreter', 'latex')  
        set(gca, 'ydir', 'normal')
    end

    for d=1:PolDetectorsStruct.number
        [Mx, My] = size(PolDetectorsStruct.mask);
        toShow = zeros(Ny,Nx);
        toShow(floor(My/2)+1:end-floor(My/2), floor(Mx/2)+1:end-floor(Mx/2)) = abs(Results_Pol{d});
        toShow(toShow==0) = nan;
        Amax = max(Results_Pol{d}(:));
        figure('Position', [400, 400, 600,900])
        imagesc(Header.x, Header.y, toShow)%, [DetectorsStruct.list(d).threshold, Amax])
        cbar = colorbar;
        set(cbar, 'color', 'w', 'location', 'south outside')
        title(PolDetectorsStruct.list(d).name, 'interpreter', 'latex')
        axis image
        xlabel('x (m)', 'interpreter', 'latex')
        ylabel('y (m)', 'interpreter', 'latex')
        set(gca, 'ydir', 'normal')
        axis image
        axis off
        set(gcf,'color','w');
        colormap('jet')
    end
    

 
    for d=1:PolDetectorsStruct.number
        [Mx, My] = size(PolDetectorsStruct.mask);
        toShow = real(Results_Pol{d});
  
      Affiche_map(Header.x(floor(Mx/2):end-floor(Mx/2)),Header.y(floor(My/2):end-floor(My/2)), toShow, PolDetectorsStruct.list(d).name)
        axis image
        set(gca, 'ydir', 'normal')
        axis image
        colormap('jet')
    end

%% ------------- END CODE ----------------
