%% Description
% Compute Pd against Bartlett distance
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
    clc;
%     clear variables;
    close all;
    rng('default');
    set(0,'defaulttextinterpreter','latex');
    addpath('..', '../Detectors')

  %% Some variables
    
    % Definitions
    markerStyle = {'s', 'd', 'o', '*', '+', '^'};
    
    % Simulation params
    rng(sum(100*clock));
    p = 3; % Size of data vector
    T = 10; % Number of Images
    N = 7; % Number of observations
    numberOfPoints = 30;%
    numberOfPointsPlot = 10;
    indexes = unique(floor(linspace(1,numberOfPoints,numberOfPointsPlot)));
    PfaVector = logspace(-3,0, numberOfPoints); % Pfa of detection test
    numberOfTrialsPfa = floor(100/min(PfaVector)); % For Monte-Carlo Trials
    Pfa_vector = (numberOfTrialsPfa:-1:1)/numberOfTrialsPfa;
    numberOfTrialsPd = 1000;
    
    
    % Change param
    t_C = 5;
    rho = 0.1; % No change
    alphaGamma = 0.3;
    betaGamma = 0.1;
    rho_C = 0.8; % Change
    beta_C = 0.3;
    gMT =  sqrt(gamrnd(alphaGamma,betaGamma,1,N));
    gMT = ones(p,1)*gMT;
    Sigma = toeplitz(rho.^(0:p-1));
    Sigma_C = toeplitz(rho_C.^(0:p-1));
    
    rhoPbbTG = [0.001,0.01, 0.1, 0.2,0.5,0.7,0.9,0.99,0.999,0.9999];%rand(T,1);
    
    % Detectors
    DetectorStruct = initTheoChangeDetectors();
    ndetectors = DetectorStruct.number;
    
    Pb = 'MT';
    
%% Pfa-lambda to select appropriate threshold

    lambda = zeros(numberOfTrialsPfa, ndetectors);
    thresholds = zeros(numberOfPoints, ndetectors);
    hbar = parfor_progressbar(numberOfTrialsPfa,'Computing Pfa-thresholds...');  %create the progress bar
    parfor trial=1:numberOfTrialsPfa

        % Generate Data Cube
        if strcmp(Pb, 'MG')
            Data = zeros(p,N,T);
            for t=1:T
                g = sqrt(gamrnd(alphaGamma,betaGamma,1,N));
                g = ones(p,1)*g;
                Data(:,:,t) = g.*GenerateGaussianData(p,N,Sigma);
            end
        elseif strcmp(Pb, 'TG')
            Data = zeros(p,N,T);
            for t=1:T
                Sigmatemp = toeplitz(rhoPbbTG(t).^[0:p-1]);
                Data(:,:,t) = gMT.*GenerateGaussianData(p,N,Sigmatemp);
            end
        else
            Data = zeros(p,N,T);
            for t=1:T
                Data(:,:,t) = gMT.*GenerateGaussianData(p,N,Sigma);
            end
        end

        for d=1:ndetectors % Loop on detectors
            detector = DetectorStruct.list(d);
            lambda(trial, d) = real(detector.detect(Data, detector.Args));
        end
        hbar.iterate(1);
        
    end
    close(hbar);


    for d=1:ndetectors
        lambdaTemp = sort(lambda(:,d));
        for i_pfa=1:numberOfPoints
            thresh = lambdaTemp(abs(PfaVector(i_pfa) - Pfa_vector) == min(abs(PfaVector(i_pfa) - Pfa_vector)));
            thresholds(i_pfa,d) = thresh(1);
        end
    end

    
%% Compute Pd
    
    Pd = zeros(numberOfPoints,ndetectors);
    h = waitbar(0,'Please wait for ROC...');      
    for i_pfa=1:numberOfPoints
        Lambdatemp = zeros(numberOfTrialsPd,ndetectors);
        parfor trial = 1:numberOfTrialsPd
            % Generation of time series
            if strcmp(Pb, 'MG')
                Data = zeros(p,N,T);
                for t=1:T
                    if t<t_C
                        R = Sigma;
                    else
                        R = Sigma_C;
                    end
                    g = sqrt(gamrnd(alphaGamma,betaGamma,1,N));
                    g = ones(p,1)*g;
                    Data(:,:,t) = g.*GenerateGaussianData(p,N,R);
                end
            elseif strcmp(Pb, 'MT')
                Data = zeros(p,N,T);
                g_C = sqrt(gamrnd(alphaGamma,beta_C,1,N));
                g_C = ones(p,1)*g_C;
                for t=1:T
                    if t<t_C
                        R = Sigma;
                        g = gMT;
                    else
                        R = Sigma_C;
                        g = g_C;
                    end
                    Data(:,:,t) = g.*GenerateGaussianData(p,N,R);
                end
            elseif strcmp(Pb, 'TG')
                Data = zeros(p,N,T);
                g_C = sqrt(gamrnd(alphaGamma,beta_C,1,N));
                g_C = ones(p,1)*g_C;
                for t=1:T
                    if t<t_C
                        g = gMT;
                    else
                        g = g_C;
                    end
                    Sigmatemp = toeplitz(rhoPbbTG(t).^[0:p-1]);
                    Data(:,:,t) = g.*GenerateGaussianData(p,N,Sigmatemp);
                end
            else
                Data = zeros(p,N,T);
                g_C = sqrt(gamrnd(alphaGamma,beta_C,1,N));
                g_C = ones(p,1)*g_C;
                for t=1:T
                    if t<t_C
                        R = Sigma;
                        g = gMT;
                    else
                        R = Sigma_C;
                        g = gMT;
                    end
                    Data(:,:,t) = g.*GenerateGaussianData(p,N,R);
                end
            end

            % Change Detection 
            for d=1:ndetectors % Loop on detectors
                detector = DetectorStruct.list(d);
                Lambdatemp(trial, d) = detector.detect(Data,detector.Args);
            end

        end
        for d=1:ndetectors
            Pd(i_pfa,d) = mean(Lambdatemp(:,d) >= thresholds(i_pfa,d));
        end
        waitbar(i_pfa/numberOfPoints, h);
    end
    close(h);
    
%% Plotting
    figure
    l = {};
    for d=1:ndetectors
        pl = semilogx(PfaVector, Pd(:,d), 'linestyle', '-', 'marker', markerStyle{d});
        pl.MarkerIndices = indexes;
        hold on
        l = cat(2,l,sprintf('%s', DetectorStruct.list(d).name));
    end
    lg = legend(l);
    set(lg, 'interpreter', 'latex', 'location', 'northeast', 'fontsize', 12)
    xlabel('$\mathrm{P}_{\mathrm{FA}}$', 'interpreter', 'latex')
    ylabel('$\mathrm{P}_{\mathrm{D}}$', 'interpreter', 'latex')
    grid
    grid minor
    textString = {sprintf('$p=%d$, $N=%d$, $T=%d$, $t_C=%d$', p,N,T,t_C)...,
        sprintf('%d Trials', numberOfTrialsPd)};
    annotation('textbox',...
    [0.54, 0.5, 0.1, 0.1], ...
    'Units','characters',...
    'String', textString,...
    'backgroundcolor', 'w', ...
    'FitBoxToText','on',...
    'interpreter', 'latex',...
    'fontsize', 12);    

