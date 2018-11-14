%% Description
% Test CFAR Behaviour of detectors
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
%     close all;
    rng(sum(100*clock));
    set(0,'defaulttextinterpreter','latex');
    addpath('..', '../ChangeDetection', '../Detectors')

%% Simulation parameters

    % Definitions
    markerStyle = {'+', 'o', '*', '.', 's', 'd', '^', 'p', 'h'};
    % M-C parameters
    numberOfTrials = 100000;
    Pfa = (numberOfTrials:-1:1)/numberOfTrials;
    % Detectors
    detectorStruct = initPolChangeDetectors();
    % Dimensions
    p = 3; % Size of data vector
    N = 81; % Number Of observations
    T = 2; % Number of dates
    
    % Plotting indexes
    numberOfPoints = 30;
    index = floor(logspace(1,log10(numberOfTrials),numberOfPoints));
    index = unique(floor(index(:)));
    
%% CFAR Matrix   
    
    Pb = 'Gaussian';  
    % Toeplitz matrix parameters
    Rau_vec = [0.001];
    Raulength = length(Rau_vec);
    % Textures gamma params
    alpha = 0.3;
    beta = 0.1;
    % Compute detection test under H0 hypothesis
    percent = 0;
    h = waitbar(percent);
    for d=1:detectorStruct.number
       waitbar(0, h, ['Computing Test for detector (', num2str(d), '/', num2str(detectorStruct.number), ') ', num2str(0),' %'])
       lambda = zeros(numberOfTrials,Raulength);
       detector = detectorStruct.list(d);
       for i_r=1:Raulength
          rau = Rau_vec(i_r);
          parfor trial=1:numberOfTrials % Trials
            % True covariance Matrix
            R = toeplitz(rau.^(0:p-1));
            % Generate observations
            if strcmp(Pb, 'MG')
                Data = zeros(p,N,T);
                for t=1:T
                    Data(:,:,t) = GenerateKdistributionData(p,N,R, alpha, beta);
                end
            elseif strcmp(Pb, 'Gaussian')
                for t=1:T
                    Data(:,:,t) = GenerateGaussianData(p,N,R);
                end
            else
                Data = zeros(p,N,T);
                g = sqrt(gamrnd(alpha,beta,1,N));
                g = ones(p,1)*g;
%                 g=1;
                for t=1:T
                    Data(:,:,t) = g.*GenerateGaussianData(p,N,R);
                end
            end
            % Compute detection test
            lambda(trial,i_r) = detector.detect(Data, detector.Args);
          end
          percent = i_r/Raulength;
          waitbar(percent, h, ['Computing Test for detector (', num2str(d), '/', num2str(detectorStruct.number), ') ', num2str(floor(100*percent)),' %']);
       end
       
       % Plotting
       lg = {};
       figure
       for i_r=1:Raulength
          lg = cat(2, lg, sprintf('$\\rho=%.2f$', Rau_vec(i_r)));
          pl = semilogy(sort(lambda(:,i_r)), Pfa, 'k', 'linewidth', 1, 'linestyle', '--', 'marker', markerStyle{i_r}, 'markersize', 10);
          pl.MarkerIndices = numberOfTrials-index+1;
          ylim([100/numberOfTrials, 1])
          hold on
       end
        hlg = legend(lg);
        set(hlg, 'interpreter', 'latex')
        if strcmp(detector.scale, 'log')
            xlabel('$\log(\lambda)$', 'interpreter', 'latex')
        else
            xlabel('$\lambda$', 'interpreter', 'latex')
        end
        ylabel('$P_{FA}$', 'interpreter', 'latex')
        title(sprintf('CFAR Matrix Behaviour for %s', detector.name), 'interpreter', 'latex')
        grid
        grid minor
    end
    close(h)

%% CFAR Matrix Pb TG
    
    Pb = 'Gaussian';  
    % Toeplitz matrix parameters
    Raulength = 4;
    Rau_vec = rand(T,Raulength);

    % Textures gamma params
    alpha = 0.3;
    beta = 0.1;
    % Compute detection test under H0 hypothesis
    percent = 0;
    h = waitbar(percent);
    for d=1:detectorStruct.number
       waitbar(0, h, ['Computing Test for detector (', num2str(d), '/', num2str(detectorStruct.number), ') ', num2str(0),' %'])
       lambda = zeros(numberOfTrials,Raulength);
       detector = detectorStruct.list(d);
       for i_r=1:Raulength
          rauvec = Rau_vec(:,i_r);
          parfor trial=1:numberOfTrials % Trials

            Data = zeros(p,N,T);
            for t=1:T
                % True covariance Matrix
                R = toeplitz(rauvec(t).^(0:p-1));
                Data(:,:,t) = GenerateGaussianData(p,N,R);
            end
                
            % Compute detection test
            lambda(trial,i_r) = detector.detect(Data, detector.Args);
          end
          percent = i_r/Raulength;
          waitbar(percent, h, ['Computing Test for detector (', num2str(d), '/', num2str(detectorStruct.number), ') ', num2str(floor(100*percent)),' %']);
       end
       
       % Plotting
       lg = {};
       figure
       for i_r=1:Raulength
          stringLg = [sprintf('$\\rho_%d=%.2f$',1,Rau_vec(1,i_r))];
          for t=2:T
              stringLg = cat(2,stringLg, sprintf(', $\\rho_%d=%.2f$',t,Rau_vec(t,i_r)));
          end
          lg = cat(2, lg, stringLg);
          pl = semilogy(sort(lambda(:,i_r)), Pfa, 'k', 'linewidth', 1, 'linestyle', '--', 'marker', markerStyle{i_r}, 'markersize', 10);
          pl.MarkerIndices = numberOfTrials-index+1;
          ylim([100/numberOfTrials, 1])
          hold on
       end
        hlg = legend(lg);
        set(hlg, 'interpreter', 'latex')
        if strcmp(detector.scale, 'log')
            xlabel('$\log(\lambda)$', 'interpreter', 'latex')
        else
            xlabel('$\lambda$', 'interpreter', 'latex')
        end
        ylabel('$P_{FA}$', 'interpreter', 'latex')
        title(sprintf('CFAR Matrix Behaviour for %s', detector.name), 'interpreter', 'latex')
        grid
        grid minor
    end
    close(h) 
    
%% CFAR texture (equal between each date)
    Pb = 'MT';
    % Toeplitz matrix parameters
    rau = 0.3;
    R = toeplitz(rau.^(0:p-1));  % True covariance Matrix
    % Textures gamma params
    alpha_vec = [0.1, 0.3,0.5,0.9];
    alphaLength = length(alpha_vec);
    beta = 0.1;
    % Compute detection test under H0 hypothesis
    percent = 0;
    h = waitbar(percent);
    for d=1:detectorStruct.number
       waitbar(0, h, ['Computing Test for detector (', num2str(d), '/', num2str(detectorStruct.number), ') ', num2str(0),' %'])
       lambda = zeros(numberOfTrials,alphaLength);
       detector = detectorStruct.list(d);
       for i_r=1:alphaLength
          alpha = alpha_vec(i_r);
          parfor trial=1:numberOfTrials % Trials
           
            % Generate observations
            if strcmp(Pb, 'MG')
                Data = zeros(p,N,T);
                for t=1:T
                    Data(:,:,t) = GenerateKdistributionData(p,N,R, alpha, beta);
                end
            else
                Data = zeros(p,N,T);
                g = sqrt(gamrnd(alpha,beta,1,N));
                g = ones(p,1)*g;
                for t=1:T
                    Data(:,:,t) = g.*GenerateGaussianData(p,N,R);
                end
            end
            % Compute detection test
            lambda(trial,i_r) = detector.detect(Data, detector.Args);
          end
          percent = i_r/alphaLength;
          waitbar(percent, h, ['Computing Test for detector (', num2str(d), '/', num2str(detectorStruct.number), ') ', num2str(floor(100*percent)),' %']);
       end
       
       % Plotting
       lg = {};
       figure
       for i_r=1:alphaLength
          lg = cat(2, lg, sprintf('$\\alpha=%.2f$', alpha_vec(i_r)));
          pl = semilogy(sort(lambda(:,i_r)), Pfa, 'k', 'linewidth', 1, 'linestyle', '--', 'marker', markerStyle{i_r}, 'markersize', 6);
          pl.MarkerIndices = numberOfTrials-index+1;
          ylim([100/numberOfTrials, 1])
          hold on
       end
        hlg = legend(lg);
        set(hlg, 'interpreter', 'latex')
        if strcmp(detector.scale, 'log')
            xlabel('$\log(\lambda)$', 'interpreter', 'latex', 'fontsize', 16)
        else
            xlabel('$\lambda$', 'interpreter', 'latex', 'fontsize', 16)
        end
        ylabel('$\mathrm{P}_{\mathrm{FA}}$', 'interpreter', 'latex', 'fontsize', 16)
        title(sprintf('CFAR texture (equal between dates) Behaviour for %s', detector.name), 'interpreter', 'latex')
        grid
        grid minor
    end
    close(h)
    
%% CFAR texture (different between each date)
    Pb = 'MG';
    % Toeplitz matrix parameters
    rau = 0.3;
    R = toeplitz(rau.^(0:p-1));  % True covariance Matrix
    % Textures gamma params
    alpha_vec = [0.1, 0.3,0.5,0.9];
    alphaLength = length(alpha_vec);
    beta = 0.1;
    % Compute detection test under H0 hypothesis
    percent = 0;
    h = waitbar(percent);
    for d=1:detectorStruct.number
       waitbar(0, h, ['Computing Test for detector (', num2str(d), '/', num2str(detectorStruct.number), ') ', num2str(0),' %'])
       lambda = zeros(numberOfTrials,alphaLength);
       detector = detectorStruct.list(d);
       for i_r=1:alphaLength
          alpha = alpha_vec(i_r);
          parfor trial=1:numberOfTrials % Trials
           
            % Generate observations
            if strcmp(Pb, 'MG')
                Data = zeros(p,N,T);
                for t=1:T
                    Data(:,:,t) = GenerateKdistributionData(p,N,R, alpha, beta);
                end
            else
                Data = zeros(p,N,T);
                g = sqrt(gamrnd(alpha,beta,1,N));
                g = ones(p,1)*g;
                for t=1:T
                    Data(:,:,t) = g.*GenerateGaussianData(p,N,R);
                end
            end
            % Compute detection test
            lambda(trial,i_r) = detector.detect(Data, detector.Args);
          end
          percent = i_r/alphaLength;
          waitbar(percent, h, ['Computing Test for detector (', num2str(d), '/', num2str(detectorStruct.number), ') ', num2str(floor(100*percent)),' %']);
       end
       
       % Plotting
       lg = {};
       figure
       for i_r=1:alphaLength
          lg = cat(2, lg, sprintf('$\\alpha=%.2f$', alpha_vec(i_r)));
          pl = semilogy(sort(lambda(:,i_r)), Pfa, 'k', 'linewidth', 1, 'linestyle', '--', 'marker', markerStyle{i_r}, 'markersize', 10);
          pl.MarkerIndices = numberOfTrials-index+1;
          ylim([100/numberOfTrials, 1])
          hold on
       end
        hlg = legend(lg);
        set(hlg, 'interpreter', 'latex')
        if strcmp(detector.scale, 'log')
            xlabel('$\log(\lambda)$', 'interpreter', 'latex')
        else
            xlabel('$\lambda$', 'interpreter', 'latex')
        end
        ylabel('$P_{FA}$', 'interpreter', 'latex')
        title(sprintf('CFAR texture (different between dates) Behaviour for %s', detector.name), 'interpreter', 'latex')
        grid
        grid minor
    end
    close(h)

%% ------------- END CODE ----------------
