%% Description
% Test the convergence properties of new fixed point estimators
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
    clear variables;
    close all;
    rng(sum(100*clock));
    set(0,'defaulttextinterpreter','latex');
    addpath('..', '../Hyperimages', '../Target Dictionary', '../Detectors')

%% Some variables
  
    % Definitions
    markerStyle = {'+', 'o', '*', '.', 's', 'd', '^', 'p', 'h'};
    
    % Simulation params
    rng(sum(100*clock));
    p = 3; % Size of data vector
    N = 25;
    T = 3; % Number of dates
    numberOfTrials = 1000;
     
    % Parameters
    rhoH0 = 0.7;
    SigmaH0 = toeplitz(rhoH0.^(0:p-1));
    rhoH1 = rand(T,1);
    alpha = 0.3;
    beta = 0.1;
    
    % Fixed-point params
    tol = 0;
    nitermax_vector = floor(logspace(0,3,30));

%% Compute Criterion

    errMTH0 = zeros(length(nitermax_vector),1);
    errMGH0 = zeros(length(nitermax_vector),1);
    errMTH1 = zeros(length(nitermax_vector),1);
    errMGH1 = zeros(length(nitermax_vector),1);
    errMCH0 = zeros(length(nitermax_vector),T);
    errMCH1 = zeros(length(nitermax_vector),T);
    h = waitbar(0,'Computing Error...');
    for iter=1:length(nitermax_vector) % Loopon the number of observations
        nitermax = nitermax_vector(iter);
        
        errMGH0temp = zeros(numberOfTrials,1);
        errMTH0temp = zeros(numberOfTrials,1);
        errMGH1temp = zeros(numberOfTrials,1);
        errMTH1temp = zeros(numberOfTrials,1);
        errMCH1temp = cell(numberOfTrials,1);
        errMCH0temp = cell(numberOfTrials,1);
        parfor trial=1:numberOfTrials % M-C Trial
                       
            % Case Same matrices
            Data = zeros(p,N,T);
            for t=1:T 
                tau = sqrt(gamrnd(alpha,beta,1,N));
                tau = ones(p,1)*tau;
                SigmaH1 = toeplitz(rhoH0.^(0:p-1));
                Data(:,:,t) = tau.*GenerateGaussianData(p,N,SigmaH1);
            end

            [R, niter, err] = TylerMT(Data, [tol, nitermax]);
            errMTH0temp(trial) = err; 
            [R, niter, err] = Tyler(reshape(Data,p,N*T), [tol, nitermax]);
            errMGH0temp(trial) = err; 
            [Rtemp, niter, err] = TylerMC(Data, [tol, nitermax]);
            errMCH0temp{trial} = err.';
            
            % Case Different matrices
            Data = zeros(p,N,T);
            for t=1:T
                rho = rhoH1(t);
                SigmaH1 = toeplitz(rhoH0.^(0:p-1));
                tau = sqrt(gamrnd(alpha,beta,1,N));
                tau = ones(p,1)*tau;
                Data(:,:,t) = tau.*GenerateGaussianData(p,N,SigmaH1);
            end
            [Rtemp, niter, err] = TylerMC(Data, [tol, nitermax]);
            errMCH1temp{trial} = err.';
            [R, niter, err] = TylerMT(Data, [tol, nitermax]);
            errMTH1temp(trial) = err;
            [R, niter, err] = Tyler(reshape(Data,p,N*T), [tol, nitermax]);
            errMGH1temp(trial) = err; 
           
        end
        errMTH0(iter) = mean(errMTH0temp);
        errMGH0(iter) = mean(errMGH0temp);
        errMCH0(iter,:) = mean(cell2mat(errMCH0temp),1);
        errMTH1(iter) = mean(errMTH1temp);
        errMGH1(iter) = mean(errMGH1temp);
        errMCH1(iter,:) = mean(cell2mat(errMCH1temp),1);
        waitbar(iter/length(nitermax_vector),h);
    end
    close(h);
    
%% Plotting

    % H0
    figure
    loglog(nitermax_vector, errMTH0, 'linestyle', '--', 'marker', 's')
    hold on
    loglog(nitermax_vector, errMGH0, 'linestyle', ':', 'marker', 'd')
    lg = {'$\hat{\Sigma}_{0}^{\mathrm{MT}}$', '$\hat{\Sigma}_{0}^{\mathrm{MG}}$'};
    for t=1:T
        loglog(nitermax_vector, errMCH0(:,t), 'linestyle', '-.', 'marker', 'o')
        lg = cat(2,lg, sprintf('$\\hat{\\Sigma}_{%d}^{\\mathrm{Tex}}$', t));
    end
    h = legend(lg);
    set(h, 'interpreter', 'latex', 'fontsize', 15)
    xlabel('$n$', 'interpreter', 'latex', 'fontsize',15)
    ylabel('$\frac{\| \hat{\Sigma}^{(n+1)} - \hat{\Sigma}^{(n)} \|}{\|\hat{\Sigma}^{(n)}\|}$', 'interpreter', 'latex', 'fontsize',20)
    grid
    grid minor
    temp = [sprintf('$\\rho_{0} = %.2f$, $\\alpha=%.2f$, $\\beta=%.2f$', rhoH0,alpha, beta)];
    textString = {sprintf('$p=%d$, $N=%d$, $T=%d$, $%d$ Trials',p,N,T,numberOfTrials), temp};
    annotation('textbox',...
    [0.25, 0.30, 0.1, 0.1], ...
    'Units','characters',...
    'String', textString,...
    'backgroundcolor', 'w', ...
    'FitBoxToText','on',...
    'interpreter', 'latex');

    % H1
    figure
    loglog(nitermax_vector, errMTH1, 'linestyle', '--', 'marker', 's')
    hold on
    loglog(nitermax_vector, errMGH1, 'linestyle', ':', 'marker', 'd')
    lg = {'$\hat{\Sigma}_{0}^{\mathrm{MT}}$', '$\hat{\Sigma}_{0}^{\mathrm{MG}}$'};
    for t=1:T
        loglog(nitermax_vector, errMCH1(:,t), 'linestyle', '-.', 'marker', 'o')
        lg = cat(2,lg, sprintf('$\\hat{\\Sigma}_{%d}^{\\mathrm{Tex}}$', t));
    end
    h = legend(lg);
    set(h, 'interpreter', 'latex', 'fontsize', 15)
    xlabel('$n$', 'interpreter', 'latex', 'fontsize',15)
    ylabel('$\frac{\| \hat{\Sigma}^{(n+1)} - \hat{\Sigma}^{(n)} \|}{\|\hat{\Sigma}^{(n)}\|}$', 'interpreter', 'latex', 'fontsize',20)
    grid
    grid minor
    temp = [sprintf('$\\rho_{1}=%.2f$', rhoH1(1))];
    for t=2:T
       temp = cat(2,temp, sprintf(', $\\rho_{%d}=%.2f$', t, rhoH1(t)));  
    end
    textString = {sprintf('$p=%d$, $N=%d$, $T=%d$, $%d$ Trials',p,N,T,numberOfTrials), temp, sprintf('$\\alpha=%.2f$, $\\beta=%.2f$',alpha, beta)};
    annotation('textbox',...
    [0.25, 0.30, 0.1, 0.1], ...
    'Units','characters',...
    'String', textString,...
    'backgroundcolor', 'w', ...
    'FitBoxToText','on',...
    'interpreter', 'latex');



%     autoArrangeFigures()
%% ------------- END CODE ----------------
