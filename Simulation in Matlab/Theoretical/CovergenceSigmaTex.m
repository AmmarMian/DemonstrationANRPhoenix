%% Description
% Test the convergence properties of Sigma_t^Tex
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
    p = 50; % Size of data vector
    N = 50;
    T = 3; % Number of dates
    numberOfTrials = 60;
    
    % Parameters
    rho = 0.3;
    Sigma = toeplitz(rho.^(0:p-1));
    
    % Fixed-point params
    tol = 0;
    nitermax_vector = floor(logspace(0,4,20));

%% Compute Criterion

    errTex = zeros(length(nitermax_vector), numberOfTrials);
    h = waitbar(0,'Computing Error...');
    for iter=1:length(nitermax_vector) % Loopon the number of observations
        nitermax = nitermax_vector(iter);
        
        errTextemp = zeros(numberOfTrials,1);
        parfor trial=1:numberOfTrials % M-C Trial
                       
            % We generate data using the same matrices
            xkt = GenerateGaussianData(p,N,Sigma);


            % We test te convergence of the fixed-point algorithm when we
            % estimate just one and we have a constant in the denominator
            err = inf;
            a = abs(randn); % Generating constant
            
            niter=1;
            R = eye(p);
            while (err > tol) && niter<=nitermax
                v = chol(R)' \ xkt;
                v = mean(v.*conj(v)); 
                aone =  v(ones(p,1),:); 
                all = a + aone;

                y = xkt ./ sqrt(all);
                Rnew = y*y'/N; % Numerator
                Rnew = p * Rnew / sum(diag(Rnew)); % Normalize by the trace
                err = norm(Rnew-R,'fro')/norm(R,'fro'); % Criterion
                R = Rnew;
                niter = niter+1;
            end 
            errTextemp(trial) = err;
            
           
        end
        errTex(iter, :) = errTextemp;
        waitbar(iter/length(nitermax_vector),h);
    end
    close(h);
    
%% Plotting

    % H0
    figure
    loglog(nitermax_vector, mean(errTex,2), 'linestyle', '--', 'marker', 's')
    hold on
    xlabel('$n$', 'interpreter', 'latex', 'fontsize',15)
    ylabel('$\frac{\| \hat{\Sigma}^{(n+1)} - \hat{\Sigma}^{(n)} \|}{\|\hat{\Sigma}^{(n)}\|}$', 'interpreter', 'latex', 'fontsize',20)
    grid
    grid minor
    
    
    figure
    boxplot(errTex.', 'labels', nitermax_vector)
    set(gca, 'xscale', 'log')
    set(gca, 'yscale', 'log')
    hold on
    xlabel('$n$', 'interpreter', 'latex', 'fontsize',15)
    ylabel('$\frac{\| \hat{\Sigma}^{(n+1)} - \hat{\Sigma}^{(n)} \|}{\|\hat{\Sigma}^{(n)}\|}$', 'interpreter', 'latex', 'fontsize',20)
    grid
    grid minor
%% ------------- END CODE ----------------
