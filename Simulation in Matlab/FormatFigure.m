%% Description
% Formatting figures in Matlab for easy integration in Latex document
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
function FormatFigure(handle, titleString, xlabelString, ylabelString, legendString)
% FormatFigure: Formatting figures in Matlab for easy integration in Latex document
% Syntax:  FormatFigure(handle)
% Inputs:
%   handle - handle of the figure
% Outputs: none

    if nargin < 5
        legendString = '';
    end
    if nargin < 4
        ylabelString = '';
    end
    if nargin < 3
        xlabelString = '';
    end
    if nargin < 2
       titleString = '';
    end

    % Title
    title(titleString)
    % X,Y Label
    xlabel(xlabelString)
    ylabel(ylabelString)
    % Legend
    if ~strcmp(legendString, '')
        legend(legendString)
    end
    % White Background
    set(handle, 'color', 'white')
    % Interpreter
    set(findall(handle,'-property','Interpreter'),'Interpreter','Latex');
    % Font
    set(findall(handle,'-property','FontName'),'FontName','Times');
    % Grid
    grid off
    grid on
    grid minor

end
%% ------------- END CODE ----------------