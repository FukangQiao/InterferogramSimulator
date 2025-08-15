% ===========================================================
% Filename:     main.m
% Date:   	 	2025-08-14
% Author:    	Fukang Qiao
% Description:  Main function for data phase simulation
% ===========================================================

%% Parameters
close all; clear; clc;
addpath('WZPUtil');

%%% Input Parameters
par.savePath = '../data/qfk';  % Data Save Path 
par.totalNum = 10;           % Total Numbers
par.sampleSize = 256;        % Size of Samples
par.multilook = [1,4];       % multilook [Azimuth, Range] [Row, Col]
par.demFolder = '../DEM/unzip'; % Folder containing multiple DEM (*.tif format)

%%% Output Parameters (Can be customize other types and implement in simulateData.m function
par.savePNGFlag = 1 ;      % Save the *.png image (1 means yes/0 means no)
par.out.origin = 1;        % Original phase
par.out.originWrapped = 1; % Wrapped of the original phase without Noise
par.out.interf = 1;        % Simulated interferogram with Noise
par.out.coherence = 1;     % Estimated coherence
par.out.deformBbox = 1;    % Deformation area Location and Label
par.out.VDRI = 1;          % Horizontal phase gradient + Residual map + Vertical phase gradient + Interferogram (4 channels)
par.out.branchCut = 1;     % Horizontal branch-cut and vertical branch-cut (2 channels)

%%% Phase Components with a certain probability(0-1) eg.CHANGE IT -1 means no phase
par.probSlop = 0.2;          % Slope phase
par.probBuilding = 0.2;      % Building phase
par.probTurbulence = 0.8;    % Atmospheric turbulence phase, i.e., fractal Perlin noises
par.probDeform = 0.8;        % Distorted two-dimensional Gaussian surface
par.probEarthquake = -1;     % Deformation caused by earthquakes
par.probWater = 0.2;         % Completely decorrelated area

par.noiseType = 0;           % 0 : Deformation-related noise (recommend)
                             % 1 : Noise with random signal-to-noise ratio(not recommended)
par.noiseSNRRange = [0.2,5]; % Need to specify when noiseType=1


%%% Parallel
par.Parallel = 0;           % 0 means no parrallel

%% Generate
generate(par);

%% Show Examples
showSamples(par);

