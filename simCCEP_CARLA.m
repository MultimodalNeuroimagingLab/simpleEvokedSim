%% File for simulating CCEP data to use for evaluating CARLA
%
%    Copyright (C) 2023 Harvey Huang
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
% 
%% Start by creating simulated parameters and channels

addpath('functions');

set(0, 'DefaultFigureRenderer', 'painters');

srate = 4800;
tt = -0.5:1/srate:1-1/srate;

nchs = 50; % number of simulated channels
ntrs = 12; % number of trials
chs = arrayfun(@(x) sprintf('ch%d', x), 1:nchs, 'UniformOutput', false)';

%% A) Create data, add stim artifacts and EP

seed = 10; % controls the seed for all subsequent steps. 

rng(seed);

% stores the data
V0 = zeros(length(tt), nchs);

Aart = 50 + rand(nchs, 1)*5; % slightly different artifact amplitudes for each channel
artifact = sin(2*pi*600*tt)';
artifact(tt < 0 | tt > 0.002) = 0;
V0 = V0 + artifact*Aart';

chsResp = 1:20;

A = 100;
sig = genRandSig(tt, length(chsResp), A);

V1 = V0;
V1(:, chsResp) = V0(:, chsResp) + sig;

figure('Position', [200, 200, 400, 800]); plotTrials(tt, V1, 80, [], [], 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Channels');

%% B) Add a systematic noise across all chs and trials, e.g. from reference.

rng(seed);

Aglobal = 0; % same amplitude of noise on all channels

sigCommon = genRandSig(tt, 1, Aglobal);
V2 = V1 + sigCommon;

figure('Position', [200, 200, 400, 800]); plotTrials(tt, V2, 80, [], [], 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Channels');

%% C) Add noise common to all channels (line noise and some brown noise from reference)

rng(seed);

% C) Add common noise to all channels at each trial
V3 = repmat(V2', 1, 1, ntrs); % ch x time points x trial
phLN = rand(ntrs, 3)*2*pi; % LN phases
LN = zeros(length(tt), ntrs);
for ii = 1:ntrs
    LN(:, ii) = 8*sin(2*pi*60*tt - phLN(ii, 1)) + 2*sin(2*pi*120*tt - phLN(ii, 2)) + 1*sin(2*pi*180*tt - phLN(ii, 3));
end

% brown noise shared across channels
BN = cumsum(0.5*randn(2*length(tt), ntrs));
BN = ieeg_highpass(BN, srate, true);
BN = BN((0.5*length(tt)+1) : 1.5*length(tt), :);

noiseCommon = LN + BN;
V3 = V3 + shiftdim(noiseCommon, -1);

figure('Position', [200, 200, 800, 800]);
subplot(1, 2, 1); plotTrials(tt, V3(:, :, 1)', 80, [], [], 'LineWidth', 1.5);
title('V at one trial');
xlabel('Time (s)'); ylabel('Channels');
subplot(1, 2, 2); plotTrials(tt, mean(V3, 3)', 80, [], [], 'LineWidth', 1.5);
title('V across trials');
xlabel('Time (s)'); ylabel('Channels');

%% D) add random brown noise across all channels

rng(seed);

noiseRand = cumsum(0.5*randn(nchs, 2*length(tt), ntrs), 2); % give double the number of time points so we can highpass it
for ii = 1:nchs
    noiseRand(ii, :, :) = ieeg_highpass(squeeze(noiseRand(ii, :, :)), srate, true);
end
noiseRand = noiseRand(:, (0.5*length(tt)+1) : 1.5*length(tt), :);
V4 = V3 + noiseRand;

figure('Position', [200, 200, 800, 800]);
subplot(1, 2, 1); plotTrials(tt, V4(:, :, 1)', 80, [], [], 'LineWidth', 1.5);
title('V at one trial');
xlabel('Time (s)'); ylabel('Channels');
subplot(1, 2, 2); plotTrials(tt, mean(V4, 3)', 80, [], [], 'LineWidth', 1.5);
title('V across trials');
xlabel('Time (s)'); ylabel('Channels');

%% Get CAR by using the channels with lowest variance that yield least correlated structure across all channels

rng(seed);

opts.vartype = 'cov';
badChs = {10, []; 12, []; 25, 1:8};

[Vout, CAR, stats] = ccep_CARLA(tt, V4, srate, badChs, opts); % use sensitive method


% number of channels used for the CAR
nCAR = length(stats.chsUsed);
cmSens = [1, 165/255, 0];
[~, nCARglob] = max(mean(stats.zMinMean, 2)); % number of channels at global maximum

% Plot average zmin across trials
figure('Position', [200, 200, 400, 300]); hold on
errorbar(mean(stats.zMinMean, 2), std(stats.zMinMean, 0, 2), 'k-o');
plot(nCARglob, mean(stats.zMinMean(nCARglob, :), 2), 'b*'); % global max as blue
if nCARglob ~= nCAR; plot(nCAR, mean(stats.zMinMean(nCAR, :), 2), '*', 'Color', cmSens); end
yline(0, 'Color', 'k');

% Plot the variance of channels, sorted in increasing order
vars = stats.vars(stats.order);
figure('Position', [600, 200, 400, 300]); hold on
plot(vars, 'k-o', 'LineWidth', 1, 'MarkerFaceColor', 'k');
xline(nCARglob + 0.5, 'Color', 'b'); % global threshold
if nCARglob ~= nCAR; xline(nCAR + 0.5, 'Color', cmSens); end % sensitive threshold
xlim([0, nchs+1]);

% Sort  and plot channels by increasing covariance, draw line at cutoff
V4MeanSorted = mean(V4(stats.order, :, :), 3);

% create colormap where responsive channels are red
respBool = antifind(chsResp, nchs);
respBool = respBool(stats.order); % logical array of where responsive channels are
cm = zeros(nchs, 3);
cm(respBool, 1) = 1; % make red
figure('Position', [200, 200, 250, 600]);
yspace = 80;
ys = plotTrials(tt, V4MeanSorted', yspace, [], cm, 'LineWidth', 1);
yline(ys(nCARglob)-yspace/2, 'Color', 'b', 'LineWidth', 1.5);
if nCARglob ~= nCAR; yline(ys(nCAR)-yspace/2, 'Color', cmSens, 'LineWidth', 1.5); end
xlim([-0.1, 0.5]); set(gca, 'xtick', [0, 0.5]);
xlabel('Time (s)'); ylabel('Channels');

