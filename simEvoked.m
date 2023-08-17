%% File for simulating Evoked potential data 
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
% Code initially simulated CCEP data
% DH edited to remove stimulation artifact
%

%% Start by creating simulated parameters and channels


set(0, 'DefaultFigureRenderer', 'painters');

srate = 500;
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



