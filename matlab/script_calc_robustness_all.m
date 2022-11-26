% Script to generate robustness data (currently only logSensitivity) 
% for all controllers in data/data_bias_control directory

% SPDX-FileCopyrightText: Copyright (C) 2020-2022 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2020-2022 SM Shermer <lw1660@gmail.com>
% SPDX-FileCopyrightText: Copyright (C) 2022 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-License-Identifier: CC-BY-SA-4.0 

if ~exist('../results/data_bias_control_robustness','dir')
    mkdir('../results/data_bias_control_robustenss');
end

FILES = dir('../data/data_bias_control/*.mat')
for k=1:length(FILES)
    robustness = calc_sensitivity(fullfile(FILES(k).folder,FILES(k).name)); 
    save(sprintf('../results/data_bias_control_robustness/%s-robustness.mat',FILES(k).name(1:end-4)),'robustness'); 
end