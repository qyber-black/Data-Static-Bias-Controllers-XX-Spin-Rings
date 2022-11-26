% script to plot log sensitivity data

% SPDX-FileCopyrightText: Copyright (C) 2020-2022 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2020-2022 SM Shermer <lw1660@gmail.com>
% SPDX-FileCopyrightText: Copyright (C) 2022 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-License-Identifier: CC-BY-SA-4.0 

FILES = dir('../results/data_bias_control_robustness/*.mat')
for k=1:length(FILES)
    load(fullfile(FILES(k).folder,FILES(k).name)); 
    figure(1), clf, 
    plot_logSens_vs_error(robustness), drawnow
    title(FILES(k).name)
    savefig(1,sprintf('../figures/log-sensitivity/%s',FILES(k).name(1:end-4))); 
    saveas(1,sprintf('../figures/log-sensitivity/%s',FILES(k).name(1:end-4)),'png'); 
end
