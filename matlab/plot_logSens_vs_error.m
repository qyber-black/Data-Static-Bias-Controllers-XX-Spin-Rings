function plot_logSens_vs_error(robustness)
% plot the log-sensitivity vs error based on robustness 
% input: robustness -- data structure generated by calc_sensitivity
% output: plot

% SPDX-FileCopyrightText: Copyright (C) 2020-2022 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2020-2022 SM Shermer <lw1660@gmail.com>
% SPDX-FileCopyrightText: Copyright (C) 2022 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-License-Identifier: CC-BY-SA-4.0 

semilogy(robustness.error,robustness.logSensitivity.Controller,'r.',robustness.error,robustness.logSensitivity.Hamiltonian,'b.'), 
xlabel('error (infidelity)'), ylabel('log sensitivity'), legend('controller','couplings')
grid on, xlim([0 0.1])
