function robustness = calc_sensitivity(file)
%
% input: file -- filename of matlab controller results, e.g.,
%        file = 'results/data_bias_control_dt-10-2.mat'
%
% output: robustness -- structure with fields
%         Error -- optimization error (here 1-fidelity)
%         Sensitivity -- vector of sensitivities for perturbations of biases and couplings
%         logSens_Hamiltonian -- average log sensitivity for coupling perturbations
%         logSens_Controller  -- average log sensitivity for controller perturbations

% SPDX-FileCopyrightText: Copyright (C) 2020-2022 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2020-2022 SM Shermer <lw1660@gmail.com>
% SPDX-FileCopyrightText: Copyright (C) 2022 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-License-Identifier: CC-BY-SA-4.0 

load(file)
N = length(Info.controls{1});
L = length(Results);
robustness.error = arrayfun(@(n)Results{n}.err,1:L);
robustness.sensitivity = Info.args.obj.bias_sensitivity(Results,Info);

for x = 1:L
robustness.log_sensitivity{x} = robustness.sensitivity{x}./robustness.error(x);
end

robustness.logSensitivity.Hamiltonian = arrayfun(@(n)norm(robustness.log_sensitivity{n}(N+1:2*N)),1:L);
robustness.logSensitivity.Controller  = arrayfun(@(n)norm(robustness.log_sensitivity{n}(1:N)),1:L);

%robustness.logSensitivity.Hamiltonian = arrayfun(@(n)norm(robustness.sensitivity{n}(N+1:2*N))/Results{n}.err,1:L);
%robustness.logSensitivity.Controller  = arrayfun(@(n)norm(robustness.sensitivity{n}(1:N))/Results{n}.err,1:L);