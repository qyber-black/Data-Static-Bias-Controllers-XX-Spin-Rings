function robustness = calc_sensitivity(file)
%
% input: file -- filename of matlab controller results, e.g.,
%        file = 'results/data_bias_control_dt-10-2.mat'
% output: robustness -- structure with fields
%         Error -- optimization error (here 1-fidelity)
%         Sensitivity -- vector of sensitivities for perturbations of biases and couplings
%         logSens_Hamiltonian -- average log sensitivity for coupling perturbations
%         logSens_Controller  -- average log sensitivity for controller perturbations
%
load(file)
N = length(Info.controls{1});
L = length(Results);
robustness.error = arrayfun(@(n)Results{n}.err,1:L);
robustness.sensitivity = Info.args.obj.bias_sensitivity(Results,Info);
robustness.logSensitivity.Hamiltonian = arrayfun(@(n)norm(robustness.sensitivity{n}(N+1:2*N))/Results{n}.err,1:L);
robustness.logSensitivity.Controller  = arrayfun(@(n)norm(robustness.sensitivity{n}(1:N))/Results{n}.err,1:L);