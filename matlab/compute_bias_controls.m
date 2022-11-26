function compute_bias_control(states)
% compute_bias_control (states)
%
% states  Readout id: 0 - instantaneous readout
%                     1 - .1 readout window
%                     2 - .1 readout window with .05 perturbation and 500 samples
%         Can be a vector for multiple options
%
% Compute static bias control figures for information transfer from spin 1
% to  spin targets in N rings for N=3:30 and targets=2:floor((N+1)/2) and
% create summary figure of shortest times for all problems.

% SPDX-FileCopyrightText: Copyright (C) 2020-2022 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2020-2022 SM Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0 
  
  if ~exist('../data/data_bias_control','dir')
    mkdir('../data/data_bias_control');
  end

  if ~exist('states','var')
    states = [0 1 2];
  end

  problems = [];
  for id = states
    s = set_bias_control_state(id);
    for N = 3:s.Nmax % Ring size
      for target = 2:floor((N+1)/2) % 1 to target transition (in ring)
        problems = [problems; id, N, target];
      end
    end
  end
  M = size(problems, 1);
  parfor l = 1:M
    id = problems(l,1);
    N = problems(l,2);
    target = problems(l,3);
    calc_bias_controls (id,N,target,size(find(problems(:,1) == id),1));
  end

end
