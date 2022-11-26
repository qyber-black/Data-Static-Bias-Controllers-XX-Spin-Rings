% Calculate bias controls for id (t/dt/dtp) for ring of size N from 1 to target. 
% If M is zero, always calculate. Otherwise M is number of completed tasks and 
% only calculate if results do not yet exist.

function calc_bias_controls(id,N,target,M)

% SPDX-FileCopyrightText: Copyright (C) 2020-2022 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2020-2022 SM Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0 

  ring = qsn.QSN ('ring', N);
  s = set_bias_control_state(id);
  name = sprintf('bias_control_%s-%d-%d', s.id_str, N, target);
  % Find bias controls
  if M == 0 || ~exist(['../data/data_bias_control/data_' name '.mat'], 'file')
    parsave(name, [], [], [], [], [], sprintf('bias_control_%s',s.id_str), M);
    if M == 0
      v = 1;
    else
      v = 0;
    end
    [Results,best,fastest,FailedRuns,Info] = ring.bias_control(1, target, s.T, s.B, s.readout, s.min_err, s.repeats, s.symm, s.initT, s.bias_init, s.noise, v);
    parsave(name, Results, best, fastest, FailedRuns, Info, sprintf('bias_control_%s',s.id_str), M);
  else
    disp(sprintf('%s: OK', name));
  end
end
