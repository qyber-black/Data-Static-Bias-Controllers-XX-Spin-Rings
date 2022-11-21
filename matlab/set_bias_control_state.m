function state = set_bias_control_state(id)
  % set_bias_control_state (states)
  %
  % id  Readout id: 0 - instantaneous readout
  %                 1 - .1 readout window
  %                 2 - .1 readout window with .05 perturbation and 500 samples
  %
  % Returns parameters for bias control for different cases/states used
  % by other functions to fix numerical experiment parameters.
  %

  if ~exist('id','var')
    error('No state given');
  end

  state.Nmax = 30;          % Max ring size

  state.T = NaN;            % Fix time?
  state.B = NaN;            % Maximum bias?
  state.initT = [500 5000]; % (1) - max time,
                            % (2) - samples for init times from chain (or empty)
  state.repeats = 2000;     % How many restarts?
  state.symm = 1;           % Enforce symmetry
  state.bias_init = 1;      % Initial biases with peak/trough/constant instead of random?
  state.noise = 0;          % Add noise to initial values (or 0)
  state.min_err = 0.001;    % Largest error acceptable for fastest solution.
  state.readout = [];       % Readout window (or empty to ignore)

  state.id_str = 'unknown';
  state.readout = [];

  if id == 0
    state.id_str = 't';
  elseif id == 1
    state.readout = [.1 0 0];      % Readout window (or empty to ignore)
                                   % readout(1) Time window for readout (0 - do not consider window)
                                   % readout(2) Perturbation (or 0 for none, time window required)
                                   % readout(3) Samples for expectation over pertubation (if pertubation not 0)
    state.min_err = 0.01;          % Largest error acceptable for fastest solution.
    state.id_str = 'dt';
  elseif id == 2
    state.readout = [.1 0.05 500]; % Readout window (or empty to ignore)
                                   % readout(1) Time window for readout (0 - do not consider window)
                                   % readout(2) Perturbation (or 0 for none, time window required)
                                   % readout(3) Samples for expectation over pertubation (if pertubation not 0)
    state.min_err = 0.1;           % Largest error acceptable for fastest solution.
    state.noise = 200;             % Add noise to initial values (or 0)
    state.id_str = 'dtp';
  else
    error(sprintf('Unknown state %d',state));
  end

end
