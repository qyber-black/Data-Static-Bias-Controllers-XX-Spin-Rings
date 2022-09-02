# Static Bias Controllers for XX Spin-1/2 Rings

The data set contains the original energy landscape controller dataset used in several papers including
https://qyber.black/spinnet/paper-feedback-control-laws and
https://qyber.black/spinnet/paper-sensiitivity-jt.

## Data files

The controller dataset is contained in the directory data-set.  The data files are in Matlab MAT v7.3 format (created with Matlab 2021b).  The naming convention is

  <code>data_D-N-M.mat</code>

where

  <code>D</code> - indicates type of readout at target time:
      <code>t </code>  controllers for readout at specific time,
      <code>dt</code>  controllers for average readout over time window;
  <code>N</code> - number of spins in the network from 3 to 20;
  <code>M</code> - information transfer from spin <code>|1></code> to <code>|M></code> for <code>M = 2,...,ceil(N/2)</code>.

Each data file contains the following variables/structures

  <code>info</code> - struct containing data set information

    <code>N </code>              - size of XX ring
    <code>in</code>              - initial state / input spin
    <code>out</code>             - readout state / output spin
    <code>fastest_min_err</code> - maximum error for optimisation result to be considered for fastest information transfer
    <code>dt</code>              - readout time window; 0 for instantaneous readout

  <code>results</code> - cell array of optimisation results, with each cell containing the following fields:

    <code>err</code>       - error/infidelity of controller
    <code>exec_time</code> - time it took to find controller with L-BFGS optimiser, on single core Xeon E7 2.7GHz
    <code>exec_flag</code> - exist flag of matlab optimiser
    <code>output</code>    - output of matlab optimiser
    <code>bias</code>      - vector of biases on spins
    <code>init_bias</code> - initial value for bias
    <code>time</code>      - readout time at spin M (with time window time+/-info.dt/2)
    <code>init_time</code> - initial value for time
    <code>results_idx_best</code>    - index for results of best controller
    <code>results_idx_fastest</code> - index for results of fastest controller (with error smaller than info.fastest_min_err).

 <code>sensitivity</code> - struct of controller sensitivities

    <code>error</code>     - errors of controllers, in order of incrasing infidelity
    <code>dpdJ_norm</code> - norm of sensitivities of controller w.r.t. uncertainties in couplings in same order than error.
    <code>taub</code>      - Kendall \tau_b correlation between error and sensitivity

## Figures

The data is also visualised in the following figures in PNG and Matlab figure files (generated with Matlab 2021b)

* <code>bias_D-N-M.{fig,png}</code>

Results for optimizing the information transfer probability from spin 1 to M for an XX-ring of N spins over the spatial biases and time.  The left column shows the biases and evolution (in blue vs. the natural evolution in red) giving the best fidelity at time T and the best error. The middle columns shows the fastest solution found with a fidelity greater than fastest_min at tme T The right column shows the overall solutions found by repeated optimization, plotting the time vs. the logarithm of the infidelity and a histogram of the logarithm of the infidelity. The bottom row shows the eigenstructure of the best and the fastest solution and their symmetries, with the eigenvectors being the columns of the matrices (in cyan, green and red rows indicate the |in> and |out> state resp.) and the corresponding eigenvalues at the bottom (in purple).

*  <code>fastest_D.{fig.png}</code>

Shortest times achieved for instantaneous transition fidelities greater than 0.999 for D=t or 0.99 for D=dt for rings of size N = 3, ..., 20 and transitions from 1 to M = 2, ..., ceil(N/2). Note that for some transitions no solution with high fidelity were found, so no fastest results are reported (see data files).  The color of the bars indicate the infidelity of the fastest solution.

*  sensitivity_D-N-M.{fig,png}

Logarithm of transfer probability (red) and logarithm of sensitivity (blue), ordered by increasing infidelity from left to right, of the 1 -> M controllers of an N-ring.

## Dataset generation

The dataset was generated using the <code>create_data_set.m</code> function in the the directory <code>matlab</matlab>.  This function basically converts the controller data in the directory <code>matlab/results</code>.  This data was generated using https://qyber.black/spinnet/code-matspinnet and the matlab functions in the directory <code>matlab</code>.

%% Analysis

The directory <code>analysis_RNC</code> contains some code previously used to analyze the dataset in https://qyber.black/spinnet/paper-sensiitivity-jt and spreadsheets summarizing the previous results.
