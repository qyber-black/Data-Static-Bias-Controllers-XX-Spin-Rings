# Static Bias Controllers for XX Spin-1/2 Rings

> SPDX-FileCopyrightText: Copyright (C) 2020-2022 SM Shermer <lw1660@gmail.com>\
> SPDX-FileCopyrightText: Copyright (C) 2020-2022 Frank C Langbein <frank@langbein.org>, Cardiff University\
> SPDX-FileCopyrightText: Copyright (C) 2022 Sean Patrick O'Neil <seanonei@usc.edu>\
> SPDX-License-Identifier: CC-BY-SA-4.0

Propagation of information encoded in spin degrees of freedom through networks of coupled spins enables important applications in spintronics and quantum information processing. This data set contains the results of applying optimal control of information propagation in networks of spin-1/2 particles with uniform nearest neighbor XX-couplings forming a ring with a single excitation in the network as simple prototype of a router for spin-based information. The control is implemented via spatially distributed potentials, which remain constant during information transfer. The limited degrees of freedom makes finding a control that maximizes the transfer probability in a short time difficult.

It had been originally computed for

[1] FC Langbein, SG Schirmer, EA Jonckheere. **Time optimal information transfer in spintronics networks**. _IEEE Conf Decision and Control_, pp. 6454-6459, 2015.
[[DOI:10.1109/CDC.2015.7403236]](http://dx.doi.org/10.1109/CDC.2015.7403236)
[[arXiv:1508.00928]](http://arxiv.org/abs/1508.00928)
[[PDF:paper]](https://langbein.org/wp-content/uploads/2015/11/CDC2015.pdf)
[[Details]](https://langbein.org/langbein2015/)

For further details see

[2] SG Schirmer, EA Jonckheere, FC Langbein. **Design of Feedback Control Laws for Information Transfer in Spintronics Networks**. _IEEE Trans Automatic Control_, **63**(8):2523-2536, 2018.
[[DOI:10.1109/TAC.2017.2777187]](https://dx.doi.org/10.1109/TAC.2017.2777187)
[[arXiv:1607.05294]](http://arxiv.org/abs/1607.05294)
[[PDF:paper]](https://langbein.org/wp-content/uploads/2016/07/FeedbackControlLaws-2.pdf)
[[Details]](https://langbein.org/design-feedback-control-laws-information-transfer-spintronics-networks/)

[3] EA Jonckheere, SG Schirmer, FC Langbein. **Jonckheere-Terpstra test for nonclassical error versus log-sensitivity relationship of quantum spin network controllers**. _Int J Robust and Nonlinear Control_, **28**(6):2383-2403, 2018.
[[DOI:10.1002/rnc.4022]](https://dx.doi.org/10.1002/rnc.4022)
[[arXiv:1612.02784]](http://arxiv.org/abs/1612.02784)
[[PDF:paper]](https://langbein.org/wp-content/uploads/2016/10/spin-transport-1.pdf)
[[Details]](https://langbein.org/jonckheere-terpstra/)

This repository is mirrored at [https://github.com/qyber-black/Data-Static-Bias-Controllers-XX-Spin-Rings](https://github.com/qyber-black/Data-Static-Bias-Controllers-XX-Spin-Rings).

## Versions 

**Version 2.0.0**

S O'Neil, EA Jonckheere, FC Langbein, SG Schirmer. **Static Bias Controllers for XX Spin-1/2 Rings**, Version 2.0.0. Data set, 2022.
[Github: [https://github.com/qyber-black/Data-Static-Bias-Controllers-XX-Spin-Rings/releases/tag/v2.0.0](https://github.com/qyber-black/Data-Static-Bias-Controllers-XX-Spin-Rings/releases/tag/v2.0.0)]
[Qyber\black: [https://qyber.black/spinnet/data-static-bias-controllers-xx-spin-rings/-/tree/v2.0.0](https://qyber.black/spinnet/data-static-bias-controllers-xx-spin-rings/-/tree/v2.0.0)]

**Version 1.0.0**

FC Langbein, SG Schirmer, EA Jonckheere. **Static Bias Controllers for XX Spin-1/2 Rings**, Version 1.0.0. Data set, 2022.
[Github: [https://github.com/qyber-black/Data-Static-Bias-Controllers-XX-Spin-Rings/releases/tag/v1.0.0](https://github.com/qyber-black/Data-Static-Bias-Controllers-XX-Spin-Rings/releases/tag/v1.0.0)]
[Qyber\black: [https://qyber.black/spinnet/data-static-bias-controllers-xx-spin-rings/-/tree/v1.0.0](https://qyber.black/spinnet/data-static-bias-controllers-xx-spin-rings/-/tree/v1.0.0)]

The **original version** has been published as
  
FC Langbein, SG Schirmer, EA Jonckheere. **Static Bias Controllers for XX Spin-1/2 Rings**. Data set, figshare, 3rd July 2016.
[[DOI:10.6084/m9.figshare.3485240.v1]](https://dx.doi.org/10.6084/m9.figshare.3485240.v1)
[[Details]](https://langbein.org/static-bias-controllers-xx-spin-12-rings)

and is also available as [[FeedbackControlLaws20160703.zip]](https://d.qyber.black/data/data-static-bias-controllers-xx-spin-rings/FeedbackControlLaws20160703.zip).

## Data files

The controller dataset is contained in the directory `data`.  A summary of the controller dataset containing the core information is available as a csv file `data/data_cup.csv`.  The controller data files including details about the optimization parameters and additional information are available in Matlab `mat` v7.3 format (created with Matlab 2021b) in `data/data_bias_control`.  The naming convention is

    data_bias_control_D-N-M.mat for transfer
    data_localization_N.mat for localization on spin 1

where

    D - indicates type of readout at target time:
      t  - controllers for readout at specific time,
      dt - controllers for average readout over time window
      dtp - controllers subject to readout time uncertainty
    N - number of spins in the network from 3 to 20+;
    M - information transfer from spin |1> to |M> for M = 2,...,ceil(N/2).

Each data file contains the following variables/structures. `Info` is a struct containing information about the optimization options.  The information of interest is mostly in `Info.args` which has subfields:

    obj - matspinnet object (if available)
    in  - input spin (usually 1 for ring)
    out - output spin (between 1 and N)
    T   - constraints on transfer time (upper limit), NaN if none
    B   - constraints on bias field strength (upper limit), NaN if none
    readout - vector [0.1000 0 0]
    min_err - optional error threshold such as 0.0100
    repeats - number optimization repeats
    symm    - 1 to impose symmetry constraints, 0 no constraints
    initT   - initial range of transfer times, e.g., [500 5000]
    bias_init - initial bias values
    noise   - 0 if none added

`Results` is a cell array of optimisation results, with each cell containing the following fields:

    err       - error/infidelity of controller
    exec_time - time it took to find controller with L-BFGS optimiser, on single core Xeon E7 2.7GHz
    exec_flag - exist flag of matlab optimiser
    output    - output of matlab optimiser
    bias      - vector of biases on spins
    init_bias - initial value for bias
    time      - readout time at spin M (with time window time+/-info.dt/2)
    init_time - initial value for time
 
In addition there are two variables

    best    - index for results of best controller
    fastest - index for results of fastest controller (with error smaller than info.fastest_min_err).

## Figures

The data is also visualised in the following figures in PNG and Matlab figure files (generated with Matlab 2021b)

    bias_D-N-M.{fig,png}

These figures reside in the `/figures/core` directory. These depict results for optimizing the information transfer probability from spin `1` to `M` for a ring of `N` spins with XXZ coupling in the single excitation subspace over the spatial biases and time.  The left column shows the biases and evolution (in blue vs. the natural evolution in red) giving the best fidelity at time `T` and the best error. The middle columns shows the fastest solution found with a fidelity greater than fastest_min at time `T`. The right column shows the overall solutions found by repeated optimization, plotting the time vs. the logarithm of the infidelity and a histogram of the logarithm of the infidelity. The bottom row shows the eigenstructure of the best and the fastest solution and their symmetries, with the eigenvectors being the columns of the matrices (in cyan, green and red rows indicate the `|in>` and `|out>` state resp.) and the corresponding eigenvalues at the bottom (in purple).
    
    fastest_D.{fig.png}

Also located in the `/figures/core directory`, these display the shortest times achieved for instantaneous transition fidelities greater than `0.999` for `D=t` or `0.99` for `D=dt` for rings of size `N = 3, ..., 20` and transitions from `1` to `M = 2, ..., ceil(N/2)`. Note that for some transitions no solution with high fidelity were found, so no fastest results are reported (see data files).  The color of the bars indicate the infidelity of the fastest solution.

    sensitivity_D-N-M.{fig,png}

Logarithm of transfer probability (red) and logarithm of sensitivity (blue), ordered by increasing infidelity from left to right, of the `1 -> M` controllers of an `N`-ring. These figures are located in the `figures/sensitivity` directory. 

    data_bias_control_D-N-M-robustness.{fig,png} or data_localization_D-M-robustness.{fig,png}

These figures are located in the `figures/log-sensitivity` directory and display a plot of the error versus log-sensitivity to Hamiltonian and controller perturbations for each of the controllers included in the `data/data_bias_control` directory.

    log_sens_figure_D-N-M_S=P.fig

These figures are located in the `figures/log-sensitivity-updated` figures directory and depict log-sensitivity versus error for the entire controller set in a given transfer `1 -> M` and for a particular perturbation `P`. For `P` less than or equal `N`, the perturbations are to the bias fields. For `P` between `N+1` and `2N`, the  perturbations are to the spin couplings (Hamiltonian perturbations). 

    log_sens_composite_D-N-M.{fig,png}

Located in the `figures/log-sensitivity-updated/composite` directory, these figures provide a log-log plot of the log-sensitivity of the fidelty error to Hamiltonian perturbations and controller perturbations versus the fidelity error. The data in these plots is filtered for controllers with a fidelity greater than `0.9`. The blue crosses depict the norm of the `N` perturbations to the couplings, and the red asteriks depcit the norm of the `N` bias-field pertrubations. 

## Dataset generation

The dataset was generated using the `create_data_set.m` function in the the directory `matlab`.  This function basically converts the controller data in the directory matlab/results.  This data was generated using https://qyber.black/spinnet/code-matspinnet and the matlab functions in the directory matlab.

## Analysis

The following routines are available in the `matlab` directory to calculate (log-)sensitivity, analyze and visualize the data.

    script_calc_robustness_all.m 

This is the main routine for calculating robustness data (sensitivity and log-sensitivity) from the 
controller data files in the `data/data_bias_control` directory and in accordance with [3]. If the
`data/data_bias_control directory` is not empty, the routine requires no input. It calls on 
`calc_sensitivity.m` for the computations and saves the files as `data_bias_control_D-N-M-robustness` 
and `data_localization_D-N-robustness` in the `results/data_bias_control_robustness` directory.

    calc_sensitivity.m 

This routine calculates the sensitivity and log-sensitivity to both Hamiltonian and
controller perturbations. It takes as input the data file provided by `script_calc_robustness_all.m` and calls on the `bias_sensitivity.m` routine from the `MatSpinNet` package to calculate the sensitivity in accordance with [3]. It then uses this data to compute the log-sensitivity. It saves all computed robustness data in a variable ”robustness” and sends to `script_calc_robustness_all.m` for saving.
    
    script_plot_log_sensitivity_all.m

This script is the main routine for plotting of robustness data from the results in the 
`results/data_bias_control_robustness` directory. It requires no input and is designed to run from 
the `/matlab directory`. It calls the `plot_logSens_vs_error.m` routine to generate plots from the 
data in `results/data_bias_control_robustness` directory. It saves the plots as MATLAB figures and 
accompanying .png files in the `figures/log-sensitivity` directory with filenames
`data_bias_control_D-N-M-robustness` or `data_localization_D-N-robustness`.

    plot_logSen_vs_error.m 

This script produces a scatter plot of log-sensitivity to Hamiltonian and controller perturbations 
versus error based on the data files in the `results/data_bias_control_robustness` directory. The 
script is designed to run based on a call from the parent `script_plot_log_sensitivity_all.m` function.

    calc_log_sens.m 
    
This is the main routine for calculation of the log-sensitivity of the error to Hamiltonian and bias field perturbations using an updated analytical formula of the Bloch-transformed system. It requires no user input but takes the controller and transfer time data from the `data_bias_control_D-N-M` or `data_localisation_D-N` data files and computes the fidelity error, sensitivity of the error to parameter variation, and log-sensitivity of the error to parameter variations On output the routine saves the data as `log_sens_data_D-N-M` in the `results/log_sens_results` directory. It also creates and saves `.fig` files for the log-sensitivity versus error in the `figures/log-sensitivity-updated/figures` directory. Finally, it also creates `.fig` and `.png` plots of the log-sensitivity versus fidelity error for high-fideity controllers (those with and error less than 0.1) and saves them in the `figures/log-sensitivity-updated/composite` directory.