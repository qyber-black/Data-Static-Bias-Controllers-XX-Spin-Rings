% SPDX-FileCopyrightText: Copyright (C) 2020-2022 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2020-2022 SM Shermer <lw1660@gmail.com>
% SPDX-FileCopyrightText: Copyright (C) 2022 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-License-Identifier: CC-BY-SA-4.0 

load results/data_bias_control_t-10-3.mat
R10 = qsn.QSN('ring',ones(1,10))
for k=1:100
  sys = Results{k};
  U = expm(1i*sys.time*(R10.H+diag(sys.bias)));
  err0(k) = sys.err;
  err1(k) = 1- abs(U(1,3))^2;
  err2(k) = 1- abs(U(1,3));
end