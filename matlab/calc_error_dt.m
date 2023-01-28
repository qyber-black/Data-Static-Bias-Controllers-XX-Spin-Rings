function err = calc_error_dt (H,Result,init,target,dt)
%
% suggested use
% R10 = qsn.QSN('ring',ones(1,10))
% err1 = cellfun(@(x)calc_error_dt(R10.H,x,1,3,0.1),Results)
% err2 = cellfun(@(x)calc_error_dt(R10.H,x,1,3,0.1*pi),Results)
% [~,ind]=sort(err1)
% plot(err1(ind)), hold on, plot(err2(ind)),
% set(gca,'yscale','log'), legend('err dt=0.1','err dt=0.1\pi')

% SPDX-FileCopyrightText: Copyright (C) 2022 SM Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0

T = Result.time;
H = H + diag(Result.bias);
N = size(H,1);
out = zeros(N,1);
in  = zeros(N,1);
out(target)= 1;
in(init)   = 1;
[V,E] = eig(H);
E = diag(E);
L = ones(N,1)*E.'-E*ones(1,N);
W = (out'*V).*(V'*in)';
err = 1-real(W*(exp(1i*L*T).*sinc(L*dt/(2*pi)))*W');
