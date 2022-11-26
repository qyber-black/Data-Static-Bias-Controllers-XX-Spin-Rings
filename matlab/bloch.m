%% Produce Bloch Transformed Matrices
function [A,Sb,lambda,V,r_in,r_out,c] = bloch(N,Hd,S,e,out_index,b)

% SPDX-FileCopyrightText: Copyright (C) 2020-2022 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2020-2022 SM Shermer <lw1660@gmail.com>
% SPDX-FileCopyrightText: Copyright (C) 2022 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-License-Identifier: CC-BY-SA-4.0 

% Define input and output 
in = e{1};
out = e{out_index}';
rho_in = e{1}*e{1}';
rho_out = e{out_index}*e{out_index}';

% Produce Bloch Representations of H (A)
for n=1:N^2
    for m=1:N^2
        A(n,m)=trace(1*i*Hd*(b{n}*b{m}-b{m}*b{n}));
    end
    r_in(n,1) = trace(b{n}*rho_in);
    r_out(n,1) = trace(b{n}*rho_out);
end

% Produce Bloch Representations of S (Sb)
for run = 1:2*N
    for n=1:N^2
        for m=1:N^2
            Sb{run}(n,m)=trace(1*i*S{run}*(b{n}*b{m}-b{m}*b{n}));
        end
    end
end

% Set c to allow e(t) = c*exp(A*t)*r_in
c(1,1:N^2-1) = -r_out(1:N^2-1,1);
c(1,N^2) = (N-1)/sqrt(N);

% Get eigvalues and eigenvectors of A
[V,L] = eig(A);
for k = 1:N^2
    lambda(k,1) = L(k,k);
end
