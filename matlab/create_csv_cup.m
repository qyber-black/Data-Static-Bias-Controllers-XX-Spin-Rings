function out = create_csv_cup(mat_file_name,out_file_name)
% reads mat_file_name, which must be a matspinnet controller file, and 
% extracts the size of the spin network, input and output spins, 
% fidelity error, time of transfer, and bias applied to each spin for each
% controllers, and saves results as a .csv file out_file_name.csv

% SPDX-FileCopyrightText: Copyright (C) 2020-2022 SM Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0 

load(mat_file_name)
if ~exist('Results','var')
    error('no results variable containing controllers')
end
% 
for k=1:length(Results) 
    tmp = Results{k}; 
    if isempty(Info.args.readout)
        dt = 0;
    else
        dt = Info.args.readout(1);
    end
    out(k,:) = [length(tmp.bias), Info.args.in, Info.args.out, dt, tmp.err, tmp.time, tmp.bias']; % tmp.init_time, tmp.init_bias']; 
end
csvwrite(sprintf('../data/%s.csv',out_file_name),out);
disp(name)
