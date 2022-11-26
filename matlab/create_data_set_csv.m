% This script generates a .csv file consolidating all the error, transfer time, 
% and controller data from the data sets from N=3 to N=20. The property of 
% each entry can be read from the table headers. C1 to C20 depict the bias 
% applied to spin 1 through 20, respectivley. 

% SPDX-FileCopyrightText: Copyright (C) 2020-2022 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2020-2022 SM Shermer <lw1660@gmail.com>
% SPDX-FileCopyrightText: Copyright (C) 2022 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-License-Identifier: CC-BY-SA-4.0 

function create_data_set_csv()

LIST = {
'data_localisation_dt-3'
'data_bias_control_t-3-2'
'data_bias_control_dt-3-2'
'data_localisation_dt-4'
'data_bias_control_t-4-2'
'data_bias_control_dt-4-2'
'data_localisation_dt-5'
'data_bias_control_t-5-2'
'data_bias_control_dt-5-2'
'data_bias_control_t-5-3'
'data_bias_control_dt-5-3'
'data_localisation_dt-6'
'data_bias_control_t-6-2'
'data_bias_control_dt-6-2'
'data_bias_control_t-6-3'
'data_bias_control_dt-6-3'
'data_localisation_dt-7'
'data_bias_control_t-7-2'
'data_bias_control_dt-7-2'
'data_bias_control_t-7-3'
'data_bias_control_dt-7-3'
'data_bias_control_t-7-4'
'data_bias_control_dt-7-4'
'data_localisation_dt-8'
'data_bias_control_t-8-2'
'data_bias_control_dt-8-2'
'data_bias_control_t-8-3'
'data_bias_control_dt-8-3'
'data_bias_control_t-8-4'
'data_bias_control_dt-8-4'
'data_localisation_dt-9'
'data_bias_control_t-9-2'
'data_bias_control_dt-9-2'
'data_bias_control_t-9-3'
'data_bias_control_dt-9-3'
'data_bias_control_t-9-4'
'data_bias_control_dt-9-4'
'data_bias_control_t-9-5'
'data_bias_control_dt-9-5'
'data_localisation_dt-10'
'data_bias_control_t-10-2'
'data_bias_control_dt-10-2'
'data_bias_control_t-10-3'
'data_bias_control_dt-10-3'
'data_bias_control_t-10-4'
'data_bias_control_dt-10-4'
'data_bias_control_t-10-5'
'data_bias_control_dt-10-5'
'data_localisation_dt-11'
'data_bias_control_t-11-2'
'data_bias_control_dt-11-2'
'data_bias_control_t-11-3'
'data_bias_control_dt-11-3'
'data_bias_control_t-11-4'
'data_bias_control_dt-11-4'
'data_bias_control_t-11-5'
'data_bias_control_dt-11-5'
'data_bias_control_t-11-6'
'data_bias_control_dt-11-6'
'data_localisation_dt-12'
'data_bias_control_t-12-2'
'data_bias_control_dt-12-2'
'data_bias_control_t-12-3'
'data_bias_control_dt-12-3'
'data_bias_control_t-12-4'
'data_bias_control_dt-12-4'
'data_bias_control_t-12-5'
'data_bias_control_dt-12-5'
'data_bias_control_t-12-6'
'data_bias_control_dt-12-6'
'data_localisation_dt-13'
'data_bias_control_t-13-2'
'data_bias_control_dt-13-2'
'data_bias_control_t-13-3'
'data_bias_control_dt-13-3'
'data_bias_control_t-13-4'
'data_bias_control_dt-13-4'
'data_bias_control_t-13-5'
'data_bias_control_dt-13-5'
'data_bias_control_t-13-6'
'data_bias_control_dt-13-6'
'data_bias_control_t-13-7'
'data_bias_control_dt-13-7'
'data_localisation_dt-14'
'data_bias_control_t-14-2'
'data_bias_control_dt-14-2'
'data_bias_control_t-14-3'
'data_bias_control_dt-14-3'
'data_bias_control_t-14-4'
'data_bias_control_dt-14-4'
'data_bias_control_t-14-5'
'data_bias_control_dt-14-5'
'data_bias_control_t-14-6'
'data_bias_control_dt-14-6'
'data_bias_control_t-14-7'
'data_bias_control_dt-14-7'
'data_localisation_dt-15'
'data_bias_control_t-15-2'
'data_bias_control_dt-15-2'
'data_bias_control_t-15-3'
'data_bias_control_dt-15-3'
'data_bias_control_t-15-4'
'data_bias_control_dt-15-4'
'data_bias_control_t-15-5'
'data_bias_control_dt-15-5'
'data_bias_control_t-15-6'
'data_bias_control_dt-15-6'
'data_bias_control_t-15-7'
'data_bias_control_dt-15-7'
'data_bias_control_t-15-8'
'data_bias_control_dt-15-8'
'data_localisation_dt-16'
'data_bias_control_t-16-2'
'data_bias_control_dt-16-2'
'data_bias_control_t-16-3'
'data_bias_control_dt-16-3'
'data_bias_control_t-16-4'
'data_bias_control_dt-16-4'
'data_bias_control_t-16-5'
'data_bias_control_dt-16-5'
'data_bias_control_t-16-6'
'data_bias_control_dt-16-6'
'data_bias_control_t-16-7'
'data_bias_control_dt-16-7'
'data_bias_control_t-16-8'
'data_bias_control_dt-16-8'
'data_localisation_dt-17'
'data_bias_control_t-17-2'
'data_bias_control_dt-17-2'
'data_bias_control_t-17-3'
'data_bias_control_dt-17-3'
'data_bias_control_t-17-4'
'data_bias_control_dt-17-4'
'data_bias_control_t-17-5'
'data_bias_control_dt-17-5'
'data_bias_control_t-17-6'
'data_bias_control_dt-17-6'
'data_bias_control_t-17-7'
'data_bias_control_dt-17-7'
'data_bias_control_t-17-8'
'data_bias_control_dt-17-8'
'data_bias_control_t-17-9'
'data_bias_control_dt-17-9'
'data_localisation_dt-18'
'data_bias_control_t-18-2'
'data_bias_control_dt-18-2'
'data_bias_control_t-18-3'
'data_bias_control_dt-18-3'
'data_bias_control_t-18-4'
'data_bias_control_dt-18-4'
'data_bias_control_t-18-5'
'data_bias_control_dt-18-5'
'data_bias_control_t-18-6'
'data_bias_control_dt-18-6'
'data_bias_control_t-18-7'
'data_bias_control_dt-18-7'
'data_bias_control_t-18-8'
'data_bias_control_dt-18-8'
'data_bias_control_t-18-9'
'data_bias_control_dt-18-9'
'data_localisation_dt-19'
'data_bias_control_t-19-2'
'data_bias_control_dt-19-2'
'data_bias_control_t-19-3'
'data_bias_control_dt-19-3'
'data_bias_control_t-19-4'
'data_bias_control_dt-19-4'
'data_bias_control_t-19-5'
'data_bias_control_dt-19-5'
'data_bias_control_t-19-6'
'data_bias_control_dt-19-6'
'data_bias_control_t-19-7'
'data_bias_control_dt-19-7'
'data_bias_control_t-19-8'
'data_bias_control_dt-19-8'
'data_bias_control_t-19-9'
'data_bias_control_dt-19-9'
'data_bias_control_t-19-10'
'data_bias_control_dt-19-10'
'data_localisation_dt-20'
'data_bias_control_t-20-2'
'data_bias_control_dt-20-2'
'data_bias_control_t-20-3'
'data_bias_control_dt-20-3'
'data_bias_control_t-20-4'
'data_bias_control_dt-20-4'
'data_bias_control_t-20-5'
'data_bias_control_dt-20-5'
'data_bias_control_t-20-6'
'data_bias_control_dt-20-6'
'data_bias_control_t-20-7'
'data_bias_control_dt-20-7'
'data_bias_control_t-20-8'
'data_bias_control_dt-20-8'
'data_bias_control_t-20-9'
'data_bias_control_dt-20-9'
'data_bias_control_t-20-10'
'data_bias_control_dt-20-10'};

%OUTFILE = '../data/cup_data_set.csv';
%fid = fopen(OUTFILE,'w') 
%fprintf(fid,'N,In,Out,dt,error,time,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19,C20,\n');
%fclose(fid);

T =[];
start = 1;

for k=1:length(LIST)
  out = create_csv_cup(sprintf('../data/data_bias_control/%s.mat',LIST{k}));
  L = max(size(out));
  W = min(size(out));
  T = [T; out zeros(L,25 - (W-1))];
  % lower precision for smaller files
  % T = [T; round(1e2*out)/1e2 zeros(L,25 - (W-1))];
  % writematrix(out,OUTFILE,'WriteMode','append')
end

C = array2table(T);
header = ["N","In","Out","dt","error","time","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20"];
C.Properties.VariableNames = header;
writetable(C,'../data/cup_data_set.csv');

function out = create_csv_cup(mat_file_name)
load(mat_file_name)
if ~exist('Results','var')
  error('no results variable containing controllers')
end 
for k=1:length(Results) 
  tmp = Results{k}; 
  if isempty(Info.args.readout)
    dt = 0;
  else
    dt = Info.args.readout(1);
  end
  out(k,:) = [length(tmp.bias), Info.args.in, Info.args.out, dt, tmp.err, tmp.time, tmp.bias']; % tmp.init_time, tmp.init_bias']; 
end
