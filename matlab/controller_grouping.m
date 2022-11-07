% load data
load  data_bias_control_dt-12-4.mat
% calc distance between controllers (ignoring readout time)
for k=1:2000, for l=k+1:2000, n(k,l)= norm(Results{k}.bias-Results{l}.bias); end; end;
% find controllers that are close 
for k=1:1999, close{k}=find(log10(n(k,k+1:2000))<1); end
% check if biases appear close, e.g. close{5} is not empty
bar(Results{5}.bias), hold on, pause
for k=1:length(close{5})
  bar(Results{5+close{5}(k)}.bias), pause
end
% visually the biases appear indistinguishable 