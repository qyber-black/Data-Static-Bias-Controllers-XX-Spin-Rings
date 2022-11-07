function [tau,taulog]=JTanalysis_vsalpha(sensitivity,s,alpha)
%
% This function does the Jonckheere-Terpstra analysis of the increasing
% error versus increasing (log)sensitivity of the performance relative to
% coupling errors in spin rings. 
% At the Matlab prompt, one of the data files from Frank should be loaded;
% for example
%
% >> load data_dt-8-3.mat
%
% This brings the sensitivity data in the workspace. The above is the data
% of a 8-spin ring, with transfer from |1> to |3>, where everything is
% "windowed" over an interval dt. 
%
% The difference between this 'version salpha' and the 'version s' is that
% the significance level alpha can be adjusted 
% (before it was set to the default value of 0.05). 
%
sensitivity
% Define basic error and sensitivity data
x=log(sensitivity.error);
y=log((sensitivity.dpdJ_norm)./(sensitivity.error));
z=log(sensitivity.dpdJ_norm);
% Plot data
[dontcare N]=size(sensitivity.error);
h=plot([1:N],z,[1:N],y,[1:N],x)
% Set color, line thickness, etc.
set(h(1),'Color',[218 231 89] ./ 255); % color yellow-green
set(h(2),'Color',[141 53 11]/255);     % color brown
set(h(3),'Color',[0 0 1]);             % color blue
set(h(3),'LineWidth',1);               % make error plot a little ticker
legend('log(sensitivity)','log(logarithmic sensitivity)','log(error)','location','SouthEast');
xlabel('controller index')
% Now we start the Jonckheere-Terpstra statitics
n=floor(N/s);
d=[sort(y([1:n]))];
for k=2:(s-1)
    d=[d sort(y([(k-1)*n+1:k*n]))];
end
d=[d sort(y([(s-1)*n+1:N]))];
g=[];
for k=1:(s-1)
    g=[g k.*ones(1,n)];
end
g=[g s.*ones(1,N-(s-1)*n)];
X=[d' g'];
jttrend(X,alpha,[1:n])
% Now the kendal tau for the log sensitivity
taulog=corr(x',y','type','kendall');
tau=corr(x',z','type','kendall');
end

