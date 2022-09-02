function JTanalysis(sensitivity)
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
n=floor(N/5);
d=[sort(y([1:n])) sort(y([n+1:2*n])) sort(y([2*n+1:3*n])) sort(y([3*n+1:4*n])) sort(y([4*n+1:N]))];
g=[ones(1,n) 2.*ones(1,n) 3.*ones(1,n) 4.*ones(1,n) 5.*ones(1,N-4*n)];
X=[d' g'];
jttrend(X)
end

