function [xhat,resnorm,residual] = bindingCurveFit(c,f,fitParams)
% bindingCurveFit(c,f,fitParams) fits a canonical binding curve and plots input data and fit curve
%
%   Inputs:
%   c = concentration series (nx1 or 1xn matrix)
%   f = signal from binding assay (nx1 or 1xn matrix)
%   fitParams = array of initialization parameters of form
%               [kd, amplitudeOfSignal, minSignal]
%   Outputs: [fitParams,resnorm,residual]
%
%   fits binding curve of form:
%           f = (c/(c+kd))*(amplitudeOfSignal-minSigna)+minSignal

xdata = c;
ydata = f;
n = length(xdata); % number of observations

% model
bindingCurve = @(x,xdata) (xdata./(xdata+x(1)))*(x(2)-x(3))+x(3);

% fit the least squares fit for this data
[xhat,resnorm,residual] = lsqcurvefit(bindingCurve, fitParams, xdata, ydata);

% generate curve using fit parameters
newx = logspace(4,ceil(log10(max(c))));
yhat = bindingCurve(xhat,newx);

% plot the simulated ydata vs xdata
semilogx(xdata, ydata,strcat('b','o'),'MarkerSize',8,'LineWidth',1);
hold on;
a = semilogx(newx, yhat, 'b', 'LineWidth', 2);
xlabel('[Protein] (nM)','FontSize',14);
ylabel('Fraction Signal','FontSize',14);
axis([1E-4 ceil(log10(max(c))) -.1 ceil(log10(max(f)))]);
legend([a],strcat('kD=',num2str(xhat(1),3)));

hold off;

end
