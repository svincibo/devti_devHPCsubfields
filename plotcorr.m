function [r, p] = plotcorr(x,y,xname,yname,titlestrm, c)

% function plotcorr(x,y,xname,yname,titlestr)
%
% plots scatter plot of x against y and calculates r and p values
% xname,yname are optional arguments for axis labels
    
[r,p] = corrcoef(cat(2, x, y), 'Rows','pairwise');

idxy = (isnan(x)+isnan(y)) > 0;

linearCoef = polyfit(x(~idxy),y(~idxy),1);
linearFit = polyval(linearCoef,x(~idxy));
% plot(x(~idxy),y(~idxy), 'bo'); 
hold on
if strcmp(c, 'b')
    plot(x(~idxy),linearFit,'b--')
elseif strcmp(c, 'c')
    plot(x(~idxy),linearFit,'c--')
elseif strcmp(c, 'r')
    plot(x(~idxy),linearFit,'r--')
elseif strcmp(c, 'g')
    plot(x(~idxy),linearFit,'g:')
else
    plot(x(~idxy),linearFit,'k:')
end

% if size(r, 1) ~= 1
% r = r(2, 1); p = p(2, 1);
% end
statstr = sprintf('r = %.3f, p = %.3f',r,p);

r = r(2, 1);
p = p(2, 1);

if exist('xname','var')
    xlabel({xname; statstr})
else
    xlabel(statstr)
end


if exist('titlestr','var'), title(titlestr), end
if exist('yname','var'), ylabel(yname), end

