
function [xfill,yfill]=fillformat(xvector, yMean, ySd)
% Return shadow to draw error. SM 2016
% [xfill,yfill]=fillformat(xvector, yMean, ySd)
% INPUT
%   xvector: dependent variable values (vector)
%   yMean: independent variable values (vector)
%   ySd: error of the independent variable in all points of xvector (vector) 
%   
% OUTPUT
%   xfill: x input for fill to draw the poligon (vector)
%   yfill: y input for fill to draw the polygon (vector)
%%
%
if size(xvector,1)>size(xvector,2)
    xvector=xvector';
end
if size(yMean,1)>size(yMean,2)
    yMean=yMean';
end
if size(ySd,1)>size(ySd,2)
    ySd=ySd';
end

% discard ymean
yMean=yMean(~isnan(yMean));
ySd=ySd(~isnan(yMean));
xvector=xvector(~isnan(yMean));

ySd(isnan(ySd))=0; % NaN to 0;

    xfill=[xvector xvector(end) fliplr(xvector) xvector(1)];
    yfill=[(yMean+ySd) (yMean(end)+ySd(end)) fliplr(yMean-ySd) (yMean(1)+ySd(1))];
end