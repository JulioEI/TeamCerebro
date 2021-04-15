function h = plotFill1(datax,datay,color1,typeDisp, n,linestyle)
% Mean and transparent error plot
% h = plotFill1(datax,datay,typeDisp)
% INPUT
%   datax: data x axis, same for all datay (1 x m).
%   datay: data matrix (n x m). Different entries must be in columns.
%   code: code for datay entries (1 x n).
%   typeDisp: Dispersion mode for datay, 'ic95' (default) or 'stdev'.

if ~exist('color1')
    color1 = [.0 .75 .75];
end

if ~exist('typeDisp')
    typeDisp = 'sem';
end

if length(datax) ~= size(datay,1)
    datay = datay';
    if length(datax) ~= size(datay,1)
        disp('Dimenssion no match');
       
    end
end

if ~exist('n')
    n =  size(datay,2);
elseif isempty(n)
    n =  size(datay,2);
end

if strcmp(typeDisp, 'ic95')
    f1 = 1.96/ sqrt(n);
elseif  strcmp(typeDisp, 'sd')
    f1 = 1;
elseif  strcmp(typeDisp, 'sem')
    f1 = 1/ sqrt(n);
end
[x1, y1] = fillformat(datax, nanmean(datay,2),...
    nanstd(datay,[],2) * f1);

if ~exist('linestyle')
    linestyle = '-';
end
%figure % membrane potential
hold on
fill(x1, y1, color1,'EdgeColor','none','faceAlpha',.2,'HandleVisibility','off');
plot(datax, nanmean(datay,2),'lineWidth',2,'color',color1,'lineStyle', linestyle);

end