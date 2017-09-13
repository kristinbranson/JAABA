function vizLineRpts(ax,data,fld1,fld2,line)
% vizLineRpts(ax,data,fld1,fld2,line)
% Visualize experimental repeats for given stat/line.
%
% data: data struct
% fld1: field for x-axis of scatterplot. Defaults to day_num.
% fld2: field for y-axis of scatterplot
% line: line to visualize

if isempty(ax)
    figure;
    ax = axes;
end

if isempty(fld1)
    fld1 = 'day_num';
end 

tfLine = strcmp({data.line_name}',line);
tmpGF = {data.bdpOrNot}';
%tmpGF = repmat({'not'},numel(tfLine),1);
tmpGF(tfLine) = {line};

axes(ax);
x = OlyDat.Analysis.getNumericScalarFieldSafe(data,fld1);
y = OlyDat.Analysis.getNumericScalarFieldSafe(data,fld2);
h = gscatter(x,y,tmpGF);
set(h(1),'MarkerEdgeColor',[0.75 0.75 0.75]);
set(h(2),'MarkerEdgeColor',[0.3 0.3 0.8]);
set(h(3),'MarkerEdgeColor',[1 0 0],'MarkerSize',14);
hLeg = legend;
set(hLeg,'interpreter','none');
xlabel(fld1,'interpreter','none','FontWeight','bold');
ylabel(fld2,'interpreter','none','FontWeight','bold');
titlestr = sprintf('%s: %d exps',line,nnz(tfLine));
title(titlestr,'interpreter','none','FontWeight','bold');
grid on;

% plot ctl mean, std in x- and y- dirs
hold on;
tfCtl = [data.tfBDP]';
xCtl = [data(tfCtl).(fld1)];
yCtl = [data(tfCtl).(fld2)];
mux = mean(xCtl);
sdx = std(xCtl);
lox = mux-2*sdx;
hix = mux+2*sdx;
muy = mean(yCtl);
sdy = std(yCtl);
loy = muy-2*sdy;
hiy = muy+2*sdy;
%plot(mux,muy,'ob','MarkerSize',14,'MarkerFaceColor',[0 0 1]);
% plot([lox hix],[loy loy],'b');
% plot([lox hix],[hiy hiy],'b');
% plot([lox lox],[loy hiy],'b');
% plot([hix hix],[loy hiy],'b');
hold off;


