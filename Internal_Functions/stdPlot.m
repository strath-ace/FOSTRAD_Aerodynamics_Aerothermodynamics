function stdPlot(Labels, Title, Legend, Annotation)
% stdPlot(Labels, Title, Legend, Annotation)
% Labels = {xlabel, ylabel, zlabel (if existing)} 
% Title = {Title}

fig = gcf;
fig.Units =  'pixels';
fig.Position = [100, 50, 850, 700];
xlabel(Labels{1}, 'Interpreter', 'Latex')
ylabel(Labels{2}, 'Interpreter', 'Latex')
ax = gca;
ax.TickLabelInterpreter = 'Latex' ;
if length(Labels) > 2 
    zlabel(Labels{3}, 'Interpreter', 'Latex')
    view([-50 30])
    if exist('Annotation','var')
        dim = [.20 .60 .3 .3];    
    end
else
    dim = [.60 .60 .3 .3];    
end
if exist('Annotation','var')
    annotation('textbox',dim,'String',Annotation,'FitBoxToText','on','FontSize', 14, 'FontName','Sans', 'Interpreter', 'Latex');
end

if exist('Legend','var')
    if ~strcmpi(Legend,['off'])
%         legend(Legend, 'Interpreter', 'Latex','Location','SouthOutside')   
        legend(Legend, 'Interpreter', 'Latex','Location','NorthEast')
    end
end
title(Title,'FontSize',18, 'Interpreter', 'Latex')
set(gca, 'FontSize', 18)

grid on;
set(ax.XLabel, 'FontSize', 18)
set(ax.YLabel, 'FontSize', 24)


return


annotation('textarrow',x,y,'String','y = x ','FontSize', 14, 'FontName','Sans', 'Interpreter', 'Latex')
%


% colorbar;
% c= colorbar;
% set(gca, 'CLim', [-MaxC MaxC]);
% set(gca, 'CLim', [-0.01 0.01]);
% ylabel(c,'Local Radius [m] ','FontSize',16)