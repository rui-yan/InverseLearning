function hFig = visualizeRS(rs, model)
    varsToPlot = char('u, m', 'v, m');
    for i = 1:size(varsToPlot,1)
        hFig(i) = figure(i);
        pdeplot(model, 'XYData', rs(:,i), 'Contour', 'on');
        axis equal;
        title(varsToPlot(i,:))
        % scale the axes to make it easier to view the contours
        axis square
    end
end