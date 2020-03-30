function fig = plot_strain(e_xx, e_yy, e_xy, X, Y, filename, epsilon)
    fig = figure(1);
    subplot(2,2,1);
    surf(reshape(e_xx,size(X)));
    colorbar();
    axis equal;
    view(2)
    subplot(2,2,2);
    surf(reshape(e_xy,size(X)));
    colorbar();
    axis equal;
    view(2)
    subplot(2,2,3);
    surf(reshape(e_yy,size(X)));
    colorbar();
    axis equal;
    view(2)
    subplot(2,2,4);
    ceval = concentration(filename, epsilon);
    surf(ceval(X,Y));
    cbar = colorbar();
    caxis([0,1]);
    axis equal;
    view(2);
end