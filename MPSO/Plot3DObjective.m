function [] = Plot3DObjective(func_handle, fname, fopt, fopt_loc, cmin, cmax, stepsize)
    
    [X, Y] = meshgrid(cmin:stepsize:cmax);

    Z = zeros(size(X));

    for i = 1 : numel(X)
        x = [X(i) Y(i)];
        Z(i) = func_handle(x);
    end

    C = del2(Z);

    fig = figure;
    set(fig, 'Visible', 'off');
    mesh(X,Y,Z,C,'FaceLighting','gouraud','LineWidth',0.5);
    hold on;
    plot3(fopt_loc(1),fopt_loc(2),fopt,'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',10);
    xlabel('x');
    ylabel('y');
    view(35,14);
    saveas(fig, strcat(fname, '.png'));
    
end