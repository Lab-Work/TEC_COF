function LH_plot3D(tScale, xScale, N,k, fd)

    surf(tScale, xScale, N',ones(size(k')),'EdgeColor', 'none');
    view([-110 30]);
    light('Position',[-15 300 14],'Style','local');
    zlabel({'vehicle count'});

    % % Create xlabel
    xlabel({'time (s)'});
    % 
    % % Create ylabel
    ylabel({'x (m)'});
end