% Yanning Li, Sep 02, 2015
% This function plot the 2D time space diagram of a junction
% If the junctin is connection, the timespace diagram of two links are
% concatenated together.
% If the junction is merge of diverge, it will plot a 1x2 subplots, left is
% the concatenated 1->2, and right is the 1->3 (merge).


function [C,h] =LH_plot2D_junc(varargin)

if nargin==11
    % merge or diverge, plot two figures
    middleX = varargin{1};
    T_cum = varargin{2};
    num_steps = varargin{3};
    tScale = varargin{4};
    xScaleLeft = varargin{5};
    xScaleRight = varargin{6};
    NLeft = varargin{7};
    kLeft = varargin{8};
    NRight = varargin{9};
    kRight = varargin{10};
    junc = varargin{11};
    
    subplot(1,2,1)
    contourf(tScale, xScaleLeft, kLeft',64,'LineColor', 'none')
    caxis([0 1]);
    hold on
    
    [~, ~] = contour(tScale, xScaleLeft, NLeft',30,'k');
    
    %color map transform
    colorbar('YTick',[0 0.5 1],'YTickLabel',{'Zero','Critical','Max'});
    % Use math symbols
%     yt_labels = {'0','\rho_c','\rho_m'};
%     cbar_pro(yt_labels);
    %colorbar('YTick',[0 fd.density(0) fd.kappa],'YTickLabel',{'Zero density','Critical density','Max density'});
    
    % Create xlabel
    xlabel('time (s)','fontsize',22);
    % Create ylabel
    if strcmp(junc.type_junc,'merge') || strcmp(junc.type_junc,'onrampjunc')
        ylabel(sprintf('Link %d to Link %d   (m)',junc.inlabel(1),junc.outlabel),'fontsize',22);
        title(sprintf('%d %d Merge to %d',junc.inlabel(1), junc.inlabel(2), junc.outlabel(1)),...
            'fontsize',24);
    elseif strcmp(junc.type_junc,'diverge')
        ylabel(sprintf('Link %d to Link %d   (m)',junc.inlabel,junc.outlabel(1)),'fontsize',22);
        title(sprintf('%d Diverge to %d %d',junc.inlabel(1), junc.outlabel(1), junc.outlabel(2)),...
            'fontsize',24);
    end
    
    % Plot T marks
    for i = 1:num_steps
        plot([T_cum(i+1) T_cum(i+1)],[middleX(1)-middleX(1)/50 middleX(1)+middleX(1)/50],...
            'k','LineWidth',2);
    end
    set(gca,'fontsize',18);
    hold off;
    
    
    subplot(1,2,2)
    contourf(tScale, xScaleRight, kRight',64,'LineColor', 'none')
    caxis([0 1]);
    hold on
    
    [C, h] = contour(tScale, xScaleRight, NRight',30,'k');
    
    %color map transform
    colorbar('YTick',[0 0.5 1],'YTickLabel',{'Zero','Critical','Max'});
    % Use math symbols
%     yt_labels = {'0','\rho_c','\rho_m'};
%     cbar_pro(yt_labels);
    %colorbar('YTick',[0 fd.density(0) fd.kappa],'YTickLabel',{'Zero density','Critical density','Max density'});
    
    % Create xlabel
    xlabel('time (s)','fontsize',22);
    % Create ylabel
    if strcmp(junc.type_junc,'merge') || strcmp(junc.type_junc,'onrampjunc')
        ylabel(sprintf('Link %d to Link %d   (m)',junc.inlabel(2),junc.outlabel),'fontsize',22);
        title(sprintf('%d %d Merge to %d',junc.inlabel(1), junc.inlabel(2), junc.outlabel(1)),...
            'fontsize',24);
    elseif strcmp(junc.type_junc,'diverge')
        ylabel(sprintf('Link %d to Link %d   (m)',junc.inlabel,junc.outlabel(2)),'fontsize',22);
        title(sprintf('%d Diverge to %d %d',junc.inlabel, junc.outlabel(1), junc.outlabel(2)),...
            'fontsize',24);
    end
    
    % Plot T marks
    % if even discretization over time, then transform to vector
    
    for i = 1:num_steps
        plot([T_cum(i+1) T_cum(i+1)],...
            [middleX(2)-middleX(2)/50 middleX(2)+middleX(2)/50],...
            'k','LineWidth',2);
    end
    set(gca,'fontsize',18) 
    hold off;
    
    
elseif nargin==8
    % one connection
    
    middleX = varargin{1};
    T_cum = varargin{2};
    num_steps = varargin{3};
    tScale = varargin{4};
    xScaleComb = varargin{5};
    NComb = varargin{6};
    kComb = varargin{7};
    junc = varargin{8};
    
    contourf(tScale, xScaleComb, kComb',64,'LineColor', 'none')
    caxis([0 1]);
    hold on
    
    [C, h] = contour(tScale, xScaleComb, NComb',30,'k');
    
    %color map transform
    colorbar('YTick',[0 0.5 1],'YTickLabel',{'Zero','Critical','Max'});
    % Use math symbols
%     yt_labels = {'0','\rho_c','\rho_m'};
%     cbar_pro(yt_labels);
    %colorbar('YTick',[0 fd.density(0) fd.kappa],'YTickLabel',{'Zero density','Critical density','Max density'});
    
    xlabel({'time (s)'});
    if strcmp(junc.type_junc,'one2one')
        ylabel('x (m)','fontsize',22);
        title(sprintf('Link %d to Link %d',junc.inlabel,junc.outlabel),'fontsize',24);
    elseif strcmp(junc.type_junc,'merge')   
        ylabel( sprintf('Link %d downflow and Link %d downflow',junc.inlabel(1),junc.inlabel(2)),...
                'fontsize',22 );
        title('Merge check priority rule','fontsize',24);
    elseif strcmp(junc.type_junc,'diverge')   
        ylabel( sprintf('Link %d upflow and Link %d upflow',junc.inlabel(1),junc.inlabel(2)),...
                'fontsize',22 );
        title('Diverge check distribution rule','fontsize',24);
    end
    
    % Plot T marks
    for i = 1:num_steps
        plot([T_cum(i+1) T_cum(i+1)],...
            [middleX-middleX/50 middleX+middleX/50],...
            'k','LineWidth',2);
    end
    set(gca,'fontsize',18) 
    hold off;
    
end

end






