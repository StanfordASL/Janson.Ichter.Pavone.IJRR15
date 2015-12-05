function plotChain(thetas, color, varargin)
    markerSize = 5;
    lineWidth = 2;
    lineStyle = '-';
    
    if length(varargin) >= 1 && ~isempty(varargin{1})
        lineWidth = varargin{1};
    end
    if length(varargin) >= 2 && ~isempty(varargin{2})
        lineStyle = varargin{2};
    end
    
    dim = size(thetas,2);
    linkLength = 1/dim/sqrt(2);
    theta = 0;
    w = ones(1,2)*0.5;
    for link = 1:dim
        theta = theta + thetas(link);
        v = w;
        w = v + [cos(theta) sin(theta)]*linkLength;
        plot([v(1) w(1)],[v(2) w(2)],'-o',...
            'Color',color,'MarkerSize',markerSize,'LineWidth',lineWidth,...
            'LineStyle',lineStyle);
    end
end
