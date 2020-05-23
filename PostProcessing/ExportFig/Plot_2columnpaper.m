function Plot_2columnpaper(h1, filename, Outpu_Extension, Resolution, x, y, String_X, String_Y, Title, Legend, Grid_Maj_Min, Figure_Type, Line_Style, Color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set up figure size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% By default, matlab prints a figure as 4 (width) x 3 (height) inches.
%% On US paper, a two-column figure is typically 3.5 inches wide and a
%% one-column figure is typically 6.5 inches wide (heights will vary).
%% Importing a figure of the default size into LaTeX and then scaling it to
%% fit makes the fonts shrink or expand, and makes the figure look bad.
%% This stuff sets up the figure area to a specific width x height.

if isempty(h1)
    h = figure;                                        % select figure pane
else
    h = h1;
end
set(gcf, 'PaperUnits', 'inches');              % units
width=3.5;                                     % width of figure
height=width*0.75;                             % 75% width keeps default ratio
set(gcf, 'PaperPositionMode', 'manual');       % turn off 'auto' so it doesn't resize upon printing
papersize = get(gcf, 'PaperSize');             % for US paper, 8.5 x 11
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
set(gcf, 'PaperPosition', [left bottom width height]);  % middle of page

%% By default, matlab surrounds a figure with extra whitespace.
%% This stuff sets up the 3 axes to contract that whitespace.
%% See http://nibot-lab.livejournal.com/73290.html

%set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
set(gca,'LooseInset',get(gca,'TightInset'));    % this works better
% Create axes
% axes1 = axes('Parent',h,'YGrid','on','XGrid','on',...
%     'Units','centimeters',...
%     'Position',[1.2 1.2 13 11],...
%     'LineWidth',0.75,...
%     'FontWeight','bold',...
%     'FontSize',12,...
%     'DataAspectRatio',[500 500 1],...
%     'ALim',[0 1]);
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
markers = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};
markers = markers(randperm(length(markers)));
L = size(y,1);
temp11 = rand(L,3);
switch Figure_Type
    case 'plot'
        for i=1:L
            if isempty(Line_Style)
                Line_Style1 = [ '-' cell2mat(markers(i))];
            else
                Line_Style1 = Line_Style;
            end
            if isempty(Color)
                Color1 = temp11(i,:);
            else
                Color1 = Color;
            end
            plot(x, y(i,:), Line_Style1, 'Color', Color1, 'LineWidth', 0.7, 'MarkerSize', 5)    % plot the computed Gaussian fit
            hold on
        end
    case 'semilogy'
        for i=1:L
            if isempty(Line_Style)
                Line_Style1 = [ '-' cell2mat(markers(i))];
            else
                Line_Style1 = Line_Style;
            end
            if isempty(Color)
                Color1 = temp11(i,:);
            else
                Color1 = Color;
            end
            semilogy(x,y(i,:), Line_Style1, 'Color', Color1, 'LineWidth',0.7, 'MarkerSize', 3)
            hold on
        end
    case 'hist'
        hist(y)
    case 'semilogx'
        for i=1:L
            if isempty(Line_Style)
                Line_Style1 = [ '-' cell2mat(markers(i))];
            else
                Line_Style1 = Line_Style;
            end
            if isempty(Color)
                Color1 = temp11(i,:);
            else
                Color1 = Color;
            end
            semilogx(x,y(i,:), [Line_Style1 Color1], 'LineWidth',0.7, 'MarkerSize', 3)
            hold on
        end
    case 'bar'
%         for i=1:L
%             if isempty(Line_Style)
%                 Line_Style1 = [ '-' cell2mat(markers(i))];
%             else
%                 Line_Style1 = Line_Style;
%             end
%             if isempty(Color)
%                 Color1 = temp11(i,:);
%             else
%                 Color1 = Color;
%             end
%             h = bar(x, y(i,:), 0.4);
%             set(h, 'FaceColor', Color1)
%             
% %             set(h, 'MarkerFaceColor','red','Marker','square')
% %             set(h, 'MarkerSize', 3)
%             hold on
%         end
        bar(x, y, 0.5)
        ax = gca;
        set(ax, 'XTick', x)
    case 'stem'
        for i=1:L
            if isempty(Line_Style)
                Line_Style1 = [ '-' cell2mat(markers(i))];
            else
                Line_Style1 = Line_Style;
            end
            if isempty(Color)
                Color1 = temp11(i,:);
            else
                Color1 = Color;
            end
            stem(x, y(i,:), [Line_Style1 Color1], 'LineWidth',0.7, 'MarkerSize', 3);
%             set(h, 'MarkerFaceColor','red','Marker','square')
%             set(h, 'MarkerSize', 3)
            hold on
        end
    otherwise
        disp('This option is used to save the figure')
end
hold off
%% Axes Labels
if ~isempty(String_X)
xlabel(String_X, 'FontWeight','bold','FontSize',6);     % x-axis label
end
if ~isempty(String_Y)
    ylabel(String_Y, 'FontWeight','bold','FontSize',6);     % y-axis label
end
if ~isempty(Title)
    title(Title,  'FontWeight','bold','FontSize',6); % title
end
set(gca, 'FontWeight','bold','FontSize',6);             % change the font for the ticks

if ~isempty(Legend)
    hleg = legend(Legend);
%     set(hleg,'Box','off')
    set(hleg, 'Location', 'SouthEast')
end
% NorthInside plot box near top
% SouthInside bottom
% EastInside right
% WestInside left
% NorthEastInside top right (default for 2-D plots)
% NorthWestInside top left
% SouthEastInside bottom right
% SouthWestInside bottom left
% NorthOutsideOutside plot box near top
% SouthOutsideOutside bottom
% EastOutsideOutside right
% WestOutsideOutside left
% NorthEastOutsideOutside top right (default for 3-D plots)
% NorthWestOutsideOutside top left
% SouthEastOutsideOutside bottom right
% SouthWestOutsideOutside bottom left
% BestLeast conflict with data in plot
% BestOutsideLeast unused space outside plot
% axis([0 12 0.0180 0.0215])
%% Grid
if ~isempty(Grid_Maj_Min)
    grid on
    if strcmp(Grid_Maj_Min, 'Min') || strcmp(Grid_Maj_Min, 'Minor')
        grid minor
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% done creating figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% print figure to file (-deps for greyscale, -depsc for color, -dpdf for pdf, etc.)
% print('-f1','-dpdf',filename);
if ~isempty(filename) || ~isempty(Outpu_Extension)
    for i=1:length(Outpu_Extension)
        Outpu_Extension_i = cell2mat(Outpu_Extension(i));
        if strcmp(Outpu_Extension_i, 'fig')
            saveas(h, [filename '.fig'])
        else
            print(h,['-r' num2str(Resolution)], ['-d' Outpu_Extension_i], filename);
        end
    end
end