%% Ver2_1_1
%           - Organized version 2_1
%
%
function MovieCarVer2_1_1(fileName, plotParam)

N     =  plotParam.N;
state =  plotParam.stateMat;
L     =  plotParam.L;       % Wheelbase


%% Folder to save the figure

currentFolder = pwd;
address =  strcat(currentFolder,'\SavedFigs\');


%% Parameters

trace   = 0;            % Trace length of loci

% x-y margin pads in the figure
xMargin = 0.3;
yMargin = 0.3;

% Name and address of the video file
fileType        = '.avi';
fullAddress     = strcat(address,fileName,fileType);

% Video settings
vid             = VideoWriter(fullAddress);
vid.Quality     = 100;      % A value between 0 to 100
vid.FrameRate   = 10;

% Color map:
cmap = repmat([255, 68, 0]./255, N,1);


%% Locus Movie

sizeFig     = [10 8];
position    = [2 2, sizeFig];
figure('Units', 'inches', 'Position', position);
axis equal
box on

X       = state(:,1:2:(2*N-1));     % x coordinates        
Y       = state(:,2:2:(2*N));       % y coordinates
Delta   = state(:,2*N+1:3*N) + state(:,3*N+1:4*N);  % Steering angles
Theta   = state(:,2*N+1:3*N);       % Heading angles

itrTot  = size(X,1);                % Number of iterations (movie frames)

open(vid);                          % Start recording

for itr = 1 : itrTot
    itr;
    hold on
    
    % Line transparancy
    alpha = linspace(0.0,0.95,trace);
    alpha = alpha(max(1,trace-itr):trace);
    
    % Plot initial positions as circles
    for i = 1 : N        
        if (itr ~= 1) && (itr < trace)
            scatter(X(1,i),Y(1,i), 200, cmap(i,:), ...
                'MarkerEdgeColor',cmap(i,:).^(alpha(1)), 'LineWidth',2);
        end
    end
    
    % Plot trajectories
    for i = 1 : N        
        jMin = max(1,itr-trace);
        for j = jMin : itr-1
            plot(X(j:j+1,i),Y(j:j+1,i),...
                'Color',cmap(i,:).^(alpha(j-jMin+1)), 'LineWidth',3); 
        end
    end
    
    % Draw cars and number them
    for i = 1 : N
        DrawCar([X(itr,i),Y(itr,i), Delta(itr,i), Theta(itr,i)], L); % Plot cars
        
        strNums = strtrim(cellstr(num2str(i,'%d'))); % Number cars
        text(X(itr,i),Y(itr,i),strNums,         ...
            'color', [0,0,0],                   ...
            'VerticalAlignment','top',          ...
            'HorizontalAlignment','right',      ...
            'FontWeight','demi','FontSize',22   );
    end
    
    % Label axes
    xlabel('x','FontWeight','demi');
    ylabel('y','FontWeight','demi');
    
    % Adjust margins
    set(gca, 'XLimMode', 'auto');
    set(gca, 'YLimMode', 'auto');
    axis equal
    xLim = get(gca,'XLim');    
    yLim = get(gca,'YLim');    
    set(gca, 'XLim', xLim + [-xMargin, xMargin]);
    set(gca, 'YLim', yLim + [-yMargin, yMargin]);
    
    hold off
    drawnow        
    writeVideo(vid, getframe(gcf));     % Write video frame
    
    if itr ~= itrTot(end), cla, end     % Do not wipe the last frame
end


close(vid);                             % Stop recording
























































