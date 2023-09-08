% Function to draw the fft of multiple sources in 3 axis, 6axis?
% input paramters: data cell array, FS vector, title, legend, AxisNames, X
% range, yrange, Filter

function rslt = showMultiAxisFFT(data, FS, titlename, dataNames, AxisNames,xlabels, ylabels,  xRange, yRange, filter)
if isempty(data)
    disp("ERROR: showMultiAxisFFT: Data Input Parameter is Empty")
end

numAxis = size(data{1},2);
numColumns = ceil(1+(numAxis-3)/3);
numSources = size(data,1);

fig = figure("Name",titlename);

if numAxis <= 3
    t = tiledlayout(numAxis,1);
elseif numAxis > 3
    t = tiledlayout(3, numColumns);
end

t.Padding = 'compact';
t.TileSpacing = 'compact';

% for loop for different Axis
for AxIndex = 1:numAxis
    tiles(AxIndex ) = nexttile;
    hold on
    %for loop for multiple sources
    for sIndex = 1:numSources
        %choose data
        plotData = data{sIndex}(:,AxIndex);
        Fs = FS(sIndex);

        %         L = size(plotData,1);
        [f, P1] = calcOneSidedFFT(plotData, [],Fs, false);
        P1 = sgolayfilt(P1,filter(1), filter(2));

        %         f = Fs *(0:(L/2))/L;
        %         P1 = P1 ./ max(P1);
        plot(f,P1)

        %      set(gca, 'YScale', 'log')
    end
    grid on
    title(AxisNames(AxIndex))
    if (i==1)
        legend(dataNames)
    end
    set(gca,'FontSize',16, 'FontName', 'Times')
    grid on
    grid minor
    xlabel(xlabels(AxIndex));
    ylabel(ylabels(AxIndex));
    if ~isempty(xRange)
        xlim(xRange)
        xticks(xRange(1):100:xRange(2))
    end
    if ~isempty(yRange)
        ylim(yRange)
    end
end
linkaxes(tiles, 'x' )
% exportgraphics(fig, 'output.pdf', 'ContentType', 'vector');
end
