%Sensor Fusion Algorithm for fusion 6D IMU data with certain Resonance
%Modes using geometric realtions and adaptive combiner
%18.04.2023 Pieter Try
%input Parameter: data(n x 6 Mat), GRMat (W x 6) ,Fs, levels, motherWavelet, thresholdType, TPlevel, HPlevel,tresholdfilterBool, LevelDependent, windowFactor, showNoise)
% GRMat = geometric resonance matrix, which contains geometrix relations of
% each resonance mode as well as the frequency and name

function rslt = waveletResonanceSensorFusion(data, Fs, GRMat, levels, motherWavelet, thresholdType, noiseMat, thresholdFactor,tresholdfilterBool, debug)
%% set Prerequisits

%conversion function
conversion = @(factors , input) (factors(1) * input );%+ factors(2) .* abs(input) .* input);

if isempty(thresholdFactor)
thresholdFactor = 1;
end
% set Motherwavelet
if isempty(motherWavelet)
    motherWavelet = 'coif5';
end
noiseWeightExp = 2;

if isempty(thresholdType)
    thresholdType = 'h';
end

% Data length N
N = size(data,1);

% Check if input data is valid
numAxis = size(data,2);
if numAxis ~= 6
    disp("waveletResonanceSensorFusion: Input Data Dimensions Invalid")
    return;
end

% check GRMAT
numResModes = size(GRMat, 1);

% set Levels as log2(N);
if (levels == 0)
    levels = ceil(log2(N));
end

% calculate Frequency Bins of Levels
lowerBound = zeros(levels+1,1);
upperBound = zeros(levels+1,1);
lastF = Fs/2;
for i = 1:levels+1
    lowerBound(i) = lastF/2;
    upperBound(i) = lastF;
    lastF =lowerBound(i);
end
freqBins = [lowerBound, upperBound];
approxBin = freqBins(end,:);
detailBins = freqBins(1:end-1,:);

% set window to only a sec
windowFactor = Fs/N;


%% Decompose Wavelets
% calculate wavelet coefficients discrete wavelet transformation of IMU
% data
IMUwavDecs = cell(numAxis,1);
for i = 1:numAxis
    IMUwavDecs{i} = mdwtdec('c',data(:,i).*tukeywin(N,windowFactor),levels,motherWavelet);
end

% create empty wavelet object
emptyWavDec  = IMUwavDecs{1};
emptyWavDec.ca = emptyWavDec.ca * 0;
for i = 1:levels
    emptyWavDec.cd{i} = emptyWavDec.cd{i} .* 0;
end

% initialize result vektor
ResModeWavDecs = cell(numResModes,1);
for i = 1:numResModes
    ResModeWavDecs{i} = emptyWavDec;
end


%% Calculate weighted mean of Aproximmation coefficients
% % estimate Noise for each data stream
% if (LevelDependent)
%     noise = ones(numAxis, 1);
%     for u = 1:numAxis
%         noise(u) =  median(abs(IMUwavDecs{u}.ca))/0.6745;
%     end
% end
%
% % calculate weight based on inverse square of noise
% noiseWeights =calcPowerNoiseWeights(noise,noiseWeightExp);
%
% % weighted average of all data streams
% weightedAverage=IMUwavDecs{1}.ca *0;
% for u = 1:numAxis
%     weightedAverage =  weightedAverage + noiseWeights(u) * IMUwavDecs{u}.ca;
% end
% xdecdenoised.ca = weightedAverage;

%% Main Loop Over Resonance Modes
for rNum = 1 : numResModes
    resonanceName = GRMat{rNum, 1};
    transformationMatrix = GRMat{rNum, 2};
    resonanceFrequency = GRMat{rNum, 3};

    % check if resonance mode lies in approximation
    if (resonanceFrequency < approxBin(2))
        disp("ERROR: Resonance Mode does not lie in approximation levels")
        return
    else
        % find wavelet level of resonance frequency
        resLvl = find((detailBins(:,1) <= resonanceFrequency) & (resonanceFrequency < detailBins(:,2)), 1);
        GRMat{rNum,4} = resLvl;

        if isempty(resLvl)
            disp("waveletResonancenSensorFusion: given Resonance Frequency not found, abort")
            return
        end

        % get the detail coefficients of all axis in a matrix
        coeffSize = size(IMUwavDecs{1}.cd{resLvl},1);
        coeffMatrix = zeros(coeffSize,6);
        for i = 1:numAxis
            coeffMatrix(:,i) = IMUwavDecs{i}.cd{resLvl};
        end
        
        Mat4Testing = coeffMatrix;

       % transform using resonance Matrix
       % and estimate Noise
        noiseSD = zeros(numAxis, 1);
        for i = 1:numAxis
            coeffMatrix(:,i) =  conversion(transformationMatrix(:,i) ,coeffMatrix(:,i));
            
            if isempty(noiseMat)
                if (i < 4) && (resLvl ~= 2)  %make noise estimation level independent for low frequency acceleration data
                    noiseSD(i) =  median(abs(conversion(transformationMatrix(:,i) ,IMUwavDecs{i}.cd{1})))/0.6745;
                else
                    noiseSD(i) =  median(abs(coeffMatrix(:,i)))/0.6745;
                end
            else
                % user predetermined noise levels
                noiseSD(i) = abs(conversion(transformationMatrix(:,i) , noiseMat(resLvl, i)));
            end
        end

%         % estimate noise for all axis with data to fuse
%         noiseSD = ones(numAxis, 1);
%         for i = 1:numAxis
% %             if (i < 4) && (resLvl > 3)%make noise estimation level independent for acceleration data
% %                 noiseSD(i) =  conversion(resonanceMatrix(:,i) ,median(abs(detrend(IMUwavDecs{i}.cd{1},0)))/0.6745);
% %             else % and level-dependent for dgyro/dt data
%                 noiseSD(i) =  median(abs(detrend(coeffMatrix(:,i),0)))/0.6745;
% %             end
%         end

        %% automatically calculate transformation value
        if (debug)
            disp("###" + resonanceName + " Frequenz: " + resonanceFrequency +" Hz")

            idxRef = find(transformationMatrix(1,:) == 1,1);
            idxOther = find((transformationMatrix(1,:) ~= 0) & (transformationMatrix(1,:) ~= 1));
            axisNames = ["ax", "ay", "az", "wx", "wy", "wz"];
            if ~isempty(idxOther)
                %calculate factor for minimum error between peaks
                for v = 1: length(idxOther)
                    % toBeConverted =  coeffMatrix(:,idxOther(v));
                    % refData = coeffMatrix(:,idxRef);
                    thresh = noiseSD(idxOther(v)) * sqrt(2*log(size(Mat4Testing(:,idxOther(v)),1)));
                    toBeConverted =  wthresh(Mat4Testing(:,idxOther(v)), 'h', thresh);

                    thresh = noiseSD(idxRef) * sqrt(2*log(size(Mat4Testing(:,idxRef),1)));
                    refData = wthresh(Mat4Testing(:,idxRef), 'h', thresh);

                    pointsWithData = (toBeConverted ~= 0 ) & (refData ~= 0);

                    rms =@(fact) sqrt(sum((conversion(fact,toBeConverted(pointsWithData)) - refData(pointsWithData)).^ 2)); %min MSE function
                    factors = fminsearch(rms,[1;1]); %find the x

                    disp(":::transform from " + axisNames(idxOther(v))+ " to " + axisNames(idxRef) + "     Number of Points: " +sum(pointsWithData))
                    disp("Faktoren:  " + factors(1) + "    " + factors(2));
                    disp("Error:  " + rms(factors))
                    disp("__________________________")

                    if (false)
                    figure("Name", resonanceName + " Frequenz: " + resonanceFrequency +" Hz  " + "transform from " + axisNames(idxOther(v))+ " to " + axisNames(idxRef));
                    plot([refData,toBeConverted])
%                     plot([refData,conversion(factors, toBeConverted)])
                    end
                end
            end
        end

        %remove all zero elements in resonanceMatrix
        indexOfEmtpy = (transformationMatrix(1,:) == 0);
        coeffMatrix(:,indexOfEmtpy) = [];
        noiseSD(indexOfEmtpy) = [];
        numFusionAxis = size(coeffMatrix,2);

        %% weighted average of all data streams
        % calculate weight based on inverse square of noise
        noiseWeights  = calcPowerNoiseWeights(noiseSD, 2);
        disp("NOISE OF COMPONENTS: " + noiseSD)
        %calculate the estimated noise of the fusion
        GRMat{rNum,5} = sqrt(sum(noiseWeights.^2 .* noiseSD.^2));

        weightedAverage = coeffMatrix(:,1) *0;
        for i = 1:numFusionAxis
            weightedAverage =  weightedAverage + coeffMatrix(:,i) .* noiseWeights(i);
        end

        ResModeWavDecs{rNum}.cd{resLvl} = weightedAverage;

        
    end
end

% disp(GRMat)

%% Universal Threshold
if (tresholdfilterBool)
    for rNum = 1:numResModes
        %         for lev = 1:levels
        lev = GRMat{rNum,4};
        if abs(sum(ResModeWavDecs{rNum}.cd{lev},'all')) > 0

%             newNoiseEst = median(abs(ResModeWavDecs{rNum}.cd{lev}() ))/0.6745;
            newNoiseStdDev = GRMat{rNum,5};

            thresh = newNoiseStdDev * sqrt(2*log(size(ResModeWavDecs{rNum}.cd{lev},1)))*thresholdFactor;
            filtered = wthresh(ResModeWavDecs{rNum}.cd{lev}, thresholdType, thresh);

            disp("Noise Estimate in Threshold Filter = " + median(abs(ResModeWavDecs{rNum}.cd{lev}-filtered))/0.6745);
            ResModeWavDecs{rNum}.cd{lev} = filtered;


            %                 if (newNoiseEst >= GRMat{rNum,5})
            %                     newNoiseEst = GRMat{rNum,5};
            %                 end
            % thresh = sqrt(2*log(size(ResModeWavDecs{rNum}.cd{lev},1)));
            % ResModeWavDecs{rNum}.cd{lev} = myThresholdFilter(ResModeWavDecs{rNum}.cd{lev}, newNoiseEst, thresholdType, thresh);

        end
        %         end
        if true
        disp("Mode " + rNum+ " Noise: "+ newNoiseStdDev)
        end
    end
end


% %% Low Pass Filter through coefficient elimination
% for tp = 1:TPlevel
%     xdecdenoised.cd{tp} = xdecdenoised.cd{tp} .*0;
% end
%
% %% High Pass Filter through coefficient elimination
% % if (HPlevel == 1)
% %     xdecdenoised.ca = xdecdenoised.ca .*0;
% if (HPlevel >= 1 )
%     %     xdecdenoised.ca = xdecdenoised.ca .*0;
%     for hp = 1:HPlevel
%         xdecdenoised.cd{end+1-hp} = xdecdenoised.cd{end+1-hp} .*0;
%     end
% end

%% DEBUGGING
showWavCoeffPre = false;
if showWavCoeffPre
    plotWavelets(xdecdenoised);
end

if false %debug
    disp("Numb Sensor Signals: " + length(IMUwavDecs));
    disp("Old Noise Variance:  "+ median(abs(IMUwavDecs{1}.cd{1}))/0.6745);
    disp("New Noise Variance:  "+ median(abs(xdecdenoised.cd{1}))/0.6745);
    disp("Quotient:   " +  (median(abs(IMUwavDecs{1}.cd{1}))/0.6745)/(median(abs(xdecdenoised.cd{1}))/0.6745))
end

%% Rücktransformation
rslt = zeros(N,numResModes);
for i = 1:numResModes
    rslt(:,i) = mdwtrec(ResModeWavDecs{i});

end
end










%% FUNCTIONS
function rslt = myThresholdFilter(input, sdev, type, threshold_value)
input = input./sdev;
threshold = ones(length(input),1) .* threshold_value ;
rslt = wthresh(input, type, threshold);
rslt = rslt .* sdev;
end


function rslt = showWavCoefficients(wav, FS)
lengthInS = wav.dataSize(1)/FS;

colormap jet;
cmap=colormap;

figure("Name","Wavelet Coefficients");
t = linspace(0,lengthInS, length(wav.ca()));
plot(t, wav.ca(),"Color",'Black');
leg = ["CA"];
hold on
for i = 1:wav.level
    Plot_color=cmap(ceil(i/wav.level*length(cmap)),:);
    t = linspace(0,lengthInS, length(wav.cd{i}));
    P(i) = plot(t, wav.cd{i},'Color', Plot_color);
    leg = [leg, "CD" + string(i)];
end
hold off
legend(P,leg)
grid
end

%% ARCHIVE






















%% Energie Filter
%
% if (sum(frequencyFilter) ~= 0)
%     if (length(frequencyFilter) ~= levels+1)
%        M = length(frequencyFilter);
%        told = 0:1/(M-1):1;
%        tnew = 0:1/(levels):1;
%        frequencyFilter = interp1(told, frequencyFilter, tnew);
%     end
%     xdecdenoised.ca = xdecdenoised.ca .* frequencyFilter(1);
%     for i = 1:levels
%         xdecdenoised.cd{i} = xdecdenoised.cd{i} .* frequencyFilter(end+1-i);
%     end
% end

% estimate new noise
%     newNoiseEst = median(abs(weightedAverage))/0.6745;
%
%     %% Universal Threshold
%     if (tresholdfilterBool)
%         thresh = sqrt(2*log(size(weightedAverage,1)))*0.95;
%         xdecdenoised.cd{lev} = myThresholdFilter(weightedAverage, newNoiseEst, thresholdType, thresh);
%     else
%         xdecdenoised.cd{lev} = weightedAverage;
%     end

%%
function rslt = plotWavelets(xdec)

lvls = xdec.level;

figure;


t = linspace(0,600, length(xdec.ca));
plot(t, xdec.ca);

grid on
title('Plots of Wavelet Coefficients')
xlabel('Time in Seconds')
ylabel('Wavelet Amplitude [mg]')

f = zeros(lvls+1,1);
lastF = 3600/2;
for freq = 1:lvls+1
    %     f(freq) = (Fs/2 - f(freq-1))/(2^(lvls-freq-1));
    f(freq) = lastF/2;
    lastF = f(freq);
end
f = flip(f);
plotnames = num2str(f);

set(gca,'FontSize',16, 'FontName', 'Palatino Linotype')
hold on
for graphs = 1:lvls

    waveCoeff = cell2mat(xdec.cd(end-graphs+1));
    t = linspace(0,600, size(waveCoeff,1));
    plot(t, waveCoeff);

    grid on
    title('Plots of Wavelet Coefficients')
    xlabel('Time in Seconds')
    ylabel('Wavelet Amplitude [mg]')

    set(gca,'FontSize',16, 'FontName', 'Palatino Linotype')
end
hold off
legend(plotnames)
end

%%
function rslt = plotCoeffEnergy(xdec1, xdec2, Fs)
lvls = xdec1.level;

energies = zeros(lvls+1,1);
energies2 = zeros(lvls+1,1);

energies(1) = sqrt(sum(abs(xdec1.ca).^2,'all'));
energies2(1) = sqrt(sum(abs(xdec2.ca).^2,'all'));

for graphs = 1:lvls
    waveCoeff = cell2mat(xdec1.cd(lvls-graphs+1));
    energies(1+graphs) = sqrt(sum(abs(waveCoeff).^2,'all'));

    waveCoeff2 = cell2mat(xdec2.cd(lvls-graphs+1));
    energies2(1+graphs) = sqrt(sum(abs(waveCoeff2).^2,'all'));
end

%
% for graphs = 1:lvls
%     waveCoeff = cell2mat(xdec1.cd(graphs));
%     energies(1+graphs) = sqrt(sum(abs(waveCoeff).^2,'all'));
%
%     waveCoeff2 = cell2mat(xdec2.cd(graphs));
%     energies2(1+graphs) = sqrt(sum(abs(waveCoeff2).^2,'all'));
% end


figure;

% f(1) = (Fs/2)/2^lvls;
f = zeros(lvls+1,1);
lastF = Fs*2;
for freq = 1:lvls+1
    %     f(freq) = (Fs/2 - f(freq-1))/(2^(lvls-freq-1));
    f(freq) = lastF/2;
    lastF = f(freq);
end
f = flip(f);

semilogx(f, energies');

hold on
semilogx(f,energies2');
hold off

grid on
title('Plot of Energy of Wavelet Coefficients')
xlabel('Frequency in [Hz]')
ylabel('Energy')
set(gca,'FontSize',16, 'FontName', 'Palatino Linotype')
end

function y = phaseShiftTimeSignal(x, shiftAmount)
    % Inputs:
    %   x: Input signal (vector or matrix)
    %   fractionalPhase: Fractional phase shift (in radians)

    % Perform circular shift to introduce fractional phase shift
    y = circshift(x, shiftAmount);
end

function rslt = calcPowerNoiseWeights(noise, exp)
noisemod = @(x) x.^(-exp);
% noisemod = @(x) -1* (factor*x);
% noisemod = @(x) -1* log(factor*x);
if any(noise == 0)
    disp("waveletResonanceSensorFusion: error calculating noise based weights")
    return
end

rslt = noisemod(noise) ./ sum(noisemod(noise),"all");
end
