function rslt = getNoisePerLvl(data, Fs, levels, motherWavelet, debug)


%% set Prerequisits

% set Motherwavelet
if isempty(motherWavelet)
    motherWavelet = 'coif5';
end


% Data length N
N = size(data,1);

% Check if input data is valid
numAxis = size(data,2);
if numAxis ~= 6
    disp("waveletResonanceSensorFusion: Input Data Dimensions Invalid")
    return;
end

% set Levels as log2(N);
if (levels == 0)
    levels = ceil(log2(N));
end

noiseMatrix = zeros(levels+1, numAxis);

%% Decompose Wavelets
% calculate wavelet coefficients discrete wavelet transformation of IMU
% data
IMUwavDecs = cell(numAxis,1);
for i = 1:numAxis
    IMUwavDecs{i} = mdwtdec('c',data(:,i),levels,motherWavelet);
end

%% Main Loop Over Levels
for i  = 1 : levels

    for u = 1:numAxis
        noiseMatrix(i,u) = median(abs(IMUwavDecs{u}.cd{i}))/0.6745;
%         noiseMatrix(i,u) = std(IMUwavDecs{u}.cd{i});
    end
end
for u = 1:numAxis
        noiseMatrix(end,u) = median(abs(IMUwavDecs{u}.ca))/0.6745;
%     noiseMatrix(end,u) = std(IMUwavDecs{u}.ca);
end
rslt = noiseMatrix;
end



