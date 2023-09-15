% calculate the signal energy for each force pulse strength and calculate
% the amplification
% Author: Pieter Try
% Last Edit: 15.09.2023


addpath functions\
thr_coeff = 1;          %"Härte" des Amplituden-Thresholds, 0.5 halbiert den Threshold
origFs = 4500;
noiseSTD = 0.00;
boolLeveldep = false;

myWavelet = "coif5";
myThresholdRule = 'Hard';
myDenoisingMethod = 'Bayes';
myNoiseEstimate = 'LevelDependent';

noiseWeightFactor = 2.0287;

type = "";

energies = [];
RMS = [];

% close all
%% Read Data

experiment = "TwinUBeamExperiment\";

path =  "C:\Users\trypi\sciebo\DAT meine Dateien\PRO_Promotion\701 Untersuchungen\TwinUBeamExperiment\";
if (~exist(path,'dir'))
    path =  "D:\Sciebo\DAT meine Dateien\PRO_Promotion\701 Untersuchungen\TwinUBeamExperiment\";
end

directory = uigetdir(path);
if directory == 0
    return;
end

filenames = funcGetFilePath(directory, ["1000"]);



% Get Data of Sensor Noise
% [noiseFile, noisePath] = uigetfile({'*.txt'}, 'Select Experiment Without Structure', path);
% noiseFile = "Noise2.txt";
% noisePath = path+experiment+"Noise\";

% figure;
% hold on
cellMat = cell(length(filenames),3); % Data, wavelet filtered Data, FS, time

parfor fnum = 1:length(filenames)
    datapath = string(directory) + "\"+ filenames(fnum);
    if  ~isnumeric(datapath)
        %% Import Data
        fileID = fopen(datapath);
        inputData = textscan(fileID,'%f %f %f %f %f %f %f %f %f %f %f %f %f ','Delimiter',';', HeaderLines = 0);
        fclose(fileID);
        %% convert Data
        inputData = cell2mat(inputData);
        % delete the first 20 datapoints
        inputData = inputData(10:end, :);
        %         time = zeros(length(inputData),1)
        %         for Z = 2:length(inputData)
        %         time(Z) = time(Z-1)+ inputData(Z,1)
        %         end
        inputData(:,1) = inputData(:,1) -inputData(1,1); %set time to zero

        time = inputData(:,1);
        sensorData = [inputData(:,3:8)];
        FS = 1/mean(diff(time))*1000*1000;
        sensorData(:,1:3) = sensorData(:,1:3) ./ (2^16-1) * 2*2;
        sensorData(:,4:6) = sensorData(:,4:6) ./ (2^16-1) * 2*250;

        sensorData = detrend(sensorData, 0);
       
        %% Calculate RMS and Energy per Impulse
        data1 = sqrt(sum(sensorData(:,1:3).^2,2));
%         data1 = sqrt(sensorData(:,1).^2+sensorData(:,2).^2+sensorData(:,3).^2);
        %         data2 = wdenoise(data1);


        interval = 0.5; %interval in seconds
        totalsteps = floor(length(data1)/floor(FS*interval));
        ratio =1;

        [D1, D2] =  calculateImpulseMag2(time, data1, interval, totalsteps,ratio);
        %         [D3, D4] =  calculateImpulseMag2(time, data2, interval, totalsteps,ratio);

        RMS(:,fnum) = D1;
        energies(:,fnum) = D2;
        %         RMSWAV(:,fnum) = D3;
        %         energiesWAV(:,fnum) = D4;
    else
        disp("FILE NOT FOUND: " + path+ filenames(fnum))
        %         return
    end
end



%% Energy

baseData = energies;

files = filenames;
filetypes = unique(extractBefore(files, length(files{1})-10));
% filetypes = filetypes(contains(filetypes,["Wstruct", "WOUT"]));
titles = filetypes;
% plotNames = ["(A) IMU Directly Mounted", "(B) IMU With Beam Structure","Signal Energy of (A) > 3\sigma_{Noise}^2","Signal Energy of (B) > 3\sigma_{Noise}^2","Amplification in Percent"];
plotNames = ["IMU Without Beam", "IMU With Beam","Amplification [%]", "Amplification With Noise Subtracted [%]"];

markers = [".","_"];
colors = ["#0072BD", "#D95319"];
fig = figure("Name","RMS MEAN and STD per Impulse");
fig.Position = [100 100 1000 600];

maxPulse = 2160;
numPul = 100;
quicksave = zeros(numPul,2);
noiseRed = zeros(2,1);

hold on
title("Energy Estimate of Vibration Magnitude With Force Pulses Off-Center")
for i = 1:length(filetypes)
    indexOfFiles = contains(files, filetypes(i));
    indexes = find(indexOfFiles);

    [numImpulses, numDim, ~] = size(baseData);

    meanVal = zeros(numImpulses, 1);
    stdVal = zeros(numImpulses, 1);
    
    for u = 1:numImpulses
        meanVal(u)  = mean(baseData(u,indexes));
        stdVal(u)  = std(baseData(u,indexes));
    end
    
    excludes = [1:11]; %,12:1:61];
    meanVal(excludes , :)  = [];
    stdVal(excludes , :)  = [];
    if size(meanVal,1 ) > numPul
        meanVal(numPul+1:end , :)  = [];
        stdVal(numPul+1:end , :)  = [];
    elseif size(meanVal,1 ) < numPul
        disp("too few datapoints");
        return;

    end
    
    xRange = [round([1:1:numPul].*maxPulse/numPul)];
    xlimits = round(linspace(21.6, maxPulse,15));
    
    noiseRed(i) = min(meanVal)*0.8;

%     errorbar(xRange, meanVal(:), stdVal(:), "Color",colors(i));
    errorbar(xRange, meanVal(:)-noiseRed(i), stdVal(:), "Color",colors(i));
    quicksave(:,i) = meanVal(:);
end

% plot(xRange, meanVal-meanValWAV-0.06);%(:)-0.06);

% xline(950, 'Color',"#D95319",'LineStyle', "--", "LineWidth",2)

% xline(1296, 'Color',"#0072BD",'LineStyle', ":", "LineWidth",2)
% yline(0.124*3)

hold off
%  xline(confX,'r')
ylabel("Acceleration Energy [g^2]")
set(gca,'FontSize',16, 'FontName', 'Times')
    set(gca, 'YScale', 'log')
%     xticklabels(xtickslab)
xlim([0 maxPulse])
xticks(xlimits)
xlabel("Force Pulse [\muNs]")


yyaxis right
hold on
% plot(xRange, quicksave(:,2)./quicksave(:,1)*100,'diamond', Color='black');
plot(xRange, (quicksave(:,2)-noiseRed(2))./(quicksave(:,1)-noiseRed(1))*100,'diamond', Color='black');
hold off
% plot(xRange, meanValWAV(:)-meanVal);

ylabel("Amplification [%]")
legend(plotNames)
hold off
grid on
grid minor
ax = gca;
ax.YAxis(2).Color = 'black';
%     set(gca, 'YScale', 'log')

% exportgraphics(fig, 'experiment_acceleration_energy_withVSwout.pdf', 'ContentType', 'vector');
