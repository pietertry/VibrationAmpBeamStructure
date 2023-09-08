%This script reads and analysis sensor data for multi-axial IMUs for the
%article: "A Vibration Sensing Device Using a 6-axis IMU and an
%Optimized Beam Structure for Activity Monitoring"
%Author: Pieter Try 08.09.2023
%pieter.try@w-hs.de
%% 
addpath functions\
thr_coeff = 0;          %"Härte" des Amplituden-Thresholds, 0.5 halbiert den Threshold
tplevel = 0;
hplevel = 8;
noiseSTD = 0.00;
boolLeveldep = true;

myWavelet = [];
myThresholdRule = 's';
myDenoisingMethod = 'Bayes';
myNoiseEstimate = 'LevelDependent';
thresholdFactor = 0.67;
noiseWeightFactor = 2.0287;

type = "";
% close all
%% Read Data

experiment = "TwinUBeamExperiment\";

path =  "C:\Users\trypi\sciebo\DAT meine Dateien\PRO_Promotion\701 Untersuchungen\";
if (~exist(path,'dir'))
    path =  "D:\Sciebo\DAT meine Dateien\PRO_Promotion\701 Untersuchungen\";
end

[woutFile, woutPath] = uigetfile({'*.txt'}, 'Select Experiment Without Structure', path+experiment);
% if woutFile == 0
%     return;
% end
[withFile, withPath] = uigetfile({'*.txt'}, 'Select Experiment With Structure', path+experiment);
% if withFile == 0
%     return;
% end

if ((woutFile == 0) & (withFile == 0))
    disp("ABORTED: No data chosen.")
    return
end

% Get Data of Sensor Noise
% [noiseFile, noisePath] = uigetfile({'*.txt'}, 'Select Experiment Without Structure', path);
noiseFile = "Noise_6667D.txt";
% noiseFile = "Noise2.txt";

noisePath = path+experiment+"Noise\";

filenames = [string(woutFile), string(withFile), string(noiseFile)];
paths = [string(woutPath), string(withPath), string(noisePath)];
filenames(filenames == "0") = [];
paths(paths == "0" ) = [];

cellMat = cell(length(filenames),5); % Data, wavelet filtered Data, FS, time

for fnum = 1:length(filenames)
    datapath = funcGetFilePath(paths(fnum), filenames(fnum));
    if  ~isnumeric(datapath)
        %% Import Data
        fileID = fopen(paths(fnum) + datapath);
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
        inputData(:,1) = inputData(:,1) - inputData(1,1); %set time to zero

        time = inputData(:,1);
        sensorData = [inputData(:,3:8)];
        FS = 1/mean(diff(time))*1000*1000;
        sensorData(:,1:3) = sensorData(:,1:3) ./ (2^16-1) * 2*2;
%%%%%%% SWITCH X AND Y BECAUSE OF SIMULATION
%         sensorData(:,1:2) = [sensorData(:,2), sensorData(:,1)];
        sensorData(:,4:6) = sensorData(:,4:6) ./ (2^16-1) * 2*250;

        %% Wavelet Filter
        if (length(sensorData) > 1)
            sensorDataWaveletFiltered = [];
            for ax = 1:6
                %                 sensorDataWaveletFiltered(:,ax) = myDWaveletFilter(detrend(sensorData(:,ax), 0), FS , 0, myWavelet, thr_coeff, 'h', noiseSTD, boolLeveldep, tplevel,hplevel, false, false, false);
                sensorDataWaveletFiltered(:,ax) = wdenoise(detrend(sensorData(:,ax), 0), DenoisingMethod="UniversalThreshold");

            end

        end

        expression = '(\d+-){5}\d+';
        impulseSettings = regexp(datapath, expression, 'match');

        expression = '(WOUTstruct|Wstruct|Noise|SmallBall|BigBall)';
        fileType = regexp(datapath, expression, 'match');

        %% Save Data
        cellMat{fnum,1} = sensorData;
        cellMat{fnum,2} = sensorDataWaveletFiltered;
        cellMat{fnum,3} = FS;
        cellMat{fnum,4} = time;
        cellMat{fnum,5} = impulseSettings;
        cellMat{fnum,6} = fileType;
        %         %% Fusion Factor
        %         accZ = detrend(cellMat{fnum,1}(:,3),0);
        %         gyroX = detrend(cellMat{fnum,1}(:,4),0);
        %         winkelbeschl = [0; diff(gyroX)];
        % %         fusionFactor(fnum) = calcFactor(winkelbeschl, accZ,"raw")
    else
        disp("FILE NOT FOUND: " + path+ filenames(fnum))
    end
end


fileTypeTable = string(cellMat(:,6));
% fileTypeTable = fileTypeTable(:,6);


%% Resonance Mode Sensor Fusion
disp("VVVV################Sensor Fusion##################VVVV")
WITHINDEX = find(fileTypeTable == 'Wstruct', 1);
% WITHINDEX = find(fileTypeTable == 'WOUTstruct', 1);

% finde noise Matrix
NOISEINDEX = find(fileTypeTable == 'Noise' ,1);
% calculate noise estimation vom prior noise measurement
noiseFS = cellMat{NOISEINDEX,3};
noise = cellMat{NOISEINDEX,1}(:,1:6);
noise(:,4:6) = mydiff(noise(:,4:6))*noiseFS;
noiseMat = getNoisePerLvl(noise, noiseFS, 0,myWavelet,false);
% noiseMat= [];
if ~isempty(WITHINDEX)

    sampRate = cellMat{WITHINDEX,3};
    IMU = detrend(cellMat{WITHINDEX,1},6);
    IMU = [IMU(:,1:3), mydiff(IMU(:,4:6))*sampRate];
    SA_IMU = timetable(seconds(cellMat{WITHINDEX,4}/1000000), IMU);


%     resModeMat = cell(1,3);
%     resModeMat(1,:) = {"ax bis az", [0,0,1,0,0,0;0,0,0,0,0,0],30};
%     resModeMat = [resModeMat; {"aZ + wX", [0,0,1,0,0,0;0,0,0,0,0,0],100}];
%     resModeMat = [resModeMat; {"aZ + wX", [0,0,1,0,0,0;0,0,0,0,0,0],118}];
%     resModeMat = [resModeMat; {"aX", [1,0,0,0,0,0;0,0,0,0,0,0],423}];
%     resModeMat = [resModeMat; {"aX", [1,0,0,0,0,0;0,0,0,0,0,0],465}];
%     resModeMat = [resModeMat; {"aX", [0,0,0,0,1,0;0,0,0,0,0,0],618}];
%     resModeMat = [resModeMat; {"aY + aZ + wX", [0,1,0,0,0,0;0,0,0,0,0,0],1025}];

resModeMat = cell(1,3);
    resModeMat(1,:) = {"ax bis az", [0,0,1,0,0,0; 0,0,0,0,0,0],103};
%     resModeMat = [resModeMat; {"aZ + wX", [0,0,1,0,0,0; 0,0,0,0,0,0],106}];
%     resModeMat = [resModeMat; {"just aX", [1,0,0,0,0,0; 0,0,0,0,0,0],390}];
%     resModeMat = [resModeMat; {"just wY", [0,0, 0,0,1,0; 0,0,0,0,0,0],528}];
    resModeMat = [resModeMat; {"aY + wX", [0,1, 0,0,0,0; 0,0,0,0,0,0],1000}];

    disp("Noise Without SF")
    baseResonance = waveletResonanceSensorFusion(IMU, sampRate, resModeMat,0,myWavelet,[],noiseMat,thresholdFactor,false,false);
    baseResonanceWav = waveletResonanceSensorFusion(IMU, sampRate, resModeMat,0,myWavelet,[],noiseMat,thresholdFactor,true,false);
%     baseResonanceWav = wdenoise(baseResonance,DenoisingMethod="UniversalThreshold", ThresholdRule="Hard");

    noiseRes = waveletResonanceSensorFusion(noise, sampRate, resModeMat,0,myWavelet,[],[],0,false,false);

%     energySum = zeros(size(baseResonanceWav));
%     for i = 1:length(baseResonanceWav)
%         for u = 1:size(baseResonanceWav, 2)
%             if i ~= 1
%                 energySum(i,u) = abs(baseResonanceWav(i,u))^2 + energySum(i-1,u);
%             else
%                 energySum(i,u) = abs(baseResonanceWav(i,u))^2;
%             end
%         end
%     end

    SA_baseResonance = timetable(seconds(cellMat{WITHINDEX,4}/1000000), baseResonance, baseResonanceWav);

%     resModeMat = cell(1,3);
%     resModeMat(1,:) = {"ax bis az", [0,0,1,0,0,0; 0,0,0,0,0,0],30};
%     resModeMat = [resModeMat; {"aZ + wX", [0,0,1,-7.144e-05,0,0; 0,0,0,0,0,0],100}];
%     resModeMat = [resModeMat; {"aZ + wX", [0,0,1,-4.83e-05,0,0; 0,0,0,0,0,0],118}];
%     resModeMat = [resModeMat; {"aX", [1,0,0,0,0,0;0,0,0,0,0,0],423}];
%     resModeMat = [resModeMat; {"aX", [1,0,0,0,0,0;0,0,0,0,0,0],465}];
%     resModeMat = [resModeMat; {"aX", [0,0,0,0,1,0;0,0,0,0,0,0],618}];
%     resModeMat = [resModeMat; {"aY + wX", [0,1, 0,-5.74e-06,0,0; 0,0,0,0,0,0],1025}];

    resModeMat = cell(1,3);
    resModeMat(1,:) = {"aZ + wX", [0,0,1,-7.2e-05,0,0; 0,0,0,0,0,0],103};
%     resModeMat = [resModeMat; {"aZ + wX", [0,0,1,-7.2e-05,0,0; 0,0,0,0,0,0],106}];
%     resModeMat = [resModeMat; {"just aX", [1,0,0,0,0,0; 0,0,0,0,0,0],390}];
%     resModeMat = [resModeMat; {"just wY", [0,0, 0,0,1,0; 0,0,0,0,0,0],528}];
    resModeMat = [resModeMat; {"aX + wX", [0,1, 0,-6.8e-06,0,0; 0,0,0,0,0,0],1000}];

    disp("Noise With SF")
    testRes = waveletResonanceSensorFusion(IMU, sampRate, resModeMat,0,myWavelet,[],noiseMat,thresholdFactor,false,false);
    testResWav = waveletResonanceSensorFusion(IMU, sampRate, resModeMat,0,myWavelet,[],noiseMat,thresholdFactor,true,false);
%         testResWav = wdenoise(testRes,DenoisingMethod="UniversalThreshold", ThresholdRule="Hard");

    
    SA_testResonance = timetable(seconds(cellMat{WITHINDEX,4}/1000000), testRes, testResWav);

    SA_IMU = timetable(seconds(cellMat{WITHINDEX,4}/1000000), IMU);


    % plot time function
    graphnames=[];
    for i = 1:size(resModeMat, 1)
        graphnames = [graphnames, "Vibration Mode: " + resModeMat{i,3}];
    end

    numFigures = size(resModeMat,1);
    fig = figure("Name","Wavelet Resonance Fusion Modes");
        t = tiledlayout(numFigures,3);
    fig.Position = [39 158 1838 803];

%     t = tiledlayout(numFigures,2);

    t.Padding = 'compact';
    t.TileSpacing = 'compact';


    legendNames = ["Without Fusion","With Fusion"];
    for v = 1: numFigures
        tiles(v) = nexttile;
        hold on
        plot(cellMat{WITHINDEX,4}/1000000, baseResonance(:,v), 'LineWidth',0.7);
        plot(cellMat{WITHINDEX,4}/1000000, testRes(:,v), 'LineWidth',0.7);
        hold off
        grid on
        grid minor
        if(v ==1)
            legend(legendNames)
        end
        set(gca,'FontSize',16, 'FontName', 'Times')
        %         title("Signal of Mode "+ v +" at "+ resModeMat{v,3} + " Hz")
        title("Unfiltered Signal at "+ resModeMat{v,3} + " Hz")
        ylabel("Acceleration [g]")
        xlabel("Time [s]")
        xlim([30.42 32.23])

        tiles(v+numFigures) = nexttile;
        hold on
        plot(cellMat{WITHINDEX,4}/1000000, baseResonanceWav(:,v), 'LineWidth',0.7);
        plot(cellMat{WITHINDEX,4}/1000000, testResWav(:,v), 'LineWidth',0.7);
        hold off
        grid on
        grid minor
%         legend(legendNames)
        set(gca,'FontSize',16, 'FontName', 'Times')
%         title("Wavelet Filtered Signal of Mode "+ v +" at :" + resModeMat{v,3} + " Hz")
        title("Wavelet Filtered Signal at " + resModeMat{v,3} + " Hz")
        ylabel("Acceleration [g]")
        xlabel("Time [s]")
        xlim([30.42 32.23])

        tiles(v+2*numFigures) = nexttile;
        hold on
        plot(cellMat{WITHINDEX,4}/1000000, baseResonance(:,v)  - baseResonanceWav(:,v), 'LineWidth',0.7);
        plot(cellMat{WITHINDEX,4}/1000000, testRes(:,v)-testResWav(:,v), 'LineWidth',0.7);
        hold off
        grid on
        grid minor
%         legend(legendNames)
        set(gca,'FontSize',16, 'FontName', 'Times')
%         title("Noise Extracted from Signal of Mode "+ v +" at :" + resModeMat{v,3} + " Hz")
        title("Noise Signal at " + resModeMat{v,3} + " Hz")
        ylabel("Acceleration [g]")
        xlabel("Time [s]")
        xlim([30.42 32.23])
    end
    linkaxes(tiles, 'x');

    % Display noise Powers and Energies
    disp("BASE RESONANCE")
    for i = 1 : size(baseResonance,2)
        disp("SIGNAL Energy: Mode " + i + ": RAW " + rms(baseResonance(:,i))^2 +" g^2  WAV: "  + rms(baseResonanceWav(:,i))^2 +" g^2");
    end


    disp("")
    disp("FUSION RESULT")
    for i = 1 : size(baseResonance,2)
        disp("SIGNAL Energy: Mode " + i + ": RAW " + rms(testRes(:,i))^2 + " g^2  WAV: " + rms(testResWav(:,i))^2 +" g^2");
    end

    disp("")
    disp("ENERGY DIFFERENCE")
    disp("Energy(Fused) - Energy(Base)")
    for i = 1 : size(baseResonance,2)
        disp("SIGNAL Energy: Mode " + i + ": RAW " + (rms(testRes(:,i))^2 - rms(baseResonance(:,i))^2) + " %  WAV: " + ((rms(testResWav(:,i))^2 - rms(baseResonanceWav(:,i))^2)/rms(baseResonanceWav(:,i))^2*100) +" %");
    end

    disp("")
    disp("ESTIMATED NOISE")
%     disp("(Energy(Fused) - Energy(Base)) - (Energy(Fused_OnlySignal)-Energy(Base_OnlySignal)) / Energy(Noise) * 100%")
%     disp("((Energy Reduction of Overall Signal) - (Energy Reduction of Useful Signal)) / Energy of Noise * 100%")

    for i = 1 : size(baseResonance,2)
        estNoiseSTDBase = std(baseResonance(:,i) - baseResonanceWav(:,i));
        estNoiseSTDTest = std(testRes(:,i) - testResWav(:,i));
        disp("SIGNAL Energy: Mode " + i + ":  " + (rms(testRes(:,i) - testResWav(:,i))^2)  +" g^2  STD: " + estNoiseSTDTest +"  " +estNoiseSTDBase+ "  NR:" + (estNoiseSTDBase-estNoiseSTDTest) +" g | "+ (estNoiseSTDBase-estNoiseSTDTest)/estNoiseSTDBase*100 + "%");
    end
end

disp("##################################")

%% Save Data for Signal Analyzer

WITHINDEX = find(fileTypeTable == 'Wstruct', 1);
if ~isempty(WITHINDEX)
    % SA_WITH = timetable(seconds(cellMat{WITHINDEX,4}/1000000), detrend(cellMat{WITHINDEX,1}(:,:),0), detrend(cellMat{WITHINDEX,2}(:,:),0));
    SA_WITH = timetable(detrend(cellMat{WITHINDEX,1}(:,:),0), detrend(cellMat{WITHINDEX,2}(:,:),0), 'SampleRate', cellMat{WITHINDEX,3});
end
WOUTINDEX = find(fileTypeTable == 'WOUTstruct' ,1);
if ~isempty(WOUTINDEX)
    % SA_WITH = timetable(seconds(cellMat{WITHINDEX,4}/1000000), detrend(cellMat{WITHINDEX,1}(:,:),0), detrend(cellMat{WITHINDEX,2}(:,:),0));
    SA_WOUT = timetable(detrend(cellMat{WOUTINDEX,1}(:,:),0), detrend(cellMat{WOUTINDEX,2}(:,:),0), 'SampleRate', cellMat{WOUTINDEX,3});
end
NOISEINDEX = find(fileTypeTable == 'Noise' ,1);
if ~isempty(NOISEINDEX)
    % SA_WITH = timetable(seconds(cellMat{WITHINDEX,4}/1000000), detrend(cellMat{WITHINDEX,1}(:,:),0), detrend(cellMat{WITHINDEX,2}(:,:),0));
    SA_NOISE = timetable(detrend(cellMat{NOISEINDEX,1}(:,:),0), detrend(cellMat{NOISEINDEX,2}(:,:),0), 'SampleRate', cellMat{NOISEINDEX,3});
end

%% Display SNR in Relation to Wavelet Filtered ACCZ

WITHINDEX = find(fileTypeTable == 'Wstruct', 1);

WOUTINDEX = find(fileTypeTable == 'WOUTstruct' ,1);

NOISEINDEX = find(fileTypeTable == 'Noise' ,1);

if ~any([isempty(WITHINDEX), isempty(WOUTINDEX), isempty(NOISEINDEX)])

    WITHFS = cellMat{WITHINDEX,3};
    WOUTFS = cellMat{WOUTINDEX,3};
    NOISEFS = cellMat{NOISEINDEX,3};

    disp("Signal Power")
    disp("ONLY CAGE")
    disp("SIGNAL ENERGY: Acceleration in X:  " + rms(cellMat{WOUTINDEX,2}(:,1))^2 +" g^2");
    disp("SIGNAL ENERGY: Acceleration in Y:  " + rms(cellMat{WOUTINDEX,2}(:,2))^2 +" g^2");
    disp("SIGNAL ENERGY: Acceleration in Z:  " + rms(cellMat{WOUTINDEX,2}(:,3))^2 +" g^2");
    disp("SIGNAL ENERGY: Rotation Rate in X:  " + rms(cellMat{WOUTINDEX,2}(:,4))^2 +" g^2");
    disp("SIGNAL ENERGY: Rotation Rate in Y:  " + rms(cellMat{WOUTINDEX,2}(:,5))^2 +" g^2");
    disp("SIGNAL ENERGY: Rotation Rate in Z:  " + rms(cellMat{WOUTINDEX,2}(:,6))^2 +" g^2");
    pythagSum1 = sqrt(cellMat{WOUTINDEX,2}(:,1).^2 + cellMat{WOUTINDEX,2}(:,2).^2 + cellMat{WOUTINDEX,2}(:,3).^2);
    sumGyro1 = sum(cellMat{WOUTINDEX,2}(:,4:6),2);
    disp("SIGNAL ENERGY: SUM Acc:  " + (rms(pythagSum1)^2) +" g^2  SUM GYRO: " + (rms(sumGyro1)^2) + " g^2");

    disp("WITH STRUCTURE")
    disp("SIGNAL ENERGY: Acceleration in X:  " + rms(cellMat{WITHINDEX,2}(:,1))^2 +" g^2");
    disp("SIGNAL ENERGY: Acceleration in Y:  " + rms(cellMat{WITHINDEX,2}(:,2))^2 +" g^2");
    disp("SIGNAL ENERGY: Acceleration in Z:  " + rms(cellMat{WITHINDEX,2}(:,3))^2 +" g^2");
    pythagSum2 = sqrt(cellMat{WITHINDEX,2}(:,1).^2 + cellMat{WITHINDEX,2}(:,2).^2 + cellMat{WITHINDEX,2}(:,3).^2);
    sumGyro2 = sum(cellMat{WITHINDEX,2}(:,4:6),2);
    disp("SIGNAL ENERGY: SUM:  " + (rms(pythagSum2)^2) +" g^2  SUM GYRO: " + (rms(sumGyro2)^2) + " g^2");
    disp("SIGNAL ENERGY DIFF PERCENT: SUM:  " + (rms(pythagSum2)^2)/(rms(pythagSum1)^2)*100 +" %  SUM GYRO: " + (rms(sumGyro2)^2)/(rms(sumGyro1)^2)*100 + " %");
end
%% Analyze Peaks
% buffer = cell(2,1);
% bufferfs = [cellMat{2,3},cellMat{1,3}];%,cellMat{1,3},cellMat{1,3}];
% buffer{1} = cellMat{2,1}(:,3);
% buffer{2} = cellMat{1,1}(:,3);
% % buffer{3} = levelIndepFusion;
% % buffer{4} = levelDepFusion;
% rslt = analyseAmpRamps(buffer, bufferfs,2, 0.005, ["without MechAmp", "with MechAmp", "with SF", "with Wave SF"]);

%% Settings
frequencyRange = [0 1300];
filterSetting = [3 301]; % envelop in Hertz

%% SHOW SENSOR FUSION RESULTS

WITHINDEX = find(fileTypeTable == 'Wstruct', 1);

WOUTINDEX = find(fileTypeTable == 'WOUTstruct' ,1);

NOISEINDEX = find(fileTypeTable == 'Noise' ,1);

if ~any([isempty(WITHINDEX), isempty(WOUTINDEX), isempty(NOISEINDEX)]) && exist("accZ", "var")

    WITHFS = cellMat{WITHINDEX,3};
    WOUTFS = cellMat{WOUTINDEX,3};
    NOISEFS = cellMat{NOISEINDEX,3};

    buffer = cell(7,1);
    bufferfs = [WITHFS,WITHFS,WITHFS,WITHFS,WITHFS, WOUTFS,NOISEFS];
    buffer{1} = accZ;
    buffer{2} = pseudoAccZ;
    buffer{3} = pseudoAccZ1;
    buffer{4} = TwoXFusion;
    buffer{5} = ThreeXFusion;
    buffer{6} = detrend(cellMat{WOUTINDEX,1}(:,3),0);
    buffer{7} = detrend(cellMat{NOISEINDEX,1}(:,3),0);

    plotNames = ["RAW AccZ", "pseudo ACCZ" ,"pseudo ACCZ1", "2 X Fusion","3 X Fusion","without Beam","RAW Sensor Noise"];
    % plotBodeDiagram(buffer, bufferfs,"With Sensor Fusion vs Without",["Without Sensorfusion", "W-SF-Level-Independent", "W-SF-Level-Dependent"], [0,550], true,filterSetting, false,true);
    % plotfftdifference(buffer, bufferfs,"Spectrum Comparison of Sensor Fusion",plotNames, [0,300], true,filterSetting, false, false, false);
    showMultiAxisFFT(buffer, bufferfs, "Spectrum Comparison of Sensor Fusion", plotNames,["Frequency [Hz]","Amplitude [g] or [dps]"],frequencyRange,[], [3 301]);

    % Draw TIME SIGNALS of SENSOR FUSION
    if (false)
        figure("Name","Wavelet Filtered");
        hold on
        plot(cellMat{WITHINDEX,4}/1000000, cellMat{WITHINDEX,2}(:,3), 'LineWidth',1);
        plot(cellMat{WITHINDEX,4}/1000000, TwoXFusionWav, 'LineWidth',1);
        plot(cellMat{WITHINDEX,4}/1000000, ThreeXFusionWav, 'LineWidth',1);
        hold off
        grid
        legend(["WAV ACCZ", "WAV 2X FUSION", "WAV 3X FUSION"])
    end
end
%% SHOW WITH VS WITHOUT RESULTs

WITHINDEX = find(fileTypeTable == 'Wstruct', 1);

WOUTINDEX = find(fileTypeTable == 'WOUTstruct' ,1);

NOISEINDEX = find(fileTypeTable == 'Noise' ,1);

if ~any([isempty(WITHINDEX), isempty(WOUTINDEX), isempty(NOISEINDEX)])

    WITHFS = cellMat{WITHINDEX,3};
    WOUTFS = cellMat{WOUTINDEX,3};
    NOISEFS = cellMat{NOISEINDEX,3};
    graphnames = ["IMU Directly Mounted", "IMU With Beam Structure"];

%     buffer = cell(2,1);
%     bufferfs = [WOUTFS, WITHFS,NOISEFS];
%     buffer{1} = sum(cellMat{WOUTINDEX,1}(:,1:3),2);
%     buffer{2} = sum(cellMat{WITHINDEX,1}(:,1:3),2);
%     buffer{3} = detrend(sum(cellMat{NOISEINDEX,1}(:,1:3),2),0);
%     graphnames = ["IMU Directly Mounted", "IMU With Beam Structure"];

%     plotfftdifference(buffer, bufferfs,"Spectrum Difference AccZ with vs without mechanical amplification",graphnames, frequencyRange, true ,filterSetting, false,true, true);

    % Draw TIME SIGNALS
    if (true)
        fig = figure("Name","Structural Vibration with Force Impulses up to 2.16 mN s");
        hold on
%         plot(cellMat{WOUTINDEX,4}/1000000, sqrt(sum(cellMat{WOUTINDEX,2}(:,1:3).^2,2)), 'LineWidth',1);
%         plot(cellMat{WITHINDEX,4}/1000000, sqrt(sum(cellMat{WITHINDEX,2}(:,1:3).^2,2)), 'LineWidth',1);
        plot(cellMat{WOUTINDEX,4}/1000000, sum(cellMat{WOUTINDEX,2}(:,1:3),2));
        plot(cellMat{WITHINDEX,4}/1000000, sum(cellMat{WITHINDEX,2}(:,1:3),2));
        hold off
        grid on
        grid minor
        xlim([0 60])
        legend(graphnames)
        ylabel("Acceleration [g]")
        xlabel("Time [s]")
        set(gca,'FontSize',16, 'FontName', 'Times New Roman')
    end
%     exportgraphics(fig, 'output.pdf', 'ContentType', 'vector');
%     %% Calculate RMS and Energy per Impulse
% 
%     WITHINDEX = find(fileTypeTable == 'Wstruct', 1);
% 
%     WOUTINDEX = find(fileTypeTable == 'WOUTstruct' ,1);
% 
%     NOISEINDEX = find(fileTypeTable == 'Noise' ,1);
%     data1 = detrend(sum(cellMat{WOUTINDEX,1}(:,1:3),2));
%     data1FS = cellMat{WOUTINDEX,3};
%     data2 = detrend(sum(cellMat{WITHINDEX,1}(:,1:3),2));
%     data2FS = cellMat{WITHINDEX,3};
% 
%     interval = 0.5; %interval in seconds
%     totalsteps = min([(floor(length(data1)/floor(data1FS*interval))-1),(floor(length(data2)/floor(data2FS*interval))-1)]);
%     energies = zeros(totalsteps,2);
%     RMS = zeros(totalsteps,2);
% 
%     
%     [RMS(:,1), energies(:,1)] =  calculateImpulseMag(data1, interval, data1FS,totalsteps,1);
%     [RMS(:,2), energies(:,2)] =  calculateImpulseMag(data2, interval, data2FS,totalsteps,1);
% 
% %     excludes = [1,3,5,7,9,10,11];
%     excludes = [1:1:10];
%     RMS(excludes,:) = [];
%     energies(excludes,:) = [];
% 
%     disp("Linearity Force Impulse: " )
%     if (true)
%         figure("Name","Energy and RMS per Impulse");
%         t = tiledlayout(2,1);
% 
%         t.Padding = 'compact';
%         t.TileSpacing = 'compact';
% 
%         nexttile
%         title("Signal Energy per Impulse");
%         hold on
%         %         bar(linspace(0,totalsteps*interval ,totalsteps),[energies(:,1),energies(:,2), energies(:,2)./energies(:,1)], 'LineWidth',1);
%         bar([energies(:,1),energies(:,2), energies(:,2)./energies(:,1)], 'LineWidth',1);
%         hold off
%         grid
%         legend(["Without Structure", "With Structure", "Amplification"])
% 
%         nexttile
%         title("RMS per Impulse");
%         hold on
%         %         bar(linspace(0,totalsteps*interval ,totalsteps),[RMS(:,1),RMS(:,2), RMS(:,2)./RMS(:,1)], 'LineWidth',1);
%         bar([RMS(:,1),RMS(:,2), RMS(:,2)./RMS(:,1)], 'LineWidth',1);
%         hold off
%         grid
%         legend(["Without Structure", "With Structure", "Amplification"])
%     end
end

%% FFT WITH VS WITHOUT

WITHINDEX = find(fileTypeTable == 'Wstruct', 1);

WOUTINDEX = find(fileTypeTable == 'WOUTstruct' ,1);

NOISEINDEX = find(fileTypeTable == 'Noise' ,1);

if ~any([isempty(WITHINDEX), isempty(WOUTINDEX), isempty(NOISEINDEX)])
    WITHFS = cellMat{WITHINDEX,3};
    WOUTFS = cellMat{WOUTINDEX,3};
    NOISEFS = cellMat{NOISEINDEX,3};

    buffer = cell(3,1);
    bufferfs = [WOUTFS,WITHFS,NOISEFS];
    buffer{1} = detrend(cellMat{WOUTINDEX,1}(:,[1,4,2,5,3,6]),0);
    buffer{2} = detrend(cellMat{WITHINDEX,1}(:,[1,4,2,5,3,6]),0);
    buffer{3} = detrend(cellMat{NOISEINDEX,1}(:,[1,4,2,5,3,6]),0);

    AxisNames = ["Acceleration X Axis","Gyro X Axis","Acceleration Y Axis","Gyro Y Axis","Acceleration Z Axis",  "Gyro Z Axis"];

    xlabels = ["Frequency [Hz]","Frequency [Hz]","Frequency [Hz]","Frequency [Hz]","Frequency [Hz]","Frequency [Hz]"];
%     ylabels = ["Acceleration [m/s^2]","Angular Rate [degree/s]","Acceleration [m/s^2]","Angular Rate [degree/s]","Acceleration [m/s^2]",  "Angular Rate [degree/s]"];
    ylabels = ["Acceleration [g]","Angular Rate [degree/s]","Acceleration [g]","Angular Rate [degree/s]","Acceleration [g]",  "Angular Rate [degree/s]"];
    customMultiAxisFFT(buffer, bufferfs, "FFT of Cage Vibration With VS Without Structure", ["IMU Without Beam","IMU With Beam", "Raw Sensor Noise"],AxisNames,xlabels, ylabels,frequencyRange,[], [1 1001]);





end

%% test
% WITHINDEX = find(fileTypeTable == 'Wstruct', 1);
% if (~isempty(WITHINDEX) && true)
%         WITHFS = cellMat{WITHINDEX,3};
%     buffer = cell(1,1);
%     bufferfs = [WITHFS];
%     buffer{1} = detrend([cellMat{WITHINDEX,1}(:,1:3),mydiff(cellMat{WITHINDEX,1}(:,4:6)).*WITHFS],0);
% 
%     AxisNames = ["Acceleration X Axis","Acceleration Y Axis","Acceleration Z Axis", "Gyro X Axis", "Gyro Y Axis", "Gyro Z Axis"];
%     xlabels = ["f [Hz]","f [Hz]","f [Hz]","f [Hz]","f [Hz]","f [Hz]"];
%     ylabels = ["Acceleration [m/s^2]","Acceleration [m/s^2]","Acceleration [m/s^2]", "Angular Rate [degree/s]", "Angular Rate [degree/s]", "Angular Rate [degree/s]"];
%     showMultiAxisFFT(buffer, bufferfs, "FFT of Cage Vibration With VS Without Structure", ["Without Structure"],AxisNames,xlabels, ylabels,frequencyRange,[], [3 301]);
% end


%% OUTDATED:
% % Level Dependent Wavelet Sensor Fusion
% % this sensor fusion is not applicable anymore
% if (false)
%     WITHINDEX = find(fileTypeTable == 'Wstruct', 1);
% 
%     if ~isempty(WITHINDEX) && false
% 
%         %calculate pseudo accZ from gyroX
%         accZ = detrend(cellMat{WITHINDEX,1}(:,3),0);
%         accY = detrend(cellMat{WITHINDEX,1}(:,2),0);
%         gyroX = detrend(cellMat{WITHINDEX,1}(:,4),0);
%         winkelbeschl = mydiff(gyroX)*cellMat{2,3}; %df/dt
%         WITHFS = cellMat{WITHINDEX,3};
% 
%         fRANGE = [90, 110];
%         % disp("scale factor:    " + calcFactor(winkelbeschl, accZ, "wavelet"));
%         % fusionFactor = calcFactor(detrend(winkelbeschl,0), detrend(accZ,0), "FFT", WITHFS, fRANGE);
%         fusionFactor = -6.8e-05;%-1.25e-04; Fusion Factor for dF/dt(GyroZ) to AccZ U0327
%         pseudoAccZ = detrend(fusionFactor * winkelbeschl,0);
%         % fusionFactor1 = calcFactor(detrend(accY,0), detrend(accZ,0), "FFT", WITHFS, fRANGE);
%         fusionFactor1 = 7.4838;% Fusion Factor AccY to AccZ U0327
%         pseudoAccZ1 = detrend(fusionFactor1 * accY,0);
% 
%         frequencyFilter = 0;
%         windFactor = WITHFS/length(accZ);
%         TwoXFusion = waveletSensorFusion([accZ, pseudoAccZ],WITHFS,0,myWavelet,'h', frequencyFilter,0, false, true, windFactor, false);
%         ThreeXFusion = waveletSensorFusion([accZ, pseudoAccZ,pseudoAccZ1],WITHFS,0,myWavelet,'h', frequencyFilter,0, false, true, windFactor, false);
%         TwoXFusionWav = waveletSensorFusion([accZ, pseudoAccZ],WITHFS,0,myWavelet,'h', frequencyFilter,0, true, true, windFactor, false);
%         ThreeXFusionWav = waveletSensorFusion([accZ, pseudoAccZ,pseudoAccZ1],WITHFS,0,myWavelet,'h', frequencyFilter,0, true, true, windFactor, false);
% 
% 
%         % waveletFiltSFAverage = waveletSensorFusion([accZ, gyroX],cellMat{1,3},0,myWavelet,'h', 0,0, true, false, false, false);
%         experimentalData =  [pseudoAccZ, TwoXFusion,ThreeXFusion];
%         experimentalData(:,4) = myDWaveletFilter(detrend(accZ, 0), WITHFS , 0, myWavelet, thr_coeff, 'h', noiseSTD, boolLeveldep, tplevel,hplevel, false, false, false);
%         experimentalData(:,5) = TwoXFusionWav;
%         experimentalData(:,6) = ThreeXFusionWav;
% 
%         SA_SF_DATA = timetable(seconds(cellMat{WITHINDEX,4}/1000000), experimentalData);
%     end
% 
% end


% Function to draw the fft of multiple sources in 3 axis, 6axis?
% input paramters: data cell array, FS vector, title, legend, AxisNames, X
% range, yrange, Filter

function rslt = customMultiAxisFFT(data, FS, titlename, dataNames, AxisNames,xlabels, ylabels,  xRange, yRange, filter)
if isempty(data)
    disp("ERROR: showMultiAxisFFT: Data Input Parameter is Empty")
end

numAxis = size(data{1},2);
numColumns = ceil(1+(numAxis-3)/3);
numSources = size(data,1);

fig = figure("Name",titlename);
fig.Position = [81 122 1791 706];

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
    if (AxIndex==6)
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

    if(any(AxIndex == [1,3,5]))
    ylim([0 4.5e-4])
    else
    ylim([0 0.01])
    end
end
linkaxes(tiles, 'x' )
exportgraphics(fig, 'spectrum experimental measurement_CenterImpulse.pdf', 'ContentType', 'vector');
end
