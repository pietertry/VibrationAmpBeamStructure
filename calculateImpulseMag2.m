function [RMS, ENERGY] = calculateImpulseMag2(time, data, interval, totalsteps, ratio)
% N = length(data);
% data = detrend(data);

if isempty(ratio) 
    ratio = 1;end

runningQuadrSum = zeros(size(data));

runningQuadrSum(1,:) = data(1,:).^2;
for i = 2: length(data)
    runningQuadrSum(i,:) = runningQuadrSum(i-1,:) +  data(i,:).^2;
end

if isempty(totalsteps)
    totalsteps = floor(time(end)/(interval*1e6))-1;
end

ENERGY = zeros(totalsteps,size(data,2));
RMS = zeros(totalsteps,size(data,2));

% figure;
for i = 1:totalsteps
    [~,stepStartIndex] = min( abs( time-((i-1)*interval *1e6)));
    stepStartIndex = stepStartIndex+1;
    [~,stepStopIndex] = min( abs( time-(((i-1)+ratio)*interval *1e6)));      

    RMS(i,:) = sqrt((runningQuadrSum(stepStopIndex,:) - runningQuadrSum(stepStartIndex,:)) / (stepStopIndex-stepStartIndex));%(interval* 1/FS));
    ENERGY(i,:) = runningQuadrSum(stepStopIndex,:) - runningQuadrSum(stepStartIndex,:);
end
end