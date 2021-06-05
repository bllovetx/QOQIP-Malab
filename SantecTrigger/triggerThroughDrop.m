clear; clc;

% import data
RepositoryPath = "E:\WorkspaceZXZQ\QOQIP-Malab\";
curDataPath = RepositoryPath + "SantecTrigger\data\";
throughData = importdata(curDataPath + "TestTriggerStep1nm100kHz_1548-1555.dat");
dropData = importdata(curDataPath + "TestTriggerStep1nm100kHz_1548-1555_Drop.dat");
% Read data
throughPower = throughData(:, 1);
throughTrigger = throughData(:, 2);
dropPower = dropData(:, 1);
dropTrigger = dropData(:, 2);

% set parameters
% (nm)
wStep = 1/100000;    % 100k sampling rate
wStart= 1548;
wStop = 1555;
initWavelength = wStart: wStep: wStop + 10 - wStep;
minTriggerStep = 1 / 2;
triggerStep = 1;
estimatedFSR = 0.8;
triggerThreshold = 0.1;

% Calibration with trigger
[throughPks, throughPksLocs] = findtriggers(throughTrigger, triggerThreshold, minTriggerStep/wStep);
[throughWavelength, throughPower] = wavelengthCalibration(throughPower, throughPksLocs, wStart, triggerStep, wStop);
[dropPks, dropPksLoc] = findtriggers(dropTrigger, triggerThreshold, minTriggerStep/wStep);
[dropWavelength, dropPower] = wavelengthCalibration(dropPower, dropPksLoc, wStart, triggerStep, wStop);

% Plot
throughDropFig = figure();
plot(throughWavelength, throughPower);
hold on;
plot(dropWavelength, dropPower);

%% Plot test
%testPks = figure();
%hold on;
%plot(initWavelength, throughTrigger)
%plot(initWavelength(throughPksLocs), throughTrigger(throughPksLocs), '*')

% assist funcs
% find trigger
function [pks, locs] = findtriggers(triggerData, threshold, minTriggerStep)
    % minTriggerStep : index
    myminTriggerStep = round(minTriggerStep);
    currentIndex = 1;
    pks = [];
    locs = [];
    while currentIndex <= length(triggerData)
        currentData = triggerData(currentIndex);
        if currentData > threshold
            locs = [locs, currentIndex];
            pks = [pks, currentData];
            currentIndex = currentIndex + myminTriggerStep;
        else
            currentIndex = currentIndex + 1;
        end
    end
end

% calibration with trigger
function [newWavelength, newData] = wavelengthCalibration(oldData, triggersInd, trigWavelengthStart, trigWavelengthStep, trigWavelengthStop)
    trigNum = length(triggersInd);
    assert(trigNum > 1, "Wrong Trigger Index");
    if size(oldData, 1) > 1
        oldData = transpose(oldData);
    end
    assert(size(oldData, 1) <= 1, "old data size error");
    lastTrig = 1;
    currentTrig = 2;
    trigWavelength = trigWavelengthStart: trigWavelengthStep:trigWavelengthStop;
    assert(length(trigWavelength) == trigNum, "missing trigger");
    newWavelength = [];
    newData = [];
    for i = 1:(trigNum -1)
        curIntervalNum = triggersInd(i+1) - triggersInd(i);
        curIntervalStart = trigWavelength(i);
        curIntervalStop  = trigWavelength(i+1);
        curIntervalStep  = (curIntervalStop - curIntervalStart) / curIntervalNum;
        newWavelength = [newWavelength, curIntervalStart: curIntervalStep: curIntervalStop - curIntervalStep];
        newData = [newData, oldData( triggersInd(i):(triggersInd(i+1)-1) )];
    end
end
