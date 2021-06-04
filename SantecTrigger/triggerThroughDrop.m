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
triggerStep = 1 / 2;
estimatedFSR = 0.8;
triggerThreshold = 0.1;

% Calibration with trigger
[throughPks, throughPksLocs] = findtriggers(throughTrigger, triggerThreshold, triggerStep/wStep);

% Plot
testPks = figure();
hold on;
plot(initWavelength, throughTrigger)
plot(initWavelength(throughPksLocs), throughTrigger(throughPksLocs), '*')

% assist funcs
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
            