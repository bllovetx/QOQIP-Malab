% init
clear;clc;

% Global Path
% Notes: 1. output folder will be created if not exist, output file will be
%        overwritten if exist
pathToWorkSpace = "E:/WorkspaceZXZQ/";
global pathToData pathToOutput;
pathToData = pathToWorkSpace + "data_processing/AMF_Rings/data_input/";% inputPath
pathToOutput = pathToWorkSpace + "data_processing/AMF_Rings/output/";  % outputPath

% Pump Data (0.1mW)
% Notes: 1. pump should be calibrated using same gain
%        2. if zeroPointData is not needed, overwrite it to
%        zeros(length(importdata(pumpFile)))
dataPower = 10; %mW
pumpPower = 0.1;%mW
pumpFile = pathToData + "pump/" + "calibration_0_1mW_Gain10.dat";
closedPumpFile = pathToData + "pump/" + "calibration_close_Gain10.dat";
global pumpData zeroPointData;
zeroPointData = importdata( closedPumpFile );
pumpData = importdata( pumpFile )*(dataPower/pumpPower) - zeroPointData;

% Ring Data
% Notes:1. all input data should have same input power(set as dataPower)
%       2. if (both Through and Drop ports)exist, Through & Drop ports should
%       have same length with length of GapDict{ii}
%       3. typeAbbr, typeName, prefixDict, GapDict, EstimatedPara should
%       have same length
%       4. filename should be arranged in the same order with gap
gap = ["200nm", "300nm", "400nm", "500nm", "600nm", "700nm"];
global typeAbbr typeName EstimatedPara;
typeAbbr = ["DBC", "DBT", "SBC", "SBT", "LBT"];
typeName = ["Double Bus Circle ",...,
            "Double Bus Track " ,...,
            "Single Bus Circle ",...,
            "Single Bus Track " ,...,
            "Single Bus Long Track "...,
];
DBC_Through =  ["DoubleBusCircle200nm_In02Through_Gain10_1569_67_z08mum",...,
                "DoubleBusCircle300nm_In06Through_Gain10_1570_z08mum",...,
                "DoubleBusCircle400nm_In10Through_Gain10_1570_z08mum",...,
                "DoubleBusCircle500nm_In02Through_Gain10_1570_z10mum",...,
                "DoubleBusCircle600nm_In06Through_Gain10_1570_z10mum",...,
                "DoubleBusCircle700nm_In10Through_Gain10_1570_z10mum_02"...,
];
DBC_Drop = ["DoubleBusCircle200nm_In02Drop_Gain10_1570_z08mum",...,
            "DoubleBusCircle300nm_In06Drop_Gain10_1570_z08mum",...,
            "DoubleBusCircle400nm_In10Drop_Gain10_1569_45_z08mum",...,
            "DoubleBusCircle500nm_In02Drop_Gain10_1570_z10mum",...,
            "DoubleBusCircle600nm_In06Drop_Gain10_1569_95_z10mum",...,
            "DoubleBusCircle700nm_In10Drop_Gain10_1569_82_z10mum_02"...,
];
DBT_Through =  ["DoubleBusTrack200nm_In02Through_Gain10_1570_z08mum",...,
                "DoubleBusTrack300nm_In06Through_Gain10_1570_z08mum",...,
                "DoubleBusTrack400nm_In10Through_Gain10_1570_z08mum",...,
                "DoubleBusTrack500nm_In02Through_Gain10_1570_z08mum",...,
                "DoubleBusTrack600nm_In06Through_Gain10_1570_z08mum",...,
                "DoubleBusTrack700nm_In10Through_Gain10_1570_z08mum"...,
];
DBT_Drop = ["DoubleBusTrack200nm_In02Drop_Gain10_1570_z08mum",...,
            "DoubleBusTrack300nm_In06Drop_Gain10_1569_53_z08mum",...,
            "DoubleBusTrack400nm_In10Drop_Gain10_1569_55_z08mum",...,
            "DoubleBusTrack500nm_In02Drop_Gain10_1569_43_z08mum",...,
            "DoubleBusTrack600nm_In06Drop_Gain10_1569_37_z08mum",...,
            "DoubleBusTrack700nm_In10Drop_Gain10_1569_77_z08mum"...,
];
SBC =  ["SingleBusCircle200nm_In01_Gain10_1570_z08mum",...,
        "SingleBusCircle300nm_In03_Gain10_1570_z08mum",...,
        "SingleBusCircle400nm_In05_Gain10_1570_z08mum",...,
        "SingleBusCircle500nm_In01_Gain10_1570_z09mum",...,
        "SingleBusCircle600nm_In03_Gain10_1570_z09mum",...,
        "SingleBusCircle700nm_In05_Gain10_1570_z09mum"...,
];
SBT =  ["SingleBusTrack200nm_In01Through_Gain10_1570_z09mum",...,
        "SingleBusTrack300nm_In03Through_Gain10_1570_z09mum",...,
        "SingleBusTrack400nm_In05Through_Gain10_1570_z09mum",...,
        "SingleBusTrack500nm_In07Through_Gain10_1570_z09mum",...,
        "SingleBusTrack600nm_In09Through_Gain10_1570_z09mum",...,
        "SingleBusTrack700nm_In11Through_Gain10_1570_z09mum"...,
];
LBT =  ["SingleBusLongTrack300nm_In01Through_Gain10_1570_z06mum",...,
        "SingleBusLongTrack400nm_In03Through_Gain10_1570_z08mum",...,
        "SingleBusLongTrack500nm_In01Through_Gain10_1570_z08mum",...,
        "SingleBusLongTrack600nm_In03Through_Gain10_1570_z08mum"...,
];
PrefixDict = { {DBC_Through, DBC_Drop}, {DBT_Through, DBT_Drop}, ...,
                {SBC, []}, {SBT, []}, {LBT, []} ...,
};
GapDict = { gap, gap, gap, gap, gap(2:5) };
% EstimatedPara: [Delta_lambda/nm, FWHM/nm]
EstimatedPara = { [0.8, 0.02], [0.8, 0.02], [0.8, 0.02], [0.8, 0.02], [0.16, 0.02] };


% args
global wavelength wavelengthStep;
wavelengthStart = 1500;     %nm
wavelengthEnd = 1630;       %nm
wavelengthStep = 100/1000000;%nm
wavelength = (wavelengthStart+wavelengthStep:wavelengthStep:wavelengthEnd)';%nm

% precalculate
groupNum = 0;% used for waitbar
for ii = 1:length(typeAbbr)
    groupNum = groupNum + length(PrefixDict{ii}{1}) + length(PrefixDict{ii}{2});
end

% process all
curNum = 0;
myWaitBar = waitbar(curNum/groupNum, "Processing Data");
for ii = 1:length(typeAbbr)
    iiDict = PrefixDict{ii};
    assert( length(iiDict) == 2, "Each PrefixDict should have two ports(Through/Drop), if not needed input empty array");
    iiThroughDict = iiDict{1};
    iiDropDict = iiDict{2};
    iiGapDict = GapDict{ii};
    iiGapNum = length(iiGapDict);
    Q_average = [];
    if ~isempty(iiThroughDict)   % has Through Port
        assert(length(iiThroughDict)==iiGapNum, typeAbbr(ii) + " Through port number conflict with gap number");
        Q_average_through = zeros(iiGapNum, 1);
        for jj = 1:iiGapNum
            Q_average_through(jj) = DropThoughProcessor(iiThroughDict(jj), "Through", ii, iiGapDict(jj));
            curNum = curNum + 1;
            waitbar(curNum/groupNum, myWaitBar, "Processing Data(" + sprintf('%2.2f', 100*curNum/groupNum) + "%)");
        end
        Q_average = [Q_average, Q_average_through];
    end
    if ~isempty(iiDropDict) % has Drop Port
        assert(length(iiDropDict)==iiGapNum, typeAbbr(ii) + " Drop port number conflict with gap number");
        Q_average_drop = zeros(iiGapNum, 1);
        for jj = 1:iiGapNum
            Q_average_drop(jj) = DropThoughProcessor(iiDropDict(jj), "Drop", ii, iiGapDict(jj));
            curNum = curNum + 1;
            waitbar(curNum/groupNum, myWaitBar, "Processing Data(" + sprintf('%2.2f', 100*curNum/groupNum) + "%)");
        end
        Q_average = [Q_average, Q_average_drop];
    end
    save(pathToOutput + typeAbbr(ii) + "/" + typeAbbr(ii) + "_Q_ave.dat",   'Q_average', '-ascii', '-double', '-tabs');
end
close(myWaitBar);


function Q_factor_ave = DropThoughProcessor(dataFilePrefix, DropThroughProp, RingTypeIndex, GapStr)
    % dataFilePrefix example: "SingleBusCircle200nm_In01_Gain10_1570_z08mum"
    % DropThroughProp: {"Drop", "Through"}
    % GapStr: "***nm"
    % RingTypeIndex: 1~5 ["DBC", "DBT", "SBC", "SBT", "LBT"]
    % !! Noting that this func call saveas, which will **overwrite** if file already exist
    % fetch data (global vars can convert to input parameters)
    global pumpData zeroPointData wavelength wavelengthStep;
    global typeAbbr typeName pathToData pathToOutput EstimatedPara;
    peakSearchPeriod = 0.9*( (EstimatedPara{RingTypeIndex}(1)) / wavelengthStep );
    fittingSamplesHalf = peakSearchPeriod/2;
    initEstimatedSigma = EstimatedPara{RingTypeIndex}(2);
    curInputPath = pathToData + typeAbbr(RingTypeIndex) + "/";
    tempData = importdata( curInputPath + dataFilePrefix + ".dat" ) - zeroPointData;
    % Settings
    curOutputPath = pathToOutput + typeAbbr(RingTypeIndex) + "/";
    if ~exist(curOutputPath, 'dir')
        mkdir(curOutputPath);
    end
    curOutputPrefix = curOutputPath + dataFilePrefix + "_";
    logPath = curOutputPath + "logfiles/";
    if ~exist(logPath, 'dir')
        mkdir(logPath);
    end
    logPrefix = logPath + dataFilePrefix + "_log_";
    % figure('visible','off');
    % set(gcf, 'position', [50, 50, 1500, 450]);
    % plot(wavelength, tempData./pumpData, 'LineWidth', 1)
    % xlabel('Wavelength [nm]', 'FontSize', 18, 'FontName', 'Arial');
    % ylabel('Transmission', 'FontSize', 18, 'FontName', 'Arial');
    % title(typeName(RingTypeIndex) + DropThroughProp + " Spectral Transimission" + " gap=" + GapStr, 'FontSize', 18, 'FontName', 'Arial');

    %% Find peaks
    % search peaks
    power = (tempData./pumpData)/(max(tempData./pumpData));
    if DropThroughProp == "Through"
        findpeakTemp = 1-power;
    elseif DropThroughProp == "Drop"
        findpeakTemp = power;
    else
        assert(false, "Wrong DropThroughProp");
    end
    [pks,locs] = findpeaks(findpeakTemp,'minpeakdistance',peakSearchPeriod);
    resonance_wavelength_pks = wavelength(locs);
    % wavelength_diff_pks = abs(diff(resonance_wavelength_pks));
    resonance_frequency_pks = 299792458./resonance_wavelength_pks;% GHz
    frequency_diff_pks = abs(diff(resonance_frequency_pks));%GHz
    % frequency_diff_pks_ave = mean(frequency_diff_pks);
    % Relative_mode_index = (1:length(resonance_frequency_pks))';
    % Plot peaks
    normalizedTransmissionFig = figure('visible','off');
    set(gcf, 'position', [50, 50, 1500, 450]);
    plot(wavelength, power, '.-', 'Color', [65,105,225]/255, 'MarkerSize', 4, 'LineWidth', 1);
    xlabel('Wavelength [nm]', 'FontSize', 24, 'FontName', 'Arial');
    ylabel('Normalized transmission', 'FontSize', 24, 'FontName', 'Arial');
    hold on
    plot(resonance_wavelength_pks, power(locs),'*','color','r');
    title(typeName(RingTypeIndex) + DropThroughProp + " Normalized Spectral Transimission" + " gap=" + GapStr, 'FontSize', 18, 'FontName', 'Arial');

    %% Fitting peaks
    % define fitshape
    myfittype1 = fittype('A + B*(sigma^2)/(4*(x-x0)^2+sigma^2)', 'dependent',{'y'},...,
        'independent', {'x'}, 'coefficients',{'A', 'B', 'sigma', 'x0'}...,
    );
    % estimate sigma
    sigmaEstimation = zeros(length(resonance_wavelength_pks)-2, 1);
    for ii=1:(length(resonance_wavelength_pks)-2)
        x = wavelength(locs(ii+1)-fittingSamplesHalf:locs(ii+1)+fittingSamplesHalf);
        y = power(locs(ii+1)-fittingSamplesHalf:locs(ii+1)+fittingSamplesHalf);
        if DropThroughProp == "Through"
            estimatedAB = [max(y), -max(y)];
        elseif DropThroughProp == "Drop"
            estimatedAB = [0, max(y)];
        end
        myfit_1 = fit(x, y, myfittype1, 'Start', [estimatedAB, initEstimatedSigma, resonance_wavelength_pks(ii+1)]);
        sigmaEstimation(ii) = abs(myfit_1.sigma);
    end
    estimatedSigma = mean(sigmaEstimation);
    
    % declare variables
    resonance_wavelength_fit = zeros(length(resonance_wavelength_pks)-2, 1);
    FWHM = zeros(length(resonance_wavelength_pks)-2, 1);
    FWHM_mu = zeros(length(resonance_wavelength_pks)-2, 1);
    Q_factor = zeros(length(resonance_wavelength_pks)-2, 1);
    Extinction_ratio = zeros(length(resonance_wavelength_pks)-2, 1);
    tau = zeros(length(resonance_wavelength_pks)-2, 1);
    % fit & plot fitting progress
    fitFig = figure('visible','off');
    for ii=1:(length(resonance_wavelength_pks)-2)
        x = wavelength(locs(ii+1)-fittingSamplesHalf:locs(ii+1)+fittingSamplesHalf);
        y = power(locs(ii+1)-fittingSamplesHalf:locs(ii+1)+fittingSamplesHalf);
        if DropThroughProp == "Through"
            estimatedAB = [max(y), -max(y)];
        elseif DropThroughProp == "Drop"
            estimatedAB = [0, max(y)];
        end
        myfit_1 = fit(x, y, myfittype1, 'Start', [estimatedAB, estimatedSigma, resonance_wavelength_pks(ii+1)]);
        
        resonance_wavelength_fit(ii) = myfit_1.x0;  % nm
        FWHM(ii) = abs(myfit_1.sigma*1000);              % pm
        Q_factor(ii) = myfit_1.x0/FWHM(ii)*1000;    
        FWHM_mu(ii) = 299792458*FWHM(ii)/(myfit_1.x0)^2; 
        if DropThroughProp == "Through"
            Extinction_ratio(ii) = abs(10*log10(min(y)/myfit_1.A));
        elseif DropThroughProp == "Drop"
            Extinction_ratio(ii) = abs(10*log10(myfit_1.A/max(y)));
        end
        Extinction_ratio(ii) = abs(10*log10(min(y)/myfit_1.A));
        tau(ii) = 1/(2*pi*FWHM_mu(ii))*1000000;
        
        plot(myfit_1, x, y)%, 'b.', 'MarkerSize', 14
        hold on
        legend off
    end
    title(typeName(RingTypeIndex) + DropThroughProp + " Fitting Progress" + " gap=" + GapStr, 'FontSize', 18, 'FontName', 'Arial');

    resonance_frequency_fit = 299792458./resonance_wavelength_fit;
    frequency_diff_fit = abs(diff(resonance_frequency_fit));
    % frequency_diff_fit_ave = mean(frequency_diff_fit);
    Q_factor_ave = mean(Q_factor);

    stackedFig = figure('visible','off');
    stackedFig.Position = [50 50 650 900];
    plot_frequency_diff_fit = [frequency_diff_fit; NaN];
    s = stackedplot(resonance_wavelength_fit, ...,
        [plot_frequency_diff_fit, Q_factor/1000, Extinction_ratio]...,
    );
    % setting stacked fig properties
    s.FontSize = 16;
    s.FontName = 'Arial';
    s.XLabel = 'Wavelength [nm]';
    s.LineStyle = 'none';
    s.Color = 'k';
    s.LineWidth = 1;
    s.Marker = 'o';
    s.MarkerSize = 6;
    s.MarkerFaceColor = 'k';
    s.MarkerEdgeColor = 'none';
    s.DisplayLabels = {'FSR [GHz]', 'Q [\times10^3]', 'Extinction ratio'};
    sax = findobj(s.NodeChildren, 'Type','Axes');
    set([sax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
    sax(1).YLabel.Interpreter = 'tex';
    sax(2).YLabel.Interpreter = 'tex';
    sax(3).YLabel.Interpreter = 'tex';
    s.Title = typeName(RingTypeIndex) + " gap = " + GapStr + " " + DropThroughProp;
    %% save output
    % image
    saveas(normalizedTransmissionFig, curOutputPrefix + "normTrans", 'fig');
    saveas(normalizedTransmissionFig, curOutputPrefix + "normTrans", 'png');
    saveas(normalizedTransmissionFig, curOutputPrefix + "normTrans", 'eps');
    saveas(fitFig, logPrefix + "fit", 'fig');
    saveas(fitFig, logPrefix + "fit", 'png');
    saveas(fitFig, logPrefix + "fit", 'eps');
    saveas(stackedFig, curOutputPrefix + "properties", 'fig');
    saveas(stackedFig, curOutputPrefix + "properties", 'png');
    saveas(stackedFig, curOutputPrefix + "properties", 'eps');
    % data
    save(curOutputPrefix +     "FSR.mat", 'resonance_wavelength_fit', 'plot_frequency_diff_fit');
    save(curOutputPrefix + "Qfactor.mat", 'resonance_wavelength_fit', 'Q_factor');
    save(curOutputPrefix +   "ratio.mat", 'resonance_wavelength_fit', 'Extinction_ratio');
    FSR_dat=    [resonance_wavelength_fit, plot_frequency_diff_fit];
    Qfactor_dat=[resonance_wavelength_fit, Q_factor               ];
    ratio_dat=  [resonance_wavelength_fit, Extinction_ratio       ];
    save(curOutputPrefix +     "FSR.dat",     'FSR_dat', '-ascii', '-double', '-tabs');
    save(curOutputPrefix + "Qfactor.dat", 'Qfactor_dat', '-ascii', '-double', '-tabs');
    save(curOutputPrefix +   "ratio.dat",   'ratio_dat', '-ascii', '-double', '-tabs');
    % delete
    delete(normalizedTransmissionFig);
    delete(fitFig);
    delete(stackedFig);
end
