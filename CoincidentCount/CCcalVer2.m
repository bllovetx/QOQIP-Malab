% This version(2) considers out of bound case
% to handle this 

% load data
%allData = importdata("data\Timetags_CH35_13dBm_DBT400_CH37-CH33_CH38-CH32_100s.dat");
allData = importdata("data\Timetags_10.0mW_30s.txt");
channel1 = [];
channel2 = [];
channel3 = [];
channel4 = [];
for i = 1:length(allData)
    switch allData(i, 1)
    case 1
        channel1(end+1) = allData(i, 2);
    case 2
        channel2(end+1) = allData(i, 2);
    case 3
        channel3(end+1) = allData(i, 2);
    case 4
        channel4(end+1) = allData(i, 2);
    end
end

delayStep = 2.5 %ns
windowWidth = 2.5 %ns
timeInterval = 30%s

delays = -20:delayStep:20;
figure
plot(delays, CoincidentCountCal(channel1, channel2, timeInterval, windowWidth, delays));
hold on
plot(delays, CoincidentCountCal(channel3, channel4, timeInterval, windowWidth, delays));







function coinCounts = CoincidentCountCal(channel_one, channel_two, T, window_ns, delays)
    %   Calculate Coincident Count of Two Channel with Given Windows
    %   and Delays
    %   Para:
    %       channel_(one/two): time stamps of two channels
    %       T: total time interval(s)
    %       window_ns: coincident window in ns
    %       delays: list of delay in ns
    %   Return:
    %       coinCounts: list of cc (for each delay)

    %   parameters:
    timeResol = 0.15625; %(ns)

    %   init
    delayNum = length(delays);
    ch1Num = length(channel_one);
    ch2Num = length(channel_two);
    ccwind = window_ns / timeResol; % time stamp
    coinCounts = zeros(1, delayNum);
    startInds = zeros(1, delayNum) + 1;
    endInds = zeros(1, delayNum) + 1;
    delayStamps = delays / timeResol;   % time stamp
    startDelays = delayStamps - ccwind/2;%time stamp
    endDelays = delayStamps + ccwind/2; % time stamp

    %   Moving windows
    TStamp = T*1000000000 / timeResol;
    ch1Start = max(0, -min(delayStamps));   % time stamp
    ch1End   = TStamp + min(0, -max(delayStamps));  % time Stamp
    i = 1;
    while channel_one(i) < ch1Start
        i = i + 1;
    end
    while i <= ch1Num && channel_one(i) < ch1End
        curStarts = channel_one(i) + startDelays;
        curEnds   = channel_one(i) + endDelays;
        for j = 1:delayNum
            % update index of first ch2-stamp after each Win_start/end
            while startInds(j) <= ch2Num && channel_two(startInds(j)) < curStarts(j)
                startInds(j) = startInds(j) + 1;
            end
            while endInds(j) <= ch2Num && channel_two(endInds(j)) < curEnds(j)
                endInds(j) = endInds(j) + 1;
            end
            % calculate ch2 num in each window
            curNum = endInds(j) - startInds(j);
            if curNum
                if curNum ~= 1
                    startInds(j)
                    endInds(j)
                    sprintf('%12f', channel_two(startInds(j)))
                    sprintf('%12f', channel_two(startInds(j)+1))
                    sprintf('%12f', channel_two(startInds(j)+2))
                end
                assert(curNum == 1, "More than one photon in one window? Should not happen!");
                coinCounts(j) = coinCounts(j) + 1;
            end
        end
        i = i + 1;
    end
end




