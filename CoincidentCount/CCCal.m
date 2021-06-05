% load data
allData = importdata("data\Timetags_CH35_13dBm_DBT400_CH37-CH33_CH38-CH32_100s.dat");
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

delays = -1000:0.1:1000;
figure
plot(delays, CoincidentCountCal(channel1, channel2, 1, delays));
hold on
plot(delays, CoincidentCountCal(channel3, channel4, 1, delays));







function coinCounts = CoincidentCountCal(channel_one, channel_two, window_ns, delays)
    %   Calculate Coincident Count of Two Channel with Given Windows
    %   and Delays
    %   Para:
    %       channel_(one/two): time stamps of two channels
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
    ccwind = window_ns / timeResol;
    coinCounts = zeros(1, delayNum);
    startInds = zeros(1, delayNum) + 1;
    endInds = zeros(1, delayNum) + 1;
    delayStamps = delays / timeResol;
    startDelays = delayStamps - ccwind/2;
    endDelays = delayStamps + ccwind/2;

    %   Moving windows
    %   TODO: is it necessary to remove data at two ends?
    for i = 1:ch1Num
        curStarts = channel_one(i) + startDelays;
        curEnds   = channel_one(i) + endDelays;
        for j = 1:delayNum
            % update index of first ch2-stamp after each Win_start/end
            while channel_two(startInds(j)) < curStarts(j) && startInds(j) < ch2Num
                startInds(j) = startInds(j) + 1;
            end
            while channel_two(endInds(j)) < curEnds(j) && endInds(j) < ch2Num
                endInds(j) = endInds(j) + 1;
            end
            % calculate ch2 num in each window
            curNum = endInds(j) - startInds(j);
            if curNum
                assert(curNum == 1, "More than one photon in one window? Should not happen!");
                coinCounts(j) = coinCounts(j) + 1;
            end
        end
    end
end




