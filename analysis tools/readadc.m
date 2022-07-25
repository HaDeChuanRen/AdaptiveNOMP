function data_cube = readadc(filename, options) 
    if nargin > 1
        numRX = options.numRX;
        numADCBits = options.numADCBits;
        numADCSamples = options.numADCSamples;
        isReal = options.isReal;             % set to 1 if real only data, 0 if complex data0    
    else
        numADCSamples = 256;    % number of ADC samples per chirp
        numADCBits = 16;        % number of ADC bits per sample
        numRX = 4;              % number of receivers
        isReal = 0;
    end    
    

    % read file
    fid = fopen(filename,'r');
    adcData = fread(fid, 'int16');
    % if 12 or 14 bits ADC per sample compensate for sign extension
    if numADCBits ~= 16
        l_max = 2^(numADCBits-1)-1;
        adcData(adcData > l_max) = adcData(adcData > l_max) - 2^numADCBits;
    end
    fclose(fid);
    fileSize = size(adcData, 1);
    % real data reshape, filesize = numADCSamples*numChirps
    if isReal
        numChirps = floor(fileSize/numADCSamples/numRX);
        fileSize_cut = numChirps * numADCSamples * numRX;
        adcData = adcData(1 : fileSize_cut);
        %create column for each chirp
        LVDS = reshape(adcData, numADCSamples*numRX, numChirps);
        %each row is data from one chirp
        LVDS = LVDS.';
    else
        % for complex data
        % filesize = 2 * numADCSamples*numChirps
        numChirps = floor(fileSize/2/numADCSamples/numRX);
        fileSize = numChirps * 2 * numADCSamples * numRX;
        adcData = adcData(1 : fileSize);
        LVDS = zeros(1, fileSize/2);
        % combine real and imaginary part into complex data
        % read in file: 2I is followed by 2Q
        counter = 1;
        for i = 1 : 4 : fileSize - 1
            LVDS(1,counter) = adcData(i) + sqrt(-1)*adcData(i+2); 
            LVDS(1,counter+1) = adcData(i+1)+sqrt(-1)*adcData(i+3);  %IQ数据合并成复数
            counter = counter + 2;
        end
        % create column for each chirp
        LVDS = reshape(LVDS, numADCSamples*numRX, numChirps);
        % each row is data from one chirp
        LVDS = LVDS.';
    end
    % organize data per RX
    adcData = zeros(numRX,numChirps*numADCSamples);
    for row = 1:numRX
        for i = 1: numChirps % 1：1024
            adcData(row, (i-1)*numADCSamples+1:i*numADCSamples) = LVDS(i, (row - 1)*numADCSamples+1:row*numADCSamples);
        end
    end
    data_cube = nan(numADCSamples, numChirps, numRX);
    % 4个RX通道数据，numChirps列，第i列为Chirp i的回波数据
    RX1data = reshape(adcData(1,:),numADCSamples,numChirps);   % RX1数据
    RX2data = reshape(adcData(2,:),numADCSamples,numChirps);   % RX2
    RX3data = reshape(adcData(3,:),numADCSamples,numChirps);   % RX3
    RX4data = reshape(adcData(4,:),numADCSamples,numChirps);   % RX4
    data_cube(:,:,1) = RX1data;
    data_cube(:,:,2) = RX2data;
    data_cube(:,:,3) = RX3data;
    data_cube(:,:,4) = RX4data;
end