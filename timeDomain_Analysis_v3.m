
clc
clc;
clear;
close all ;
addpath 'subfunction\'
%% parameter config
% filename = '.\data\Hex_otic_data_start0_8prb_mcs0.dat';% prb 0-8 symbol 1-13
% slotTable = [9 10 11 12 13 14 15 16 17 18 19 0 1 2 3 4 5];%slotId6开始

filename = '.\data\log2\dg\Hex_otic_data.f000.dat';% prb 0-8 symbol 1-13
slotTable = [16 18 0 2 4 6 8 10 12 14 16 18 0]; %use air-slot
% filename = '.\data\log2\dg\Hex_otic_data.f667.dat';% prb 0-8 symbol 1-13
% slotTable = [0 2 4 6 8 10 12 14 16 18 0]; %use air-slot

RB_S = 0;
NumPRB = 273;
realSymSeq = [  1 2   4 5 6 7 8 9 10 11 12 13]; %0-13
Modulation = '256QAM';
% slotTable = [9 10 11 12 13 14 15 16 17 18 19 0 1 2 3 4 5];

%%  parameter fix
NSCID = 0;
NIDNSCID = 60;
Nfft = 4096;
AntNum = 1;
numLayer = 1;
dmrsPos = 3;%0-13
STOFlag = 1;
delay = 0;%800;% constell worse, need check
%%  parameter calc
% subCarrierNumOneSymbol = NumPRB*12;%0-3275
symbolNum = size(realSymSeq,2);
symbolS = realSymSeq(1);%0-13
symbolE = realSymSeq(symbolNum);%0-13
%% carrier pusch
[carrier, pusch] = genCarrierAndPusch(symbolS, symbolE, RB_S, NumPRB, Modulation, dmrsPos, NSCID, NIDNSCID);
%% read otic dat file
fid = fopen(filename,'r');
C = textscan(fid,'%s %s');
dataIQ = C{1,2};
[timePointNum, col] = size(dataIQ);
% numPerSlot = 61440; %352+288*13+4096*14;
timeI = zeros(1,timePointNum);
timeQ = zeros(1,timePointNum);
timedata = zeros(1,timePointNum);
for Idx = 1:1:timePointNum
    tmp = dataIQ{Idx,1};
    tmp1 = erase(tmp,'0x');
    timeQ(1, Idx) = hex2dec(tmp1(1:4));
    timeI(1, Idx) = hex2dec(tmp1(5:8));   
    timedata(1, Idx) = fix2float([timeQ(1, Idx) timeI(1, Idx)], 15, 1 );
end
%画时域
% figure;plot(real(timedata(1:2*61440)));grid on; 
figure;plot(abs(timedata(1*61440+1:6*61440)));grid on;
DL = 1;
timedata1 = zeros(1,timePointNum);
[phase_STO, a_angel] = STO_phase_compensation(DL);

%% time to freq
% for slotId = 1: 1 %slotNum
    slotId = 1;
    freqdata3D = zeros(4096,14,AntNum);
    for antIdx = 1:AntNum
        figure;   
        for symIdx = 1:14 
            if symIdx == 1
                cplen = 352 + delay;              
            else
                cplen = cplen + 288;                
            end
            symStart = (slotId-1)*61440 + cplen + (symIdx - 1)*4096  + 1;
            symEnd = (slotId-1)*61440 + cplen + symIdx*4096;             
            x = symStart:symEnd;
            % STO phase compensation
            if  STOFlag == 1 
                timedata1(x) = timedata(x) * exp(-1j*a_angel(symIdx));                  
            else
                 timedata1(x) = timedata(x);
            end
            % 画时域
%             figure;plot(abs(timedata(x)));grid on; % symbol timedomain
%             figure;plot(abs(timedata(1:61440)));grid on; %slot
                    
            fft_data =  fft(timedata1(x));           
            fft_result_shift = fftshift(fft_data);            
            freqdata3D(:,symIdx,1) = fft_result_shift.'; % freqdata after fftshift
            % 画频域  
            y = (-4096/2+1:1:4096/2)/4096*100;
            subplot(3,5,symIdx);plot(y, 10*log10(abs(fft_result_shift).^2));grid on;
            title(['Spectrum - sym (' num2str(symIdx) ') ant (' num2str(AntNum) ')']);
            xlabel("Hz");  ylabel("Spectrum"); 
        end
    end
    %% plot rxgrid
    rxgrid = freqdata3D(410+1:4096-410,:,:);
    pow_rxgrid =abs(rxgrid).^2;
    Sre = (1/Nfft.^2) * pow_rxgrid;
    for antIdx = 1:AntNum
        sre_db_1 = pow2db(Sre(:, :, antIdx));
        figure;
        imagesc([0 13], [0 273 * 12 - 1], sre_db_1 );
        colorbar
        axis xy;
        title(['Resource Grid (Antenna ' num2str(antIdx) ') - PUSCH']);
        xlabel("OFDM Symbol");
        ylabel("Subcarrier");    
    end

%%  channel est   % rxgrid[scidx,symbolIdx,rxAntIdx];   
    carrier.NSlot = slotTable(slotId);
    dmrsLayerIndices = nrPUSCHDMRSIndices(carrier,pusch); 
    dmrsLayerSymbolsLocal = nrPUSCHDMRS(carrier, pusch);
    [estChannelGrid,noiseEst] = nrChannelEstimate(rxgrid,dmrsLayerIndices,dmrsLayerSymbolsLocal);
    
%%  equalization
    [puschIndices,info,ptrsInd] = nrPUSCHIndices(carrier,pusch);
    [puschRx,puschHest] = nrExtractResources(puschIndices,rxgrid,estChannelGrid);
    [puschEq, csi] = nrEqualizeMMSE(puschRx,puschHest,noiseEst);
    % normalization
    if numLayer == 2
        sum_power_puschEq = sum(sum(abs(puschEq).^2));
    else
        sum_power_puschEq = sum(abs(puschEq).^2);
    end
    avg_power_puschEq = sum_power_puschEq /(length(puschEq(:,1)) * numLayer);
    for idx_layer = 1:numLayer
        puschEq(:,idx_layer) = puschEq(:,idx_layer)/sqrt(avg_power_puschEq);
    end
    % plot sym constellation and EVM   
    plotConstellationAndEVM(puschEq, carrier, pusch, realSymSeq, symbolNum, puschIndices, rxgrid);

    %% Demod and Descramble 
    pusch.TransmissionScheme  = 'nonCodebook' ;
    [ulschLLRs, rxSymbol] = nrPUSCHDecode(carrier,pusch, puschEq, noiseEst);
    
    csi = nrLayerDemap(csi);
    Qm = size(ulschLLRs,1) / size(rxSymbol,1);
    csi = reshape(repmat(csi{1}.', Qm, 1), [], 1);
    ulschLLRs = ulschLLRs .* csi;

    %% CRC Decode
    [decbits, blkerr] = decodefunction(ulschLLRs, pusch);


    %end for slot
