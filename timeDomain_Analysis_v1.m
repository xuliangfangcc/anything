
clc
clear
close all 


%% parameter config
NumPRB = 8;
symbolS = 1;%0-13
symbolE = 13;%0-13
symbolNum = 12;
subCarrierS = 0;%0-3275
subCarrierNumOneSymbol = NumPRB*12;%0-3275
modulation = 'QPSK';
%%  parameter fix
Nfft = 4096;
AntNum = 1;
numLayer = 1;
dmrsPos = 3;%0-13
STOFlag = 1;
delay = 0;%800;% constell worse, need check
%% carrier pusch
carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = 30;
carrier.NSizeGrid = NumPRB;

pusch = nrPUSCHConfig;
pusch.PRBSet = subCarrierS/ 12 : NumPRB - 1;
pusch.SymbolAllocation = [symbolS symbolNum];
pusch.NSizeBWP = NumPRB;
pusch.NStartBWP = floor(subCarrierS / 12);
pusch.Modulation = modulation;

%% read otic dat file
% filename = '.\data\Hex_otic_data_pds_pdc.dat';
% filename = '.\data\Hex_otic_data_pds.dat';% prb 0-273 symbol 1-12
filename = '.\data\Hex_otic_data_start0_8prb_mcs0.dat';% prb 0-8 symbol 1-13
%  filename = '.\data\Hex_otic_data_start0_8prb_mcs27.dat';% prb 0-8 symbol 1-13

fid = fopen(filename,'r');
C = textscan(fid,'%s %s');
dataIQ = C{1,2};
[timePointNum, col] = size(dataIQ);

numPerSlot = 61440; %352+288*13+4096*14;
slotNum = floor(timePointNum / numPerSlot);

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
figure;plot(abs(timedata(1*61440:6*61440)));grid on;
DL = 1;
timedata1 = zeros(1,timePointNum);
[phase_STO, a_angel] = STO_phase_compensation(DL);

%% time to freq
% for slotId = 1: 1 %slotNum
    slotId = 6;
    freqdata3D = zeros(4096,14,AntNum);
    for antIdx = 1:1%AntNum
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
    dmrsLayerIndices = (dmrsPos*3276 + subCarrierS + 1):2:(dmrsPos*3276 + subCarrierNumOneSymbol - 1);
    dmrsLayerIndices = dmrsLayerIndices.';
    dmrsLayerSymbols = rxgrid(1:2:subCarrierNumOneSymbol,dmrsPos + 1,1);   
   
   [estChannelGrid,noiseEst] = nrChannelEstimate(rxgrid,dmrsLayerIndices,dmrsLayerSymbols);
  
%% phase-compensation 
    y1 = interp1(1:NumPRB*12/2, dmrsLayerSymbols.', 1:0.5:(NumPRB*12/2), 'nearest');%'spline', 'nearest','linear','cubic'
    y1(NumPRB*12) = y1(NumPRB*12-1);
    phase = angle(y1(1:NumPRB*12));
%     figure; plot(phase,'*');grid on
%     figure;plot(unwrap(angle(dmrsLayerSymbols(1:48,1))),'*');
%     figure;plot(unwrap(angle(y1(1:96))),'*');
    rxgrid_modify = zeros(3276,14,1);
    symbolReq = [ 1 2  4 5 6 7 8 9 10 11 12 13];
    for reqIdx = 1:12
        symIdx = symbolReq(reqIdx) + 1;
        for scIdx = 1:NumPRB*12  
            rxgrid_modify(scIdx,symIdx,1) = rxgrid(scIdx,symIdx,1)*exp(-1j*(phase(scIdx)+pi/4)); 
        end
    end
  
%     figure;plot(unwrap(angle(rxgrid(1:96,2,1))),'*');
%     figure;plot(unwrap(angle(rxgrid_modify(1:96,2,1))),'*');    
%     figure;plot(real(rxgrid_modify(1:96,2,1)), imag(rxgrid_modify(1:96,2,1)),'*');    

%%  equalization
    for symIdx = symbolS+1:symbolE+1
        if symIdx == (symbolS + 1)
            puschIndices_S = symbolS * 3276 + subCarrierS + 1;
            puschIndices_E = puschIndices_S + subCarrierNumOneSymbol - 1;
            puschIndices = (puschIndices_S:puschIndices_E);
        elseif symIdx == (dmrsPos + 1)
            puschIndices_S = puschIndices_S + 3276;
            puschIndices_E = puschIndices_E + 3276;
        else            
            puschIndices_S = puschIndices_S + 3276;
            puschIndices_E = puschIndices_E + 3276;
            puschIndices = [puschIndices  puschIndices_S:puschIndices_E];
        end
    end
    puschIndices = puschIndices.';    
    [puschRx,puschHest] = nrExtractResources(puschIndices,rxgrid_modify,estChannelGrid);
    
    [puschEq, csi] = nrEqualizeMMSE(puschRx,puschHest,noiseEst);
    % normalization
    if numLayer == 2
        sum_power_puschEq = sum(sum(puschEq).^2);
    else
        sum_power_puschEq = sum(puschEq.^2);
    end
    avg_power_puschEq = sum_power_puschEq /(length(puschEq(:,1)) * numLayer);
    for idx_layer = 1:numLayer
        puschEq(:,idx_layer) = puschEq(:,idx_layer)/sqrt(avg_power_puschEq);
    end
    % plot sym constellation
    [dataSubCarri, numLayer]  = size(puschEq);
    symSeq = 1: dataSubCarri/subCarrierNumOneSymbol;
    realSymSeq = [1 2 4 5 6 7 8 9 10 11 12 13]; %0-13
    
    for idx_layer = 1: numLayer
        figure;
        for symIdx = 1: dataSubCarri/subCarrierNumOneSymbol                
                subplot(3,5,realSymSeq(symIdx)+1 )
                x = subCarrierNumOneSymbol*(symIdx-1)+1:subCarrierNumOneSymbol*symIdx;
                plot(real(puschEq(x,idx_layer)), imag(puschEq(x,idx_layer)),'.');grid on;
                title(['constellation sym-' num2str(realSymSeq(symIdx))]);            
        end
    end

%end for slot
