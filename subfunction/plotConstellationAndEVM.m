function plotConstellationAndEVM(puschEq, carrier, pusch, realSymSeq, symbolNum, puschIndices, rxgrid)
    %% plot sym constellation
    subCarrierNumOneSymbol = size(pusch.PRBSet, 2) * 12;
    for idx_layer = 1: pusch.NumLayers
        figure;
        for symIdx = 1:symbolNum               
                subplot(3,5,realSymSeq(symIdx)+1 )
                x =  subCarrierNumOneSymbol *(symIdx-1)+1:subCarrierNumOneSymbol*symIdx;
                plot(real(puschEq(x,idx_layer)), imag(puschEq(x,idx_layer)),'.');grid on;
                title(['constellation sym-' num2str(realSymSeq(symIdx)) ' slot' num2str(carrier.NSlot)]);            
        end
    end
    %% plot sym EVM

    [refConstellation, modBits, evmBase] = getConstellationRefPoint(pusch.Modulation);
    err_abs = zeros(length(puschEq(:,1)), 2^modBits);
    for layerId = 1 : pusch.NumLayers
        puschEqPerLayer = puschEq(:, layerId);
        for pointIdx = 1 : 2^modBits
            err_abs(:,pointIdx) = abs(puschEqPerLayer - refConstellation(pointIdx) );
        end
        for scIdx = 1 : length(puschEqPerLayer)
            [minvalue, minIndex(scIdx, layerId)] = min(err_abs(scIdx,:));
        end
    end
    puschRefSymbol = zeros(size(puschEqPerLayer,1), pusch.NumLayers);
    for layerId = 1 : pusch.NumLayers
        for scIdx = 1 : length(puschEqPerLayer)
            puschRefSymbol(scIdx, layerId) = refConstellation(minIndex(scIdx, layerId));
        end
    end

%
persistent slotEVM;
persistent rbEVM;
persistent evmPerSlot;

slotEVM = comm.EVM;
rbEVM = comm.EVM;
% evmPerSlot = NaN(carrier.NSlot, pusch.NumLayers);
evmPerSlot = slotEVM(puschRefSymbol, puschEq);
figure;
subplot(2,1,1);
plot(evmPerSlot, 'o-') ; hold on;
plot(evmBase, 'r-.');grid on;
xlabel(['slot' num2str(carrier.NSlot)]);ylabel('EVM(%)');
% legend('layer' + (1:pusch.NumLayers), 'Location','eastoutside');
title('EVM per layer');


siz = size(rxgrid);% 3276*14*2
[k,~,p] = ind2sub(siz,puschIndices);
subs = k; %allsym
scNumOneSymbol = siz(1);% sc one sym
evmPerSc = NaN(scNumOneSymbol, pusch.NumLayers);

for layerIdx = 1:pusch.NumLayers
    for sc = unique(subs).'
        this = (subs == sc & p == layerIdx);
        evmPerSc(sc,layerIdx) = rbEVM(puschRefSymbol(this), puschEq(this));
    end
end
subplot(2,1,2);
plot(0:(scNumOneSymbol -1 ), evmPerSc(1:end, 1),'x-');grid on; hold on;
if size(evmPerSc, 2) >  1
    plot(0:(scNumOneSymbol - 1), evmPerSc(1:end, 2), 'x-'); grid on; hold on;
end
plot([0 (scNumOneSymbol - 1)], [evmBase, evmBase], 'r-.');
xlabel('subcarriers');ylabel('evm(%)');hold off;

%% plot evm per symbol
evmPerSymbol = NaN(scNumOneSymbol, symbolNum, pusch.NumLayers);
figure;
for layerIdx = 1:pusch.NumLayers
    for symIdx = 1:symbolNum    
        scIndex = scNumOneSymbol * (symIdx - 1) + 1 : scNumOneSymbol * symIdx;
        puschRefSymbolPerSym = puschRefSymbol(scIndex);
        puschEqPerSym = puschEq(scIndex);
        for sc = 1:scNumOneSymbol
            evmPerSymbol(sc,symIdx,layerIdx) = rbEVM(puschRefSymbolPerSym(sc), puschEqPerSym(sc));
        end
        subplot(3,5,realSymSeq(symIdx)+1 );
        plot(evmPerSymbol(:,symIdx,layerIdx), '*-');grid on; hold on;
        plot([0 (scNumOneSymbol - 1)], [evmBase, evmBase], 'r-.');
        xlabel('subcarriers');ylabel('evm(%)');hold off;
        title(['EVM sym-' num2str(realSymSeq(symIdx)) ' slot' num2str(carrier.NSlot)]);            
    end
end



end