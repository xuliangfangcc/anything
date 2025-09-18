function   [decbits, blkerr] = decodefunction(ulschLLRs, pusch)
%% config para from codeword DL TTI REQ
    TBS = 33822;
    TargetCodeRate = 9480 / 2^14;

    decodeULSCH = nrULSCHDecoder;

    decodeULSCH.TransportBlockLength = TBS;
    decodeULSCH.MultipleHARQProcesses = false;
    decodeULSCH.LimitedBufferRateRecovery = false;    % PrbSymbRateMatchPattern not support now
    decodeULSCH.TargetCodeRate = TargetCodeRate;
    decodeULSCH.LDPCDecodingAlgorithm = 'Normalized min-sum';
    decodeULSCH.MaximumLDPCIterationCount = 16 ;

    RedundancyVersion = 0; %     harqEntity.RedundancyVersion = 0;
%     HARQProcessID = 0; %     harqEntity.HARQProcessID = 0;
    [decbits, blkerr] = decodeULSCH(ulschLLRs, pusch.Modulation, pusch.NumLayers, RedundancyVersion);
    
    if blkerr == 1
        disp(' crc false');
    else
        disp('crc true');
    end



end