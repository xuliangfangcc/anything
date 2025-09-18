function [refPoints, modBits, evmBase] = getConstellationRefPoint(mod)
    switch mod
        case 'BPSK'
            nPts = 2;
            modBits = 1;
        case 'QPSK'
            nPts = 4;
            modBits = 2;
            evmBase = 17.5;
        case '16QAM'
            nPts = 16;
            modBits = 4;
            evmBase = 12.5;
        case '64QAM'
            nPts = 64; 
            modBits = 6;
            evmBase = 8;
        case '256QAM'
            nPts = 256;
            modBits = 8;
            evmBase = 3.5;
    end
    
    binaryValuesMat = de2bi(0 : nPts - 1, log2(nPts), 'left-msb');
    binaryValues_trans = binaryValuesMat' ;
    binaryValues = binaryValues_trans(:);
    refPoints = nrSymbolModulate(binaryValues, mod);
end