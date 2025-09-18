function [carrier, pusch] = genCarrierAndPusch(symbolS, symbolE, RB_S, NumPRB, Modulation, dmrsPos, NSCID, NIDNSCID)
    carrier = nrCarrierConfig;

    carrier.SubcarrierSpacing = 30;
    carrier.NSizeGrid = 273;
    carrier.NStartGrid = 0;
    carrier.NCellID = 0;
    carrier.NSlot = 0;
    carrier.NFrame = 1;

    pusch = nrPUSCHConfig;

    pusch.PRBSet = RB_S : RB_S + NumPRB - 1;
    pusch.SymbolAllocation = [symbolS symbolE];
    pusch.NSizeBWP = 273;
    pusch.NStartBWP = RB_S;
    pusch.Modulation = Modulation;%'QPSK';
    pusch.MappingType = 'A';
    pusch.NumLayers = 1;   

    pusch.DMRS.DMRSConfigurationType = 1; % 1-->type0; 2:type1
    pusch.DMRS.DMRSTypeAPosition = dmrsPos;
    pusch.DMRS.DMRSAdditionalPosition = 0;
    pusch.DMRS.DMRSLength = 1;
    pusch.DMRS.NSCID = NSCID;
    pusch.DMRS.NIDNSCID = NIDNSCID;
    pusch.DMRS.NumCDMGroupsWithoutData = 2;

%     pusch.DMRS.DMRSPortSet = [0 3]; % 
    if pusch.DMRS.DMRSLength == 1 
        if pusch.DMRS.DMRSConfigurationType == 1  % tmp this way
            pusch.DMRS.DMRSPortSet = [0 3]; % 2 port [0 3] ; 1 port [0]
        elseif pusch.DMRS.DMRSConfigurationType == 2
            pusch.DMRS.DMRSPortSet = [0 5];  
        end
    elseif pusch.DMRS.DMRSLength == 2
        if pusch.DMRS.DMRSConfigurationType == 1  % tmp this way
            pusch.DMRS.DMRSPortSet = [0 7];
        elseif pusch.DMRS.DMRSConfigurationType == 2
            pusch.DMRS.DMRSPortSet = [0 11];  
        end
    end
    if pusch.NumLayers == 1
        pusch.DMRS.DMRSPortSet = 0;
    end
end