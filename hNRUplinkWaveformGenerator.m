%hNRUplinkWaveformGenerator Generate 5G NR uplink waveform
%   [WAVEFORM,GRIDSET,INFO] = hNRUplinkWaveformGenerator(WAVECONFIG) 
%   generates a 5G NR uplink waveform WAVEFORM with the input WAVECONFIG
%   parameters.
%
%   WAVECONFIG is a structure of structures with the following fields:
%   NCellID          - Physical layer cell identity. It is in range 0 to
%                      1007
%   ChannelBandwidth - Channel bandwidth in MHz
%   FrequencyRange   - Frequency range ('FR1','FR2'). The frequency range
%                      and the channel bandwidth must be a valid ones as
%                      per TS 38.101-1 for FR1 and 38.101-2 for FR2.
%   NumSubframes     - Number of 1ms subframes in generated waveform
%   DisplayGrids     - Display the grids after signal generation. It is
%                      either 0 or 1
%   Carriers         - Carrier(s) configuration. A structure array with the
%                      following fields:
%         SubcarrierSpacing - Subcarrier spacing configuration in kHz.
%                             Possible configurations 15,30,60,120 and 240
%                             for normal CP, and 60 for extended CP
%         NRB               - Number of resource blocks
%         RBStart           - Start of resource block for the carrier
%   BWP              - Bandwidth part configuration. A structure array with
%                      the following fields:
%         SubcarrierSpacing - Subcarrier spacing configuration in kHz.
%                             Possible configurations 15,30,60,120 and 240
%                             for normal CP, and 60 for extended CP. Note
%                             for each BWP with separate subcarrier
%                             spacing, there has to be a carrier with the
%                             same subcarrier spacing configuration
%         CyclicPrefix      - Cyclic prefix ('normal','extended')
%         NRB               - Number of resource blocks per bandwidth part
%         RBOffset          - Position of BWP in the carrier
%   PUCCH            - Physical uplink control channel configuration. A
%                      structure array with the following fields:
%         Enable               - Enable/Disable the PUCCH instance. It is
%                                either 0 or 1
%         BWP                  - Bandwidth part identity. It is indicated
%                                to which bandwidth part the PUCCH is
%                                configured
%         Power                - Power scaling (in dB)
%         AllocatedSlots       - Allocated slots with the period. A row
%                                vector indicating the slots within the
%                                repetition period
%         AllocatedPeriod      - Allocated slot period. Use empty ([]) for
%                                no repetition
%         RNTI                 - Radio network temporary identifier. It is
%                                in range 0 to 65535
%         NID                  - Scrambling identity. It is in range
%                                0 to 1023. Use empty ([]) to use
%                                physical layer cell identity NCellID
%         HoppingId            - Hopping identity. It is in range 0 to
%                                1023. It is used for formats 0/1/3/4. Use
%                                empty ([]) to use physical layer cell
%                                identity NCellID
%         NIDDMRS              - Scrambling identity. It is in range
%                                0 to 65535. Use empty ([]) to use for
%                                physical layer cell identity NCellID. It
%                                is used only for format 2
%         PowerDMRS            - Power scaling of DM-RS (in dB)
%         DedicatedResource    - Enable/Disable the dedicated resource. It
%                                is either 0 or 1. When the value is set to
%                                0, the PUCCH resource is based on the
%                                resource index as per TS 38.213, Section
%                                9.2.1
%         ResourceIndex        - Resource index for common resource. It is
%                                in range 0 to 15. The value is considered
%                                only when DedicatedResource is set to 0
%         StartPRB             - Index of first PRB prior to frequency
%                                hopping or for no frequency hopping
%         SecondHopPRB         - Index of first PRB after frequency hopping
%         IntraSlotFreqHopping - Intra-slot frequency hopping
%                                configuration. It is one of the set
%                                {'enabled',disabled'} 
%         GroupHopping         - Group hopping configuration. It is one of
%                                the set {'enable','disable','neither'}
%         PUCCHFormat          - PUCCH format number (0...4)
%         StartSymbol          - Starting symbol index in the slot
%         NrOfSymbols          - Number of OFDM symbols allocated for
%                                PUCCH. The sum of StartSymbol and
%                                NrOfSymbols should be less than the number
%                                of symbols in a slot
%         InitialCS            - Initial cyclic shift. It is applicable for
%                                formats 0 and 1. It is in range 0 to 11
%         OCCI                 - Orthogonal cover code index. It is
%                                applicable for formats 1 and 4. For format
%                                1, the range is in 0 to 6 and for format
%                                4, it must be less than the spreading
%                                factor
%         Modulation           - Modulation scheme. It is applicable for
%                                formats 3 and 4. It is either 'pi/2-BPSK'
%                                or 'QPSK'
%         NrOfRB               - Number of resource blocks. It is
%                                applicable for formats 2 and 3. It is one
%                                of the set {1,2,3,4,5,6,8,9,10,12,15,16}
%         SpreadingFactor      - Spreading factor. It is applicable for
%                                format 4. The value is either 2 or 4
%         AdditionalDMRS       - Additional DM-RS configuration. It is
%                                applicable for formats 3 and 4. It is
%                                either 0 or 1
%         NrOfSlots            - Number of slots for repetition. It is one
%                                of the set (1,2,4,8). Use 1 for no
%                                repetition. For any other value, the
%                                values in the AllocatedSlots are treated
%                                as the start of the slot for the
%                                repetition. These slots are then repeated
%                                over the period
%         InterSlotFreqHopping - Inter-slot frequency hopping
%                                configuration. It is one of the set
%                                {'enabled',disabled'}. When the value is
%                                set to 'enabled', intra-slot frequency
%                                hopping is considered as disabled
%         MaxCodeRate          - Maximum code rate. It is one of the set
%                                {0.08, 0.15, 0.25, 0.35, 0.45, 0.6, 0.8}.
%                                This parameter is used when there is
%                                multiplexing of UCI-part1 and UCI-part2
%         LenACK               - Number of hybrid automatic repeat request
%                                acknowledgment (HARQ-ACK) bits. For
%                                formats 0 and 1, the maximum value can be
%                                2. For other formats, there is no
%                                restriction. Set the value to 0, for no
%                                HARQ-ACK transmission
%         LenSR                - Number of scheduling request (SR) bits.
%                                For formats 0 and 1, the maximum value can
%                                be 1. Set the value to 0, for no SR
%                                transmission
%         LenCSI1              - Number of CSI part 1 bits. It is
%                                applicable for formats 2,3 and 4. Set the
%                                value to 0, for no CSI part 1
%                                transmission
%         LenCSI2              - Number of CSI part 2 bits. It is
%                                applicable for format 3 and 4. Set the
%                                value to 0, for no CSI part 2
%                                transmission. This value is ignored if
%                                there is no CSI part 1 bits
%       Note that the UCI multiplexing happens on PUCCH when LenCSI2 is not
%       zero for formats 3 and 4.
%         DataSource           - UCI data source. Use one of the following
%                                standard PN sequences: 'PN9-ITU', 'PN9',
%                                'PN11', 'PN15', 'PN23'. The seed for the
%                                generator can be specified using a cell
%                                array in the form {'PN9',seed}. If no seed
%                                is specified, the generator is initialized
%                                with all ones
%   PUSCH    - Physical uplink shared channel configuration. A structure
%              array with the following fields:
%         Enable               - Enable/Disable the PUSCH configuration. It
%                                is either 0 or 1
%         BWP                  - Bandwidth part identity. It is indicated
%                                to which bandwidth part PUSCH is
%                                configured
%         Power                - Power scaling (in dB)
%         TargetCodeRate       - Code rate used to calculate transport
%                                block sizes
%         Xoh_PUSCH            - Overhead parameter. It is one of the set
%                                {0,6,12,18}
%         TxScheme             - Transmission scheme. It is one of the set
%                                {'codebook','nonCodebook'}. When set to
%                                codebook, the precoding matrix is applied
%                                based on the number of layers and the
%                                antenna ports
%         Modulation           - Modulation scheme. It is one of the set
%                                {'pi/2-BPSK','QPSK','16QAM','64QAM','256QAM'}
%         NLayers              - Number of layers. It is in range 1 to 4
%         NAntennaPorts        - Number of antenna ports. It is one of the
%                                set {1,2,4}
%         RVSequence           - Redundancy version sequence. it is applied
%                                cyclically across the PUSCH allocation
%                                sequence
%         IntraSlotFreqHopping - Intra-slot frequency hopping
%                                configuration. It is one of the set
%                                {'enabled',disabled'}
%         TransformPrecoding   - Transform precoding flag. It is either 0
%                                or 1. When set to 1, DFT operation is
%                                applied before precoding and it is only
%                                used for single layer
%         TPMI                 - Transmitted precoding matrix indicator.
%                                The range depends on the number of layers
%                                and antenna ports
%         GroupHopping         - Group hopping configuration. It is one of
%                                the set {'enable','disable','neither'}
%         RBOffset             - Resource block offset for second hop
%         InterSlotFreqHopping - Inter-slot frequency hopping
%                                configuration. It is one of the set
%                                {'enabled',disabled'}. When the value is
%                                set to 'enabled', intra-slot frequency
%                                hopping is considered as disabled
%         NID                  - Scrambling identity. It is in range
%                                0 to 1023. Use empty ([]) to use physical
%                                layer cell identity
%         RNTI                 - Radio network temporary identifier
%         DataSource           - Transport block data source. Use one of
%                                the following standard PN sequences:
%                                'PN9-ITU', 'PN9', 'PN11', 'PN15', 'PN23'.
%                                The seed for the generator can be
%                                specified using a cell array in the form
%                                {'PN9',seed}. If no seed is specified, the
%                                generator is initialized with all ones
%         PUSCHMappingType     - PUSCH mapping type. It is either 'A' or
%                                'B'. The number of symbols in a slot and
%                                the start symbol depends on the mapping
%                                type
%         AllocatedSymbols     - Symbols in a slot. It needs to be a
%                                contiguous allocation. For PUSCH mapping
%                                type 'A', the start symbol must be zero
%                                and the length can be from 4 to 14 (for
%                                normal CP) and up to 12 (for extended CP)
%         AllocatedSlots       - Slots in a frame used for PUSCH
%         AllocatedPeriod      - Allocation period in slots. Use empty for
%                                no repetition
%         AllocatedPRB         - PRB allocation
%         PortSet              - DM-RS ports to use for the layers
%         DMRSTypeAPosition    - DM-RS symbol position for mapping type
%                                'A'. It is either 2 or 3
%         DMRSLength           - Number of front-loaded DM-RS symbols. It
%                                is either 1 (single symbol) or 2 (double
%                                symbol)
%         DMRSAdditionalPosition  - Additional DM-RS symbols positions. It
%                                   is in range 0 to 3. Value of zero
%                                   indicates no additional DM-RS symbols
%         DMRSConfigurationType   - DM-RS configuration type. It is either
%                                   1 or 2. The number of subcarriers
%                                   allocated for DM-RS depends on the type
%         NumCDMGroupsWithoutData - Number of DM-RS CDM groups without
%                                   data. It is value in range 1 to 3. This
%                                   helps in reserving the resource
%                                   elements in the DM-RS carrying symbols
%                                   to avoid data transmission
%         NIDNSCID             - Scrambling identity for CP-OFDM DMRS. It
%                                is in range 0 to 65535. Use empty ([]) to
%                                use physical layer cell identity NCellID
%         NSCID                - CP-OFDM DM-RS scrambling initialization.
%                                It is either 0 or 1
%         NRSID                - Scrambling identity for DFT-s-OFDM DM-RS.
%                                It is in range 0 to 1007. Use empty ([])
%                                to use the physical layer cell identity
%                                NCellID
%         PowerDMRS            - Power boosting for DM-RS (in dB)
%         DisableULSCH         - Disable UL-SCH on overlapping slots of
%                                PUSCH and PUCCH. It is either 0 or 1
%         BetaOffsetACK        - Rate matching offset for HARQ-ACK
%         BetaOffsetCSI1       - Rate matching offset for CSI part 1
%         BetaOffsetCSI2       - Rate matching offset for CSI part 2
%         ScalingFactor        - Scaling factor
%         EnablePTRS           - Enable the PT-RS configuration (0,1)
%         PTRSFrequencyDensity - Frequency density of PT-RS (2,4)
%         PTRSTimeDensity      - Time density of PT-RS (1,2,4)
%         PTRSNumSamples       - Number of PT-RS samples (2,4). It is used
%                                only for DFT-s-OFDM
%         PTRSNumGroups        - Number of PT-RS groups (2,4,8). It is used
%                                only for DFT-s-OFDM
%         PTRSREOffset         - Resource element offset ('00','01','10','11')
%         PTRSPortSet          - Antenna ports of PT-RS
%         PTRSNID              - Scrambling identity of PT-RS (0...1007).
%                                It is used only for DFT-s-OFDM
%         PowerPTRS            - Power boosting for PT-RS (dB)
%   SRS     - Sounding reference signal configuration. A structure array 
%             with the following fields:
%         Enable            - Enable/Disable the SRS configuration. It is
%                             either 0 or 1
%         BWP               - Bandwidth part identity. It is indicated
%                             to which bandwidth part SRS is configured
%         NumSRSPorts       - Number of SRS ports (1,2,4)
%         NumSRSSymbols     - Number of SRS symbols in a slot (1,2,4)
%         SymbolStart       - First OFDM symbol allocated to the SRS within
%                             a slot. It ranges (8...13) for normal cyclic
%                             prefix (CP) and (6...11) for extended CP
%         AllocatedSlots    - Allocated slots within the period. A row 
%                             vector indicating the slots within the
%                             repetition period. To configure periodic
%                             resource types and inter-slot
%                             frequency-hopping, AllocatedSlots must be a
%                             scalar lower than AllocatedPeriod. Otherwise,
%                             the SRS resource type is aperiodic (triggered
%                             by DCI)
%         AllocatedPeriod   - Allocated slot period. Use empty ([]) for
%                             no repetition
%         FreqStart         - Frequency position of the SRS in BWP in RBs
%         NRRC              - Additional offset from FreqStart specified in
%                             blocks of 4 PRBs (0..67)
%         CSRS              - Bandwidth configuration index C_SRS (0...63). 
%                             It controls the SRS allocated bandwidth and
%                             frequency hopping as defined in TS 38.211
%                             Table 6.4.1.4.3-1
%         BSRS              - Bandwidth configuration index B_SRS (0...3). 
%                             It controls the SRS allocated bandwidth and
%                             frequency hopping as defined in TS 38.211
%                             Table 6.4.1.4.3-1
%         BHop              - Frequency hopping (0..3). Set BHop < BSRS 
%                             to enable frequency hopping
%         KTC               - Comb number sets the subcarrier density (2,4) 
%         KBarTC            - Subcarrier offset (0..KTC-1)
%         CyclicShift       - Cyclic shift number (0...NCSmax-1). The 
%                             maximum number of cyclic shifts NCSmax = 12
%                             if KTC = 4 and NCSmax = 8 if KTC = 2
%         Repetition        - Repetition factor (1,2,4). Number of equal
%                             consecutive SRS symbols in a slot
%         GroupSeqHopping   - Group or sequence hopping 
%                             ('neither','groupHopping','sequenceHopping')
%         NSRSID            - SRS scrambling identity 
%
%   The parameters DisableULSCH, BetaOffsetACK, BetaOffsetCSI1,
%   BetaOffsetCSI2 and ScalingFactor are used to transmit UCI information
%   on PUSCH, whenever there is a overlapping in the slots of PUCCH and
%   PUSCH. The parameters can be set as per TS 38.212, Section
%   6.3.2.4.
%
%   For PUCCH resource with format as 1 or 3 or 4, repetition of slots and
%   inter-slot frequency hopping can be configured. In other cases, these
%   parameters are ignored.
%
%   Whenever there is a PUSCH and PUCCH in a same slot, only PUSCH is
%   transmitted with UCI information on PUSCH. In this case, the UCI
%   information is transmitted based on the PUSCH allocation.
%
%   For PUSCH multi-slot transmission, inter-slot frequency hopping can be
%   enabled, the starting PRB is the changed in each slot in the frame. The
%   slots with even indices start with RBStart and with odd indices starts
%   at sum of RBStart and RBOffset modulo NRB of associated BWP.
%
%   Copyright 2019-2020 The MathWorks, Inc.

function [waveform,gridset,winfo] = hNRUplinkWaveformGenerator(waveconfig)


    % Unbundle the channel specific parameter structures for easier access
    carriers = waveconfig.Carriers;
    bwp = waveconfig.BWP;
    pucch = waveconfig.PUCCH;
    pusch = waveconfig.PUSCH;
    srs = [];
    
    if isfield(waveconfig,'SRS')
        srs = waveconfig.SRS;
    end
    
    % Defaulting for the grid plotting
    if ~isfield(waveconfig,'DisplayGrids')
        waveconfig.DisplayGrids = 0;
    end

    % Define the instrumentation info variable 
    waveinfo = struct('PUSCH',[]); 
    
    % Check if NCellId is in valid range 0...1007
    if (waveconfig.NCellID < 0) || (waveconfig.NCellID > 1007)
        error('The NCellID must be in range 0 to 1007.');
    end

    % Cross-check the BWP and SCS carrier configurations
    carrierscs = [carriers.SubcarrierSpacing];
    for bp=1:length(waveconfig.BWP)
        % Map it into a SCS specific carrier level RE grid
        cidx = find(bwp(bp).SubcarrierSpacing == carrierscs,1);
        if isempty(cidx)
            error('A SCS specific carrier configuration for SCS = %d kHz has not been defined. This carrier definition is required for BWP %d.',bwp(bp).SubcarrierSpacing,bp);
        end
        % Record the carrier index associated with the BWP
        bwp(bp).CarrierIdx = cidx;
    end

    % Create BWP PRB resource grids
    ResourceGrids = arrayfun(@(bp)zeros(bp.NRB,waveconfig.NumSubframes*1*symbolsPerSlot(bp)*fix(bp.SubcarrierSpacing/15)),...
                                    bwp,'UniformOutput',false);

    % Create BWP subcarrier resource grids
    % Size ALL BWP RE grids by the number of layers/ports in the enabled PUSCH
    numBWP = length(bwp);
    maxlayers = ones(1,numBWP);
    for i =1:numBWP
        n = ones(1,length(pusch)+length(srs));
        for j = 1:length(pusch)
            if pusch(j).Enable
                if pusch(j).BWP == i
                   if (strcmpi(pusch(j).TxScheme,'codebook'))
                       numPorts = pusch(j).NAntennaPorts;
                       if isempty(numPorts) || numPorts <= 0
                           warning('For codebook based transmission, the number of antenna ports cannot be empty or non-positive. It is treated as 1 for PUSCH %d.',j);
                           pusch(j).NAntennaPorts = 1;
                           numPorts = 1;
                       end
                       n(j) = numPorts;
                   else
                       numLayers = pusch(j).NLayers;
                       if isempty(numLayers) || numLayers <= 0
                           warning('The number of layers cannot be empty or non-positive. It is treated as one for PUSCH %d.',j);
                           pusch(j).NLayers = 1;
                           numLayers = 1;
                       end
                       n(j) = numLayers;
                   end
                end
            end
        end
        
        % Handle number of SRS ports
        NumPUSCH = length(pusch);
        for j = 1:length(srs)
            if srs(j).Enable
                if srs(j).BWP == i
                    n(NumPUSCH+j) = srs(j).NumSRSPorts;
                end
            end
        end

        maxlayers(i) = max(n);
    end
    ResourceElementGrids = arrayfun(@(bp)zeros(bp.NRB*12,waveconfig.NumSubframes*symbolsPerSlot(bp)*fix(bp.SubcarrierSpacing/15),max(maxlayers)),...
        bwp,'UniformOutput',false);

    % Define channel plotting ID markers
    chplevel.PUCCH = 1.3;
    chplevel.PUSCH = 0.8;
    if ~isempty(srs)
        chplevel.SRS = 1.5;
    end
    
    % Get the overlapped slots of PUCCH and PUSCH of a specific RNTI
    [overlapSlots,rnti] = getOverlapSlots(pucch,pusch,bwp,waveconfig.NumSubframes);

    % Process the set of PUCCH transmission sequences
    for nch = 1:length(pucch)

        % Get a copy of the current PUCCH channel parameters
        ch = pucch(nch);

        % Only process configuration if enabled
        if ~ch.Enable
            continue;
        end

        % Check the referenced BWP and PUCCH indices
        checkIndex('PUCCH',pucch,nch,'BWP',bwp);

        % Establish whether transport coding is enabled
        uciCoding = ~isfield(ch,'EnableCoding') || ch.EnableCoding;

        % Get the number of symbols per slot for the associated BWP (CP dependent)
        symbperslot = symbolsPerSlot(bwp(ch.BWP));

        % Get the PUCCH resource parameters
        if ~ch.DedicatedResource && (symbperslot == 14) % Common resource is used only for normal CP

            % If dedicated resource is disabled and normal cyclic prefix,
            % the resource configuration parameters will be accessed based
            % on the resource index (Section 9.2.1, TS 38.213)
            % The format of the PUCCH common resource is 0 or 1.
            resInfo = pucchCommonResource(ch.ResourceIndex,bwp(ch.BWP).NRB);

            % Assign the multi-slot parameters and scrambling identity
            resInfo.NrOfSlots = ch.NrOfSlots;
            resInfo.InterSlotFreqHopping = ch.InterSlotFreqHopping;
            if ~isempty(ch.HoppingId)
                resInfo.HoppingId = ch.HoppingId;
            else
                % If NID is not configured, physical layer cell identity
                % NCellID is used
                resInfo.HoppingId = waveconfig.NCellID;
            end
        else
            % The configuration parameters provided in the PUCCH channel
            % will be accessed.
            resInfo = ch;
            if (ch.PUCCHFormat ~= 0 && ch.PUCCHFormat ~= 1) && isempty(ch.NID)
                % If NID is not configured, physical layer cell identity
                % NCellID is used. It is used for formats 2/3/4.
                resInfo.NID = waveconfig.NCellID;
            end
            if ch.PUCCHFormat ~=2 && isempty(ch.HoppingId)
                % If HoppingId is not configured, physical layer cell
                % identity NCellID is used.
                resInfo.HoppingId = waveconfig.NCellID;
            end
            if ch.PUCCHFormat == 2 && isempty(ch.NIDDMRS)
                % If NIDDMRS is not configured, physical layer cell
                % identity NCellID is used.
                resInfo.NIDDMRS = waveconfig.NCellID;
            end
        end

        % Check whether inter-slot frequency hopping is enabled for
        % repetition of slots
        interSlotFreqHopping = 0;
        if strcmpi(resInfo.InterSlotFreqHopping,'enabled') && resInfo.NrOfSlots > 1 ...
            && (resInfo.PUCCHFormat == 1 || resInfo.PUCCHFormat == 3 || resInfo.PUCCHFormat == 4)
            resInfo.IntraSlotFreqHopping = 'disabled';
            interSlotFreqHopping = 1;
        end

        % Ensure the number of OFDM symbols allocated for the specific
        % PUCCH format doesn't exceed the possible symbol allocation.
        if resInfo.PUCCHFormat == 0 || resInfo.PUCCHFormat == 2
            % For format 0 and 2, the maximum number of OFDM symbols that
            % can be allocated is 2.
            if resInfo.NrOfSymbols < 1 || resInfo.NrOfSymbols > 2
                error('The number of symbols allocated for PUCCH %d format %d must be either 1 or 2',nch,resInfo.PUCCHFormat);
            end
        else
            % For formats 1/3/4, the maximum number of OFDM symbols that
            % can be allocated is 14 (normal CP) or 12 (extended CP).
            if resInfo.NrOfSymbols < 4 || resInfo.NrOfSymbols > symbperslot
                error('The number of symbols allocated for PUCCH %d format %d must be in range 4 and %d',nch,resInfo.PUCCHFormat,symbperslot);
            end
        end

        % Get the allocated PRB, symbols and slots
        allocatedPRB = prbValues(resInfo);
        allocatedSymbols = resInfo.StartSymbol:resInfo.StartSymbol+resInfo.NrOfSymbols-1;
        allocatedSlots = cell2mat(getControlAllocatedSlots(ch,bwp,waveconfig.NumSubframes));

        % Ensure the allocated PRB is within the bandwidth part
        if any(max(allocatedPRB) >= bwp(ch.BWP).NRB)
            error('The allocated PRB indices (0-based, largest value = %d) for PUCCH %d exceed the NRB (%d) for BWP %d.',max(allocatedPRB(:)), nch, bwp(ch.BWP).NRB, ch.BWP);
        end

        % Mark the PUCCH sequences for display in the BWP grids, for visualization purposes only
        if ~interSlotFreqHopping
            ResourceGrids{ch.BWP}(1+allocatedPRB(:,1)',...
                reshape(1+symbperslot*allocatedSlots + allocatedSymbols(1:floor(resInfo.NrOfSymbols/2))',1,[])) = chplevel.PUCCH;
            if strcmpi(resInfo.IntraSlotFreqHopping,'enabled')
                ResourceGrids{ch.BWP}(1+allocatedPRB(:,2)',...
                    reshape(1+symbperslot*allocatedSlots+allocatedSymbols(floor(resInfo.NrOfSymbols/2)+1:end)',1,[])) = chplevel.PUCCH;
            else
                ResourceGrids{ch.BWP}(1+allocatedPRB(:,1)',...
                    reshape(1+symbperslot*allocatedSlots+allocatedSymbols(floor(resInfo.NrOfSymbols/2)+1:end)',1,[])) = chplevel.PUCCH;
            end
        else
            ResourceGrids{ch.BWP}(1+allocatedPRB(:,1)',...
                reshape(1+symbperslot*allocatedSlots(1:2:end) + allocatedSymbols',1,[])) = chplevel.PUCCH;
            ResourceGrids{ch.BWP}(1+(allocatedPRB(:,2)+resInfo.SecondHopPRB)',...
                reshape(1+symbperslot*allocatedSlots(2:2:end) + allocatedSymbols',1,[])) = chplevel.PUCCH;
        end

        % Create a data source for this PUCCH sequence
        datasource = hVectorDataSource(ch.DataSource); 

        % Loop over all the PUCCH transmission occasions and write the encoded
        % UCI payloads into the resource elements of the associated PUCCH instances
        count = 0;
        for index = 1:length(allocatedSlots)

            % Check for overlapping slots of PUCCH and PUSCH, ignore the
            % PUCCH transmission for the overlapped slot and transmit the
            % UCI information of PUCCH on PUSCH.
            if any(allocatedSlots(index) == overlapSlots{ch.BWP,ch.RNTI == rnti{ch.BWP}})
                count = count+1;
                continue;
            end

            % Get the slot-oriented PUCCH indices, DM-RS indices and DM-RS symbol values
            nslot = mod(allocatedSlots(index),10*(bwp(ch.BWP).SubcarrierSpacing/15));      % Slot number in the radio frame
            resInfo.NSlot = nslot;
            [indicesInfo, dmrssym] = ...
                hPUCCHResources(struct('CyclicPrefix',bwp(ch.BWP).CyclicPrefix,...
                'RBOffset',carriers(bwp(ch.BWP).CarrierIdx).RBStart+bwp(ch.BWP).RBOffset,'SubcarrierSpacing',bwp(ch.BWP).SubcarrierSpacing),resInfo);

            % Get the HARQ-ACK information based on LenACK field
            if ch.LenACK > 0
                ack = datasource.getPacket(ch.LenACK);
            else
                ack = [];
            end

            % Get the SR information based on LenSR field
            if ch.LenSR > 0
                sr = datasource.getPacket(ch.LenSR);
            else
                sr = [];
            end

            if uciCoding
                % Encode UCI bits for PUCCH formats 2/3/4
                if resInfo.PUCCHFormat == 2 || resInfo.PUCCHFormat == 3 || resInfo.PUCCHFormat == 4

                    if ch.LenCSI1
                        % Get the CSI part 1 payload
                        csi1 = datasource.getPacket(ch.LenCSI1);

                        % Combine HARQ-ACK, SR and CSI part1 information to UCI part 1
                        uciBits = [ack;sr;csi1];

                        % Length of UCI part 1
                        lenUCI1 = length(uciBits);
                        if lenUCI1 <= 2
                            error('The number of UCI bits (HARQ-ACK and SR and CSI) transmitted on PUCCH %d formats (2/3/4) must be more than 2.',nch);
                        end
                        L = getCRC(lenUCI1);         % Get the CRC length for UCI part 1

                        % Get CSI part 2 bits, if present
                        if ch.LenCSI2
                            % If the CSI part 2 is less than 3 bits, then the zeros
                            % are appended to make the length as 3 and assign to
                            % UCI part 2
                            uciBits2 = [datasource.getPacket(ch.LenCSI2); zeros(3-ch.LenCSI2,1)];
                            L2 = getCRC(length(uciBits2));
                        end

                        % Get the rate matching value for each UCI part
                        G1 = indicesInfo.G;
                        G2 = 0;
                            % UCI multiplexing happens for format 3 and 4,
                            % provided no repetition of slots
                        if ch.LenCSI2 && (resInfo.PUCCHFormat == 3 || resInfo.PUCCHFormat == 4) && ~(resInfo.NrOfSlots > 1)
                            qm = sum(strcmpi(resInfo.Modulation,{'pi/2-BPSK','QPSK'}).*[1 2]);
                            G1 = min(indicesInfo.G,ceil((lenUCI1 + L)/ch.MaxCodeRate/qm)*qm);
                            G2 = indicesInfo.G - G1;
                        end

                        % Check if the input bits are greater than the bit
                        % capacity of PUCCH resource
                        if lenUCI1+L >= G1
                            error('The sum of number of UCI part 1 (HARQ-ACK, SR, CSI part 1) bits (%d) and the CRC bits (%d) must be less than the bit capacity (%d) for PUCCH %d.',lenUCI1,L,G1,nch);
                        end

                        if G2
                            if ch.LenCSI2+L2 >= G2
                                error('The sum of number of CSI part 2 bits (%d) and the CRC bits (%d) must be less than the bit capacity (%d) for PUCCH %d.',ch.LenCSI2,L2,G2,nch);
                            end
                        end

                        % Encode the UCI payload to match the PUCCH bit
                        % capacity
                        codedUCI1 = nrUCIEncode(uciBits,G1);
                        if G2 && resInfo.PUCCHFormat ~= 2
                            % Encode UCI part 2, for format 3 and 4
                            codedUCI2 = nrUCIEncode(uciBits2,G2);

                            % Multiplex the encoded UCI part 1 and UCI part 2,
                            % assign to a codeword.
                            if resInfo.PUCCHFormat == 3
                                nRBOrSf = resInfo.NrOfRB;
                            else
                                nRBOrSf = resInfo.SpreadingFactor;
                            end
                            dmrsSymInd = indicesInfo.AllocatedDMRSSymbols-resInfo.StartSymbol;
                            codeword = hMultiplexPUCCH(codedUCI1,codedUCI2,resInfo.Modulation,resInfo.NrOfSymbols,...
                                dmrsSymInd,resInfo.PUCCHFormat,nRBOrSf);
                        else
                            % If no UCI part 2, codeword is the encoded UCI
                            % part 1.
                            codeword = codedUCI1;
                        end
                    else
                        % If no CSI part 1 on formats 2/3/4, codeword contains
                        % UCI information of HARQ-ACK and SR provided the input
                        % length is greater than 2 bits.
                        if length(ack) + length(sr) > 2
                            uciBits = [ack;sr];
                            lenUCI1 = length(uciBits);
                            L = getCRC(lenUCI1);
                            if lenUCI1+L >= indicesInfo.G
                                error('The sum of number of UCI (HARQ-ACK, SR) bits (%d) and the CRC bits (%d) must be less than the bit capacity (%d) for PUCCH %d.',lenUCI1,L,indicesInfo.G,nch);
                            end
                            codeword = nrUCIEncode(uciBits,indicesInfo.G);
                        else
                            error('The number of UCI bits (HARQ-ACK and SR) transmitted on PUCCH (%d) formats (2/3/4) must be more than 2.',nch);
                        end
                    end
                end
            else
                % Get the PUCCH codeword directly from the data source
                codeword = datasource.getPacket(indicesInfo.G);
            end

            % Process the UCI codeword on PUCCH
            switch resInfo.PUCCHFormat
                case 0
                    % PUCCH format 0
                    hoppingId = resInfo.HoppingId;           % Scrambling identity (hoppingId or NCellID, depending)
                    symbols = nrPUCCH0(ack,sr,[resInfo.StartSymbol resInfo.NrOfSymbols],...
                        bwp(ch.BWP).CyclicPrefix,nslot,hoppingId,resInfo.GroupHopping,resInfo.InitialCS,resInfo.IntraSlotFreqHopping);
                case 1
                    % PUCCH format 1
                    hoppingId = resInfo.HoppingId;           % Scrambling identity (hoppingId or NCellID, depending)
                    symbols = nrPUCCH1(ack,sr,[resInfo.StartSymbol resInfo.NrOfSymbols],...
                        bwp(ch.BWP).CyclicPrefix,nslot,hoppingId,ch.GroupHopping,resInfo.InitialCS,resInfo.IntraSlotFreqHopping,resInfo.OCCI);
                case 2
                    % PUCCH format 2
                    symbols = nrPUCCH2(codeword,resInfo.NID,ch.RNTI);
                case 3
                    % PUCCH format 3
                    symbols = nrPUCCH3(codeword,resInfo.Modulation,resInfo.NID,ch.RNTI,resInfo.NrOfRB);
                otherwise
                    % PUCCH format 4
                    symbols = nrPUCCH4(codeword,resInfo.Modulation,resInfo.NID,ch.RNTI,resInfo.SpreadingFactor,resInfo.OCCI);
            end
            symbols = reshape(symbols,[],length(indicesInfo.AllocatedSymbols)); % Reshape the column vector into the number of OFDM symbols for UCI information

            % Combine PUCCH with existing grid
            offset = 1+(reshape(symbperslot*allocatedSlots(index) + indicesInfo.AllocatedSymbols',1,[])*12*bwp(ch.BWP).NRB);
            repucch = indicesInfo.AllocatedSubcarriers;
            if interSlotFreqHopping
                count = mod(count,resInfo.NrOfSlots);
                if mod(count,2)==1
                    repucch = indicesInfo.AllocatedSubcarriers + resInfo.SecondHopPRB*12;
                end
            end
            ResourceElementGrids{ch.BWP}(offset+repucch) = ResourceElementGrids{ch.BWP}(offset+repucch) + symbols*db2mag(ch.Power); 

            % Combine PUCCH DM-RS with the grid
            offset = 1+(reshape(symbperslot*allocatedSlots(index) + indicesInfo.AllocatedDMRSSymbols',1,[])*12*bwp(ch.BWP).NRB);
            redmrs = indicesInfo.AllocatedDMRSSubcarriers;
            if interSlotFreqHopping
                count = mod(count,resInfo.NrOfSlots);
                if mod(count,2)==1
                    redmrs = indicesInfo.AllocatedDMRSSubcarriers + resInfo.SecondHopPRB*12;
                end
                count = count+1;
            end
            ResourceElementGrids{ch.BWP}(offset+redmrs) = ResourceElementGrids{ch.BWP}(offset+redmrs) + dmrssym*db2mag(ch.Power+ch.PowerDMRS);
        end
    % End of PUCCH sequence processing
    end

    % Process the set of PUSCH transmission sequences
    % Create a single UL-SCH channel processing object for use with all the PUSCH sequences
    ulsch = nrULSCH('MultipleHARQProcesses',false);

    for nch = 1:length(pusch)

        % Get a copy of the current PUSCH channel parameters
        ch = pusch(nch);

        % Only process configuration if enabled
        if ~ch.Enable
            continue;
        end

        % Establish whether transport coding is enabled
        ulschCoding = ~isfield(ch,'EnableCoding') || ch.EnableCoding;

        % Check the referenced BWP index
        checkIndex('PUSCH',pusch,nch,'BWP',bwp);

        % Get the allocated slots for PUSCH
        allocatedSlots = cell2mat(getULSCHAllocatedSlots(ch,bwp,waveconfig.NumSubframes));

        % Ensure the allocated PRB is within the bandwidth part
        ch.AllocatedPRB = sort(ch.AllocatedPRB);
        if any(ch.AllocatedPRB >= bwp(ch.BWP).NRB)
            error('The allocated PRB indices (0-based, largest value = %d) for PUSCH %d exceed the NRB (%d) for BWP %d.',max(ch.AllocatedPRB),nch,bwp(ch.BWP).NRB,ch.BWP);
        end
        if strcmpi(ch.IntraSlotFreqHopping,'enabled') || strcmpi(ch.InterSlotFreqHopping,'enabled')
            if isempty(ch.RBOffset) % Replace the empty value with 0
                ch.RBOffset = 0;
            end
            secondHopPRB = mod(ch.AllocatedPRB(1)+ch.RBOffset,bwp(ch.BWP).NRB);
            secondHopPRBSet = secondHopPRB+ch.AllocatedPRB-ch.AllocatedPRB(1);
            ch.RBOffset = secondHopPRB; % Replace RB offset with the calculated secondHopPRB to pass into hPUSCHResources
        end
        % Ensure that the allocated symbols for the slot are within a slot for the BWP CP
        symbperslot = symbolsPerSlot(bwp(ch.BWP));
        slotsymbs = ch.AllocatedSymbols(ch.AllocatedSymbols < symbperslot);
        if length(slotsymbs) ~= length(ch.AllocatedSymbols)
            warning('The slot-wise symbol allocation for PUSCH %d in BWP %d includes 0-based symbol indices which fall outside a slot (0...%d). Using only symbols within a slot.',nch,ch.BWP,symbperslot-1);
            ch.AllocatedSymbols = slotsymbs;
        end

        % Ensure the allocated symbols are contiguous
        if any(diff(ch.AllocatedSymbols) ~= 1)
            error('The allocated symbols for PUSCH %d must be contiguous.',nch);
        end

        % Ensure the PUSCH allocation starts with symbol 0 and has minimum
        % length of 4 for mapping type 'A', as per TS 38.214, Section 6.1.2
        if strcmpi(ch.PUSCHMappingType,'A') 
            if ch.AllocatedSymbols(1) ~= 0
                error('For PUSCH mapping type A, the starting symbol must be 0.');
            end
            if length(ch.AllocatedSymbols) < 4
                error('For PUSCH mapping type A, the minimum number of allocated symbols must be 4.');
            end
        end

        % Display related PRB level processing
        %
        % Calculate the *PRB* linear indices of all the PUSCH instances, primarily
        % for display purposes here.
        % This is performed by marking the allocated PRB in an empty PRB grid 
        % for the BWP in the entire waveform period, subtracting out the reserved
        % part then find the indices that have been used

        % Create an empty BWP spanning the length of the waveform
        rgrid = zeros(size(ResourceGrids{ch.BWP}));

        % Check if inter-slot frequency hopping is enabled
        interSlotFreqHopping = 0;
        if strcmpi(ch.InterSlotFreqHopping,'enabled')
            ch.IntraSlotFreqHopping = 'disabled';
            interSlotFreqHopping = 1;
        end

        % Mark the PRB/symbols associated with all the PUSCH instances in this sequence
        if ~interSlotFreqHopping
            if ~strcmpi(ch.IntraSlotFreqHopping,'enabled')
                for ns=allocatedSlots
                  rgrid(1+ch.AllocatedPRB,1+symbperslot*ns+ch.AllocatedSymbols) = 1;
                end
            else
                for ns=allocatedSlots
                  rgrid(1+ch.AllocatedPRB,1+symbperslot*ns+ch.AllocatedSymbols(1:floor(length(ch.AllocatedSymbols)/2))) = 1;
                  rgrid(1+secondHopPRBSet,1+symbperslot*ns+ch.AllocatedSymbols(floor(length(ch.AllocatedSymbols)/2)+1:end)) = 1;
                end
            end
        else
            for ns=allocatedSlots
                if mod(ns,2)==0
                  rgrid(1+ch.AllocatedPRB,1+symbperslot*ns+ch.AllocatedSymbols) = 1;
                else
                  rgrid(1+secondHopPRBSet,1+symbperslot*ns+ch.AllocatedSymbols) = 1;
                end
            end
        end

        % Identify all the indices that remain
        puschindices = find(rgrid);

        % Mark the used PUSCH locations in the PRB grid
        ResourceGrids{ch.BWP}(puschindices) = ResourceGrids{ch.BWP}(puschindices)+chplevel.PUSCH;

        % Waveform generation RE level processing
        %
        % The hPUSCHResources uses a slot-level set of parameters so map the
        % relevant parameter from the waveform level down to the slot level
        nrb = bwp(ch.BWP).NRB;
        ch.PRBSet = ch.AllocatedPRB;
        ch.SymbolSet = ch.AllocatedSymbols;
        ch.PRBRefPoint = carriers(bwp(ch.BWP).CarrierIdx).RBStart + bwp(ch.BWP).RBOffset;

        % Update the number of layers based on number of DM-RS ports when
        % PortSet is provided
        if isfield(ch,'PortSet') && ~isempty(ch.PortSet)
            ch.NLayers = numel(unique(ch.PortSet));
        end

        % Check for scrambling identities and configure physical layer cell
        % identity NCellID, if necessary
        if isempty(ch.NID)
            ch.NID = waveconfig.NCellID;
        end
        if isempty(ch.NIDNSCID)
            ch.NIDNSCID = waveconfig.NCellID;
        end
        if isempty(ch.NRSID)
            ch.NRSID = waveconfig.NCellID;
        end

        % Create a data source for this PUSCH sequence
        datasource = hVectorDataSource(ch.DataSource);

        % Configure the UL-SCH processing object for this PUSCH sequence
        if ulschCoding
            ulsch.TargetCodeRate = ch.TargetCodeRate;
        end

        % Get the PUCCH structures within the same BWP
        if ~isempty(overlapSlots{ch.BWP,ch.RNTI == rnti{ch.BWP}})
            id = cell2mat(arrayfun(@(x) getfield(pucch,{x},'BWP'),1:length(pucch),'UniformOutput',false)) == ch.BWP;
            for i = 1:length(id)
                if id(i) && pucch(i).Enable
                    if pucch(i).RNTI == ch.RNTI
                        pucchId = i;
                        pucchSlots = cell2mat(getControlAllocatedSlots(pucch(pucchId),bwp(pucch(pucchId).BWP),waveconfig.NumSubframes));
                        if any(any(reshape(pucchSlots,1,[]) == reshape(overlapSlots{ch.BWP,ch.RNTI == rnti{ch.BWP}},[],1)))
                            break;
                        end
                    end
                end
            end
        end

        % Storage for PUSCH instance information
        datastore = [];
        
        % Initialize modinfo
        modinfo = struct('Gd',0,'G',0,'NREPerPRB',0,'DMRSSymbolSet',[],'CDMGroups',[],'PTRSSymbolSet',[],'CDMLengths',[]);
            
        % Loop over all the allocated slots
        for i = 1:length(allocatedSlots)

            % Get current slot number
            s = allocatedSlots(i);
            ch.NSlot = s;

            % Create an empty slot grid to contain a single PUSCH instance
            slotgrid = zeros(12*nrb,symbperslot,maxlayers(ch.BWP));

            % Get the slot-oriented PUSCH indices, DM-RS indices and DM-RS symbol values  
            [puschREindices,dmrsREindices,dmrsSymbols,ptrsREindices,ptrsSymbols,modinfo] = ...
                hPUSCHResources(struct('NRB',bwp(ch.BWP).NRB,'CyclicPrefix',bwp(ch.BWP).CyclicPrefix,'SubcarrierSpacing',bwp(ch.BWP).SubcarrierSpacing),ch);

            if ulschCoding
                % Get the RV value for this transmission instance
                rvidx = mod(i-1,length(ch.RVSequence))+1;
                rv = ch.RVSequence(rvidx);

                % For the first RV in a sequence, get a new transport block from 
                % the data source and pass it to the UL-SCH processing
                if rvidx == 1
                   trblksize = nrTBS(ch.Modulation,ch.NLayers,length(ch.PRBSet),modinfo.NREPerPRB,ch.TargetCodeRate,ch.Xoh_PUSCH);
                   trblk = datasource.getPacket(trblksize);
                   setTransportBlock(ulsch,trblk);
                end

                % Overlap slots with PUCCH - transmit PUSCH with UCI
                if any(allocatedSlots(i)==overlapSlots{ch.BWP,ch.RNTI == rnti{ch.BWP}})

                    if strcmpi(ch.PUSCHMappingType,'B') &&...
                            ((strcmpi(ch.IntraSlotFreqHopping,'enabled') && length(ch.AllocatedSymbols) < 3) ||...
                            length(ch.AllocatedSymbols) == 1)
                        disp(['The UCI information is ignored in the overlapping slot number ', num2str(allocatedSlots(i)),...
                            ' for PUSCH ',num2str(nch),' mapping type B configuration.']);
                        codeword = ulsch(ch.Modulation,ch.NLayers,modinfo.G,rv);

                    else
                        % Get the lengths of UCI bits of the overlapped PUCCH
                        uciInfo.OACK = pucch(pucchId).LenACK;
                        uciInfo.OCSI1 = pucch(pucchId).LenCSI1;
                        uciInfo.OCSI2 = pucch(pucchId).LenCSI2;

                        if ~uciInfo.OCSI1 && uciInfo.OCSI2
                            % When there is no CSI part 1, set CSI part 2
                            % to zero
                            uciInfo.OCSI2 = 0;
                        end

                        % Get the datasource of UCI bits for the overlapped
                        % PUCCH
                        datasourcePUCCH = hVectorDataSource(pucch(pucchId).DataSource);

                        % Get the proper UCI payload for the specified
                        % configuration
                        if ch.DisableULSCH
                            % Set UL-SCH block size to 0
                            trblksize = 0;
                        end
                        % Get the pusch object from the input structure
                        puschObj = getPUSCHObject(ch);
                        rmInfo = nrULSCHInfo(puschObj,ch.TargetCodeRate,trblksize,uciInfo.OACK,uciInfo.OCSI1,uciInfo.OCSI2);

                        if pucch(pucchId).EnableCoding
                            % Perform UCI encoding based on the
                            % configuration of PUCCH
                            if uciInfo.OACK
                                ack = uciLen(uciInfo.OACK,rmInfo.GACK,puschObj.Modulation,'HARQ-ACK',datasourcePUCCH,nch);
                                uciInfo.OACK = length(ack);
                            else
                                ack = [];
                            end
                            if uciInfo.OCSI1
                                csi1 = uciLen(uciInfo.OCSI1,rmInfo.GCSI1,puschObj.Modulation,'CSI part 1',datasourcePUCCH,nch);
                            else
                                csi1 = [];
                            end
                            if uciInfo.OCSI2
                                csi2 = uciLen(uciInfo.OCSI2,rmInfo.GCSI2,puschObj.Modulation,'CSI part 2',datasourcePUCCH,nch);
                            else
                                csi2 = [];
                            end

                            % Encode UCI
                            codedACK  = nrUCIEncode(ack,rmInfo.GACK,puschObj.Modulation);
                            codedCSI1 = nrUCIEncode(csi1,rmInfo.GCSI1,puschObj.Modulation);
                            codedCSI2 = nrUCIEncode(csi2,rmInfo.GCSI2,puschObj.Modulation);
                        else
                            % If UCI coding is disabled and UL-SCH coding
                            % is enabled, get the coded bits directly from
                            % the data source.
                            if rmInfo.GACK
                                codedACK = datasource.getPacket(rmInfo.GACK);
                            else
                                codedACK = [];
                            end
                            if rmInfo.GCSI1
                                codedCSI1 = datasource.getPacket(rmInfo.GCSI1);
                            else
                                codedCSI1 = [];
                            end
                            if rmInfo.GCSI2
                                codedCSI2 = datasource.getPacket(rmInfo.GCSI2);
                            else
                                codedCSI2 = [];
                            end
                        end

                        % Encode UL-SCH
                        if ~ch.DisableULSCH
                            codedULSCH = ulsch(puschObj.Modulation,puschObj.NumLayers,rmInfo.GULSCH,rv);
                        else
                            codedULSCH = [];
                        end

                        % Multiplex UL-SCH and UCI to create a codeword
                        codeword = nrULSCHMultiplex(puschObj,ch.TargetCodeRate,trblksize,codedULSCH,codedACK,codedCSI1,codedCSI2);

                        % Display the message
                        disp(['The UCI information is transmitted on PUSCH for slot number ', num2str(allocatedSlots(i)),...
                            ' due to overlap of PUCCH ',num2str(pucchId),' and PUSCH ',num2str(nch),' in the bandwidth part ',num2str(ch.BWP),' for RNTI ',num2str(ch.RNTI),'.']);
                    end
                else
                    % UL-SCH processing to create a codeword
                    codeword = ulsch(ch.Modulation,ch.NLayers,modinfo.G,rv);
                end
            else
                % If transport coding is not enabled, get the codeword
                % directly from the data source
                codeword = datasource.getPacket(modinfo.G);
                rv = [];
                trblk = [];
            end

            % PUSCH processing to create the PUSCH symbols
            nID = ch.NID;
            nLayers = ch.NLayers;
            if ~ch.TransformPrecoding || isempty(ptrsREindices)
                symbols = nrPUSCH(codeword,ch.Modulation,nLayers,nID,ch.RNTI,ch.TransformPrecoding,length(ch.PRBSet),ch.TxScheme,ch.NAntennaPorts,ch.TPMI);
            else
                % Scrambling
                scrambled = nrPUSCHScramble(codeword,nID,ch.RNTI);
                % Symbol modulation
                modulated = nrSymbolModulate(scrambled,ch.Modulation);
                % Layer mapping
                layered = nrLayerMap(modulated,nLayers);
                % PT-RS insertion with PUSCH data
                puschSymbolsPerLayer = modinfo.Gd+ size(ptrsREindices,1);
                layeredWithPTRS = zeros(puschSymbolsPerLayer,nLayers,'like',layered);
                layeredWithPTRS(ptrsREindices) = ptrsSymbols;
                puschInd = setdiff((1:nLayers*puschSymbolsPerLayer)',ptrsREindices);
                layeredWithPTRS(puschInd) = layered;
                % Transform precoding
                transformed = nrTransformPrecode(layeredWithPTRS,length(ch.PRBSet));
                % MIMO precoding
                if strcmpi(ch.TxScheme,'codebook')
                     W = nrPUSCHCodebook(nLayers,ch.NAntennaPorts,ch.TPMI,ch.TransformPrecoding);
                else
                    W = eye(nLayers);
                end
                symbols = transformed * W;
            end

            % Write the PUSCH, DM-RS (precoded) and PT-RS (precoded)
            % symbols in the slot grid PUSCH
            if ~isempty(symbols)
                slotgrid(puschREindices) = symbols*db2mag(ch.Power);
            end
            slotgrid(dmrsREindices) = dmrsSymbols*db2mag(ch.Power+ch.PowerDMRS);
            ptpower = 1;
            if ~ch.TransformPrecoding && ~isempty(ptrsREindices)
                ptpower = db2mag(ch.Power+ch.PowerPTRS);
                slotgrid(double(ptrsREindices)) = ptrsSymbols*ptpower;
            end

            % Combine PUSCH instance with the rest of the BWP grid
            ResourceElementGrids{ch.BWP}(:,s*symbperslot+(1:symbperslot),1:maxlayers(ch.BWP)) = ResourceElementGrids{ch.BWP}(:,s*symbperslot+(1:symbperslot),1:maxlayers(ch.BWP)) + slotgrid; 

            % Capture resource info for this PUSCH instance
            datastore(i).NSlot = ch.NSlot; %#ok<AGROW>
            datastore(i).TransportBlockSize = length(trblk); %#ok<AGROW>
            datastore(i).TransportBlock = trblk;  %#ok<AGROW>
            datastore(i).RV = rv;  %#ok<AGROW>
            datastore(i).Codeword = codeword;  %#ok<AGROW>
            datastore(i).G = modinfo.G; %#ok<AGROW>
            datastore(i).Gd = modinfo.Gd; %#ok<AGROW>
            datastore(i).ChannelIndices = puschREindices; %#ok<AGROW>
            datastore(i).ChannelSymbols = symbols*db2mag(ch.Power); %#ok<AGROW>
            datastore(i).DMRSSymbolSet = modinfo.DMRSSymbolSet;  %#ok<AGROW>
            datastore(i).DMRSIndices = dmrsREindices; %#ok<AGROW>
            datastore(i).DMRSSymbols = dmrsSymbols*db2mag(ch.Power+ch.PowerDMRS); %#ok<AGROW>
            datastore(i).PTRSSymbolSet = modinfo.PTRSSymbolSet;  %#ok<AGROW>
            datastore(i).PTRSIndices = ptrsREindices; %#ok<AGROW>
            datastore(i).PTRSSymbols = ptrsSymbols*ptpower; %#ok<AGROW>         
                
        end

        % Capture all resources info for this PUSCH sequence
        name = "";
        if isfield(pusch(nch),'Name')
            name = pusch(nch).Name;
        end
        waveinfo.PUSCH(nch).Name = name;
        waveinfo.PUSCH(nch).CDMLengths = modinfo.CDMLengths;
        waveinfo.PUSCH(nch).Resources = datastore;
        
    % End of PUSCH sequence processing
    end
    
    % Process the set of SRS transmissions     
    if ~isempty(srs)
        
        for nch = 1:length(srs)       
            % Get a copy of the current SRS and BWP parameters
            ch = srs(nch);
            bp = bwp(ch.BWP);        

            % Only process configuration if SRS is enabled
            if ~ch.Enable
                continue;
            end
            
            SRSAllocatedSlots = expandbyperiod(ch.AllocatedSlots,ch.AllocatedPeriod,waveconfig.NumSubframes,bp.SubcarrierSpacing);

            % Check the referenced BWP index
            checkIndex('SRS',srs,nch,'BWP',bwp);

            symbperslot = symbolsPerSlot(bp);

            % Display-related PRB and waveform generation RE level processing
            %
            % Create an empty BWP grid spanning the length of the waveform
            rbgrid = zeros(size(ResourceGrids{ch.BWP}));

            for s = SRSAllocatedSlots            
                % Get the carrier configuration object
                carrierObj = getCarrierObject(bp);

                % Create an empty slot grid to contain a single SRS instance
                slotgrid = nrResourceGrid(carrierObj,maxlayers(ch.BWP));

                % Get the SRS configuration object
                srsObj = getSRSObject(ch);

                % Set slot number
                carrierObj.NSlot = s;

                % Generate SRS indices and sequence 
                srsInd = nrSRSIndices(carrierObj,srsObj);
                slotgrid(srsInd) = nrSRS(carrierObj,srsObj);

                % Combine SRS instance with the rest of the RE carrier grid
                ResourceElementGrids{ch.BWP}(:,s*symbperslot+(1:symbperslot),1:maxlayers(ch.BWP)) = ResourceElementGrids{ch.BWP}(:,s*symbperslot+(1:symbperslot),1:maxlayers(ch.BWP)) + slotgrid; 

                % Mark the PRB/symbols associated with all the SRS instances in this sequence
                rbgrid(unique( ceil( double(srsInd(1:end/ch.NumSRSPorts) + bp.NRB *12 *symbperslot*s )/12 ) )) = 1;            
            end
            % Identify the indices of the marked RBs
            srsInd = find(rbgrid);

            % Mark the used SRS locations in the PRB grid
            ResourceGrids{ch.BWP}(srsInd) = ResourceGrids{ch.BWP}(srsInd)+chplevel.SRS;
        end 

        % End of SRS sequence processing
    end
    
    % Create a new figure to display the plots
    % Map the BWPs into carrier sized PRB grids for display
    if waveconfig.DisplayGrids
        figure;
        for bp = 1:length(ResourceGrids)

            % Mark the unused RE in the overall BWP, relative to the
            % carrier, so that it is easier to see with respect to the
            % complete carrier layout
            bgrid = ResourceGrids{bp};
            cgrid = zeros(carriers((bwp(bp).CarrierIdx)).NRB, size(bgrid,2));
            bgrid(bgrid==0) = 0.15;

            % Write the BWP into the grid representing the carrier
            cgrid(bwp(bp).RBOffset + (1:size(bgrid,1)),:) = bgrid;

            % Plot the PRB BWP grid (relative to the carrier)
            cscaling = 40;
            subplot(length(ResourceGrids),1,bp)
            im = image(cscaling*cgrid); axis xy; title(sprintf('BWP %d in Carrier (SCS=%dkHz). PUCCH, PUSCH and SRS locations',bp,bwp(bp).SubcarrierSpacing)); xlabel('Symbols'); ylabel('Carrier RB');
            cmap = parula(64);
            colormap(im.Parent,cmap);
            % Add a channel legend to the first BWP plot (applies to all)
            if bp == 1
                % Extract channel names and color marker levels
                fnames = fieldnames(chplevel);
                chpval = struct2cell(chplevel);
                clevels = cscaling*[chpval{:}];
                N = length(clevels);
                L = line(ones(N),ones(N), 'LineWidth',8);                   % Generate lines
                % Index the color map and associated the selected colors with the lines
                set(L,{'color'},mat2cell(cmap( min(1+clevels,length(cmap) ),:),ones(1,N),3));   % Set the colors according to cmap
                % Create legend
                legend(fnames{:});
            end

        end
    end

    % Initialize output variables for the baseband waveform and info structure
    waveform = 0;
    gridset = struct('ResourceGridBWP',{},'ResourceGridInCarrier',{},'Info',{});

    % Establish the maximum carrier SCS configured and the associated k0 subtrahend 
    [maxcarrierscs,maxidx] = max(carrierscs); %#ok<ASGLU>
    k0offset = (carriers(maxidx).RBStart + carriers(maxidx).NRB/2)*12*(carriers(maxidx).SubcarrierSpacing/15);

    % Calculate the maximum OFDM sampling rate used across the configured SCS carriers
    sr = @(x)getfield(nrOFDMInfo(x.NRB,x.SubcarrierSpacing),'SampleRate');
    maxsr = max(arrayfun(sr,carriers));

    % Modulate all the BWP grids and combine all into a single baseband waveform matrix
    for bp = 1:length(ResourceElementGrids)

        % Get the current BWP RE grid
        bgrid = ResourceElementGrids{bp};

        % Configure a carrier config object associated with the BWP numerology
        carrierObj = nrCarrierConfig;
        carrierObj.NCellID = waveconfig.NCellID;
        carrierObj.SubcarrierSpacing = bwp(bp).SubcarrierSpacing;
        carrierObj.CyclicPrefix = bwp(bp).CyclicPrefix;
        carrierObj.NSizeGrid = carriers(bwp(bp).CarrierIdx).NRB;
        carrierObj.NStartGrid = carriers(bwp(bp).CarrierIdx).RBStart;
        carrierObj.NSlot = 0;

        % Check BWP dimensions relative to SCS carrier
        if (bwp(bp).RBOffset+bwp(bp).NRB) > carrierObj.NSizeGrid
            error('BWP %d (NRB = %d and RBOffset = %d @ %d kHz SCS) is outside of the SCS carrier bandwidth (NRB = %d).',bp,bwp(bp).NRB,bwp(bp).RBOffset,bwp(bp).SubcarrierSpacing,carrierObj.NSizeGrid);
        end

        % Create empty SCS carrier grid and assign in the BWP
        cgrid = nrResourceGrid(carrierObj,size(bgrid,3));
        cgrid = repmat(cgrid,1,waveconfig.NumSubframes*carrierObj.SubcarrierSpacing/15);
        cgrid(12*bwp(bp).RBOffset + (1:size(bgrid,1)),:,:) = bgrid;

        % If the input structure waveconfig contains a Windowing field, use
        % it in the OFDM modulation. Use the default value for windowing
        % otherwise.
        windowing = [];
        if isfield(waveconfig,'Windowing')
            windowing = waveconfig.Windowing;
        end

        % Modulate the entire grid
        [bwpwave,minfo] = nrOFDMModulate(carrierObj,cgrid,'Windowing',windowing,'SampleRate',maxsr);

        % Apply numerology dependent k0 offset, if required
        k0 = (carrierObj.NStartGrid + carrierObj.NSizeGrid/2)*12 - (k0offset/(carrierObj.SubcarrierSpacing/15));
        if k0~=0
           t = (0:size(bwpwave,1)-1)' / minfo.SampleRate;
           bwpwave = bwpwave .* exp(1j*2*pi*k0*carrierObj.SubcarrierSpacing*1e3*t);
        end

        % Add additional fields to the info
        minfo.NSubcarriers = size(cgrid,1);
        minfo.SubcarrierSpacing = carrierObj.SubcarrierSpacing;
        minfo.SymbolsPerSubframe = minfo.SymbolsPerSlot*minfo.SlotsPerSubframe;
        minfo.SamplesPerSubframe = minfo.SampleRate/1000;
        minfo.SubframePeriod = minfo.SamplesPerSubframe/minfo.SampleRate;
        minfo.k0 = k0;
        minfo.SamplingRate = minfo.SampleRate;

        % Combine this BWP with the rest of the waveform
        waveform = waveform + bwpwave;

        % Capture the intermediate grids and modulation info
        gridset(bp).ResourceGridBWP = bgrid;
        gridset(bp).ResourceGridInCarrier = cgrid;
        gridset(bp).Info = minfo;

    end

    % Display *subcarrier* resource grids
    if waveconfig.DisplayGrids

        % Create a new figure to display the subcarrier plots
        figure;
        plotCarriers(waveconfig,gridset);

        % Create a new figure to display the subcarrier plots
        figure;
        % Modulate all the BWP grids and combine all into a single baseband waveform matrix
        cmap = parula(64);
        for bp = 1:length(ResourceElementGrids)
            % Plot the resource element grid (scaled complex magnitude)
            subplot(length(ResourceElementGrids),1,bp)
            im = image(40*abs(gridset(bp).ResourceGridInCarrier(:,:,1))); axis xy;
            colormap(im.Parent,cmap);
            title(sprintf('BWP %d in Carrier (SCS=%dkHz)',bp,bwp(bp).SubcarrierSpacing)); xlabel('Symbols'); ylabel('Subcarriers');
        end
    end
    
    winfo.WaveformResources = waveinfo;
    
% End of main function
end

% Expand 's' with respect to period 'd', up to length 'nsf', optional accounting for the SCS
function sp = expandbyperiod(s,p,nsf,scs)

    if nargin > 3
        % Expand s by period p for ts length
        ts = nsf*1*fix(scs/15);
    else
        ts = nsf;
    end
    % Is the period is empty then the pattern doesn't repeat, so doesn't need extending
    if isempty(p)
        sp = s;
    else
        sp = reshape(s(s<p),[],1)+p*(0:ceil(ts/p)-1);
    end
    sp = reshape(sp(sp < ts),1,[]);            % Trim any excess
end

% Establish the number of symbols per slot from the cyclic prefix
function symbs = symbolsPerSlot(config)

    if isstruct(config)
        config = config.CyclicPrefix;
    end
    cpoptions = {'Normal','Extended'};
    symbs = sum(strcmpi(config,cpoptions) .* [14 12]);

end

% Check a set of named parameter index values against the length of the things 
% that they index
function checkIndex(chname,chset,nch,varargin)
    for p = 1:2:length(varargin)
        paramname = varargin{p};
        paramset = varargin{p+1};
        chindex = chset(nch).(paramname);
        plength = length(paramset);
        if (chindex < 1) || (chindex > plength)
            error('For %s %d, the %s index (%d) must be between 1 and the number of %s defined (%d)',chname,nch,paramname,chindex,paramname,plength);
        end
    end
end

function info = pucchCommonResource(resourceIndex,nSizeBWP)

%   INFO = pucchCommonResource(RESOURCEINDEX,NSIZEBWP) returns
%   the structural information of PUCCH resources for normal cyclic prefix
%   based on the resource index RESOURCEINDEX and size of bandwidth part
%   NSIZEBWP. The structure INFO contains the following fields:
%   PUCCHFormat  - PUCCH format configuration (0 or 1)
%   SymIndex     - Starting symbol index of PUCCH allocation in a slot
%   nPUCCHSym    - Number of OFDM symbols allocated for PUCCH
%   RBOffset     - PRB offset in the BWP
%   StartPRB     - Index of first PRB prior to frequency hopping
%   SecondHopPRB - Index of first PRB after frequency hopping
%   InitialCS    - Initial cyclic shift
%   OCCI         - Orthogonal cover code index. 0 for format 1, [] for
%                  format 0.

    % PUCCH Format number
    if resourceIndex < 3
        format = 0;
        occi = [];
    else
        format = 1;
        occi = 0;
    end

    % Symbol allocation
    if resourceIndex < 3
        symIndex = 12;
        nPUCCHSym = 2;
    elseif resourceIndex < 7
        symIndex = 10;
        nPUCCHSym = 4;
    elseif resourceIndex < 11
        symIndex = 4;
        nPUCCHSym = 10;
    else
        symIndex = 0;
        nPUCCHSym = 14;
    end

    % Initial cyclic shift
    csIndexesTable = {{0 3},{0 4 8},{0 4 8},{0 6},{0 3 6 9},{0 3 6 9},...
                      {0 3 6 9},{0 6},{0 3 6 9},{0 3 6 9},{0 3 6 9},{0 6},...
                      {0 3 6 9},{0 3 6 9},{0 3 6 9},{0 3 6 9}};
    nCSList = arrayfun(@(x)(length(csIndexesTable{x})), 1:length(csIndexesTable), 'UniformOutput', false);
    nCS = nCSList{resourceIndex+1};
    initialCS = csIndexesTable{resourceIndex+1}{mod(mod(resourceIndex,8),nCS)+1};

    % Resource block index for each hop
    rbOffsetList = [0 0 3 0 0 2 4 0 0 2 4 0 0 2 4 floor(nSizeBWP/4)];
    rbOffset = rbOffsetList(resourceIndex+1);
    if resourceIndex < 8
        startPRB = rbOffset + floor(resourceIndex/nCS);
        secondHopPRB = nSizeBWP-1-rbOffset-floor(resourceIndex/nCS);
    else
        startPRB = nSizeBWP-1-rbOffset-floor((resourceIndex-8)/nCS);
        secondHopPRB = rbOffset + floor((resourceIndex-8)/nCS);
    end

    % Combine information into a structure
    info.PUCCHFormat = format;
    info.StartSymbol = symIndex-1; % 0-based
    info.NrOfSymbols = nPUCCHSym;
    info.RBOffset = rbOffset;
    info.GroupHopping = 'enable';
    info.IntraSlotFreqHopping = 'enabled';
    info.StartPRB = startPRB;
    info.SecondHopPRB = secondHopPRB;
    info.InitialCS = initialCS;
    info.OCCI = occi;
end

function cas = getControlAllocatedSlots(pucch,bwp,NumSubframes)
    % Allocated slots for control information (PUCCH)
    cas = cell(length(pucch),1);
    for i = 1:length(pucch)
        ch = pucch(i);
        if ch.Enable
            if ch.NrOfSlots > 1 ...
                        && (ch.PUCCHFormat == 1 || ch.PUCCHFormat == 3 || ch.PUCCHFormat == 4)
                allocatedSlots = expandbyperiod(unique(ch.AllocatedSlots + (0:ch.NrOfSlots-1)')',ch.AllocatedPeriod,NumSubframes,bwp(ch.BWP).SubcarrierSpacing);
            else
                allocatedSlots = expandbyperiod(ch.AllocatedSlots,ch.AllocatedPeriod,NumSubframes,bwp(ch.BWP).SubcarrierSpacing);
            end
            cas{i} = allocatedSlots;
        end
    end
end

function uas = getULSCHAllocatedSlots(pusch,bwp,NumSubframes)
    % Allocated slots for UL-SCH information (PUSCH)
    uas = cell(length(pusch),1);
    for i = 1:length(pusch)
        ch = pusch(i);
        if ch.Enable
            uas{i} = expandbyperiod(ch.AllocatedSlots,ch.AllocatedPeriod,NumSubframes,bwp(ch.BWP).SubcarrierSpacing);
        end
    end
end

function [overlapSlots,rnti] = getOverlapSlots(pucch,pusch,bwp,NumSubframes)

    % Overlap slots of PUCCH and PUSCH in a BWP
    rnti = cell(length(bwp),1);
    for i = 1:length(bwp)
        % PUCCH indices associated with BWP 'i'
        cIndex = false(1,length(pucch));
        cRNTI = -1*ones(1,length(pucch));
        count = 0;
        for ipucch = 1:length(pucch)
            ch = pucch(ipucch);
            if ch.Enable
                cIndex(ipucch) = (ch.BWP == i);
                if cIndex(ipucch)
                    if ch.RNTI < 0
                        error('For PUCCH %d in BWP %d, the RNTI (%d) must be in range 0 to 65535.',ipucch,i,ch.RNTI);
                    end
                    cRNTI(ipucch) = ch.RNTI;
                    count = count+1;
                end
            end
        end
        % PUSCH indices associated with BWP 'i'
        uIndex = false(1,length(pusch));
        uRNTI = -1*ones(1,length(pusch));
        count = 0;
        for ipusch = 1:length(pusch)
            ch = pusch(ipusch);
            if ch.Enable
                uIndex(ipusch) = (ch.BWP == i);
                if uIndex(ipusch)
                    if ch.RNTI < 0
                        error('For PUSCH %d in BWP %d, the RNTI (%d) must be in range 0 to 65535.',ipusch,i,ch.RNTI);
                    end
                    uRNTI(ipusch) = ch.RNTI;
                    count = count+1;
                end
            end
        end
        rnti{i} = union(cRNTI,uRNTI);

        % Get the overlapping slots
        for j = 1:length(rnti{i})
            pucchIdx = (rnti{i}(j) == cRNTI);
            puschIdx = (rnti{i}(j) == uRNTI);

            if ~(any(pucchIdx) && any(puschIdx)) || (rnti{i}(j) == -1)
                overlapSlots{i,j} = [];    %#ok<AGROW>
            else
                cch = pucch(pucchIdx);
                cSlotsAll = [];
                for iCCH = 1:length(cch)
                    if cch(iCCH).NrOfSlots > 1 ...
                            && (cch(iCCH).PUCCHFormat == 1 || cch(iCCH).PUCCHFormat == 3 || cch(iCCH).PUCCHFormat == 4)
                        cSlots = expandbyperiod(unique(cch(iCCH).AllocatedSlots + (0:cch(iCCH).NrOfSlots-1)')',cch(iCCH).AllocatedPeriod,NumSubframes,bwp(cch(iCCH).BWP).SubcarrierSpacing);
                    else
                        cSlots = expandbyperiod(cch(iCCH).AllocatedSlots,cch(iCCH).AllocatedPeriod,NumSubframes,bwp(cch(iCCH).BWP).SubcarrierSpacing);
                    end
                    cSlotsAll = [cSlotsAll cSlots]; %#ok<AGROW>
                end
                sch = pusch(puschIdx);
                uSlotsAll = [];
                for iSCH = 1:length(sch)
                    uSlots = expandbyperiod(sch(iSCH).AllocatedSlots,sch(iSCH).AllocatedPeriod,NumSubframes,bwp(sch(iSCH).BWP).SubcarrierSpacing);
                    uSlotsAll = [uSlotsAll uSlots]; %#ok<AGROW>
                end
                overlapSlots{i,j} = intersect(cSlotsAll,uSlotsAll);   %#ok<AGROW>
            end
        end
    end
end

function L = getCRC(oUCI,fl)
% CRC bits for UCI information for input length oUCI (default).
% Also returns the number of CRC bits added when fl is 1, oUCI is treated
% as the rate matched length.
    if oUCI > 19
        L = 11;
    elseif oUCI > 11
        L = 6;
    else
        L = 0;
    end
    if nargin == 2 && fl
        if oUCI > 25
            L = 11;
        elseif oUCI > 17
            L = 6;
        else
            L = 0;
        end
    end
end

function prbIndices = prbValues(resInfo)

    % Check with intra-slot frequency hopping
    if (strcmpi(resInfo.IntraSlotFreqHopping,'enabled') && (isempty(resInfo.SecondHopPRB) || ~isfield(resInfo,'SecondHopPRB')))...
            || strcmpi(resInfo.IntraSlotFreqHopping,'disabled')
        resInfo.SecondHopPRB = resInfo.StartPRB;
    end

    % Initialize number of resource blocks
    nrOfRB = 1;
    if resInfo.PUCCHFormat == 2 || resInfo.PUCCHFormat == 3
        nrOfRB = resInfo.NrOfRB;
    end

    prbIndices = [resInfo.StartPRB resInfo.SecondHopPRB]+(0:nrOfRB-1)';
end

function v = uciLen(i,r,m,t,d,n)
%  Returns the UCI packet v, based on the input length i, rate match
%  length r, modulation scheme m, the type t and the data source d, when
%  UCI information is to be multiplexed on PUSCH n.

    if (i == 1 && strcmpi(m,'pi/2-BPSK')) ...
        && (i == r)
        v = d.getPacket(i);
    else
        if i + getCRC(i) >= r
            mUCI = r-getCRC(r,1)-1;
            if mUCI > 0
                warning('For the specified configuration, the number of %s bits (%d) is reduced to %d bits to transmit on PUSCH %d',...
                    t,i,mUCI,n);
                v = d.getPacket(mUCI);
            else
                warning('The %s is not transmitted on the specified PUSCH %d configuration',t,n);
                v = [];
            end
        else
            v = d.getPacket(i);
        end
    end
end

% Plot the SCS carriers
function plotCarriers(waveconfig,gridset)

    if (isfield(waveconfig,'ChannelBandwidth'))
        cbw = waveconfig.ChannelBandwidth;
    else
        cbw = [];
    end

    hold on;
    ylimits = [-numel(gridset) numel(gridset)*2.5];
    xlabel('Frequency (MHz)');

    carriers = waveconfig.Carriers;
    [~,ix] = sort([carriers.SubcarrierSpacing]);
    carriers = carriers(ix);

    scsString = mat2str([carriers.SubcarrierSpacing]);
    nrbString = mat2str([carriers.NRB]);
    if (~isempty(cbw))
        titleStr = [num2str(cbw) 'MHz channel, '];
    else
        titleStr = [];
    end
    title([titleStr ' NRB=' nrbString ', SCS=' scsString 'kHz']);

    blue = [0 0.447 0.741];
    red = [0.850 0.325 0.098];
    orange = [0.929 0.694 0.125];
    green = [0.466 0.674 0.188];

    bwpinfo = [gridset.Info];
    bwpscs = [bwpinfo.SubcarrierSpacing];

    xlim_upper = 0;
    
    for i = 1:numel(carriers)

        NRB = carriers(i).NRB;
        SCS = carriers(i).SubcarrierSpacing;
        RBStart = carriers(i).RBStart;

        bwpidx = find(SCS == bwpscs,1);
        k0 = bwpinfo(bwpidx).k0;

        if (~isempty(cbw) && isfield(waveconfig,'FrequencyRange'))
            guardband = getGuardband(waveconfig.FrequencyRange,cbw,SCS);
        else
            guardband = [];
        end

        SCS = SCS / 1e3; % SCS now in MHz

        ypos = numel(carriers) - i;

        for rb = 0:(NRB-1)

            f = ((rb*12) + k0 - (NRB*12/2)) * SCS;
            r = rectangle('Position',[f ypos 12*SCS 1]);
            r.FaceColor = orange;

            if (rb==0)
                point_a = f - (RBStart*12*SCS);
            end
            if (rb==(NRB-1))
                xlim_upper = max(xlim_upper,f + 12*SCS);
            end

        end

        if (~isempty(guardband))
            p_guardband = plot(ones(1,2) * (-cbw/2 + guardband), ypos + [0 1],'Color',red,'LineWidth',2);
            plot(ones(1,2) * (cbw/2 - guardband), ypos + [0 1],'Color',red,'LineWidth',2);
        end
        
        p_point_a = plot(ones(1,2) * point_a, ypos + [-0.2 1.2],'-..','Color',green,'LineWidth',2);
        
        p_k0 = plot(ones(1,2) * k0 * SCS,ypos + [0.1 0.9],'Color',blue,'LineWidth',2);

        if (i==numel(carriers))
            p_f0 = plot(ones(1,2) * k0 * SCS,ylimits,'k:');
        end

    end

    xlim_lower = point_a;
    if (~isempty(cbw))
        p_channel = plot([-cbw/2 -cbw/2],ylimits,'k--');
        plot([cbw/2 cbw/2],ylimits,'k--');
        xlim_lower = min(xlim_lower,-cbw/2);
        xlim_upper = max(xlim_upper,cbw/2);
    end
    ax = p_point_a.Parent;
    span = xlim_upper - xlim_lower;
    ax.XLim = [xlim_lower-(span*0.05) xlim_upper+(span*0.05)];
    ax.YLim = ylimits;
    ax.YTick = [];

    p = [];
    legends = {};
    if (~isempty(cbw))
        p = [p p_channel];
        legends = [legends 'Channel edges'];
    end
    if (~isempty(guardband))
        p = [p p_guardband];
        legends = [legends 'Guardband edges'];
    end
    p = [p p_point_a p_k0 p_f0];
    legends = [legends 'Point A' 'k_0' 'f_0'];
    legend(p,legends);

end

% Section 5.3.3, TS 38.101-1 (FR1), TS 38.101-2 (FR2)
% Minimum guardband and transmission bandwidth configuration
function guardband = getGuardband(fr,cbw,scs)

    % Table 5.3.3-1: Minimum guardband [kHz] (FR1) (TS 38.101-1)
    cbwFR1    = [    5     10    15     20     25     30     40     50     60     80     90    100];
    guardsFR1 = [242.5  312.5 382.5  452.5  522.5  592.5  552.5  692.5    NaN    NaN    NaN    NaN; ...
                 505.0  665.0 645.0  805.0  785.0  945.0  905.0 1045.0  825.0  925.0  885.0  845.0; ...
                   NaN 1010.0 990.0 1330.0 1310.0 1290.0 1610.0 1570.0 1530.0 1450.0 1410.0 1370.0];
    scsFR1 = [15 30 60].';

    % Table 5.3.3-1: Minimum guardband [kHz] (FR2) (TS 38.101-2)
    cbwFR2    = [  50  100  200  400];
    guardsFR2 = [1210 2450 4930  NaN; ...
                 1900 2420 4900 9860];
    scsFR2 = [60 120].';

    % Return value in MHz
    if (strcmpi(fr,'FR1'))
        guardband = guardsFR1(scsFR1==scs,cbwFR1==cbw) / 1e3;
    else % FR2
        guardband = guardsFR2(scsFR2==scs,cbwFR2==cbw) / 1e3;
    end

end

function carrier = getCarrierObject(ue)
%getCarrierObject Provides the carrier configuration object
%
%   CARRIER = getCarrierObject(UE) provides the carrier configuration
%   object CARRIER , given the input structure UE-wide settings UE.

    % Get the carrier configuration with the input ue
    carrier = nrCarrierConfig;
    carrier.SubcarrierSpacing = ue.SubcarrierSpacing;
    carrier.CyclicPrefix = ue.CyclicPrefix;
    carrier.NSizeGrid = ue.NRB;
end

function srs = getSRSObject(chs)
%getSRSObject Provides the sounding reference signal object
%
%   SRS = getSRSObject(CHS) provides the sounding reference signal
%   configuration object SRS, given the input structure SRS-specific
%   configuration CHS.

    % Get the srs configuration object with the chs input
    srs = nrSRSConfig;
    % Common field names to both chs structure and srs object
    fieldnames = {'NumSRSPorts','NumSRSSymbols','SymbolStart','KTC',...
        'KBarTC','CyclicShift','NRRC','CSRS','BSRS','BHop','Repetition',...
        'GroupSeqHopping','NSRSID'};
    % Assign the values of common field names to respective properties
    for i = 1:numel(fieldnames)
        srs.(fieldnames{i}) = chs.(fieldnames{i});
    end
    % Assign the value of FreqStart field to FrequencyStart property
    srs.FrequencyStart = chs.FreqStart;

    % When AllocatedPeriod is empty or the number of allocated slots per
    % period is more than one, set the resource type to 'aperiodic'
    % (transmission triggered by DCI). For single-slot allocation per
    % period, set resource type to 'periodic'
    numAllocatedSlots = length(chs.AllocatedSlots);
    if isempty(chs.AllocatedPeriod) || numAllocatedSlots > 1
        srs.ResourceType = 'aperiodic';
        srs.SRSPeriod = [1 0];
    elseif numAllocatedSlots == 1
        srs.ResourceType = 'periodic';
        srs.SRSPeriod = [chs.AllocatedPeriod(1) chs.AllocatedSlots];
    end

end

function pusch = getPUSCHObject(chs)
%getPUSCHObject Provides the PUSCH configuration object
%
%   PUSCH = getPUSCHObject(CHS) provides the physical uplink shared channel
%   object PUSCH, given the input structure channel specific transmission
%   configuration CHS.

    % Get the pusch configuration object with the chs input
    pusch = nrPUSCHConfig;
    dmrs = nrPUSCHDMRSConfig;
    dmrs = passign(chs,dmrs,'PortSet','DMRSPortSet');
    numDMRSPorts = numel(dmrs.DMRSPortSet);
    if numDMRSPorts
        % Assign the number of layers to the number of DM-RS antenna ports
        pusch.NumLayers = numDMRSPorts;
    else
        % Get the number of layers, when DM-RS antenna port set is empty
        pusch = passign(chs,pusch,'NLayers','NumLayers');
    end
    pusch.Modulation = chs.Modulation;
    pusch = passign(chs,pusch,'PUSCHMappingType','MappingType');
    % Get SymbolAllocation value, depending on the values of SymbolSet
    % field
    if ~isfield(chs,'SymbolSet')
        pusch.SymbolAllocation = [0 carrier.SymbolsPerSlot];
    elseif isempty(chs.SymbolSet)
        pusch.SymbolAllocation = [0 0];
    else
        [lb,ub] = bounds(chs.SymbolSet);
        pusch.SymbolAllocation = [lb ub-lb+1];
    end
    pusch.PRBSet = chs.PRBSet;
    pusch = passign(chs,pusch,'TransformPrecoding');
    pusch = passign(chs,pusch,'TxScheme','TransmissionScheme');
    pusch = passign(chs,pusch,'TPMI');
    pusch = passign(chs,pusch,'NAntennaPorts','NumAntennaPorts');
    % Get FrequencyHopping value, depending on the values of both
    % InterSlotFreqHopping and IntraSlotFreqHopping fields
    if isfield(chs,'InterSlotFreqHopping') && strcmpi(chs.InterSlotFreqHopping,'enabled')
        pusch.FrequencyHopping = 'interSlot';
    elseif isfield(chs,'IntraSlotFreqHopping') && strcmpi(chs.IntraSlotFreqHopping,'enabled')
        pusch.FrequencyHopping = 'intraSlot';
    end
    pusch = passign(chs,pusch,'RBOffset','SecondHopStartPRB');
    pusch = passign(chs,pusch,'RNTI');

    % Set DM-RS parameters
    dmrs = passign(chs,dmrs,'DMRSConfigurationType','DMRSConfigurationType',~pusch.TransformPrecoding);
    dmrs = passign(chs,dmrs,'DMRSTypeAPosition');
    dmrs = passign(chs,dmrs,'DMRSAdditionalPosition');
    dmrs = passign(chs,dmrs,'DMRSLength');
    dmrs = passign(chs,dmrs,'DMRSSymbolSet','CustomSymbolSet');
    dmrs = passign(chs,dmrs,'NIDNSCID');
    dmrs = passign(chs,dmrs,'NSCID','NSCID',~pusch.TransformPrecoding);
    % Get GroupHopping and SequenceHopping property values, depending on
    % GroupHopping field
    if isfield(chs,'GroupHopping') && pusch.TransformPrecoding
        if strcmpi(chs.GroupHopping,'enable')
            dmrs.GroupHopping = 1;
        elseif strcmpi(chs.GroupHopping,'disable')
            dmrs.SequenceHopping = 1;
        end
    end
    dmrs = passign(chs,dmrs,'NRSID');
    % Get the value of NumCDMGroupsWithoutData property depending on the
    % value of field NumCDMGroupsWithoutData
    if isfield(chs,'NumCDMGroupsWithoutData') && ~isempty(chs.NumCDMGroupsWithoutData) ...
            && ~pusch.TransformPrecoding && chs.NumCDMGroupsWithoutData
        dmrs.NumCDMGroupsWithoutData = chs.NumCDMGroupsWithoutData;
    else
        % When the field NumCDMGroupsWithoutData is not present or is set
        % to empty or zero, assign the value to 1 or 2, depending on
        % transform precoding disabled or enabled
        dmrs.NumCDMGroupsWithoutData = 1 + pusch.TransformPrecoding;
    end
    pusch.DMRS = dmrs;

    % Set PT-RS parameters
    pusch = passign(chs,pusch,'EnablePTRS');
    ptrs = nrPUSCHPTRSConfig;
    ptrs = passign(chs,ptrs,'PTRSTimeDensity','TimeDensity');
    ptrs = passign(chs,ptrs,'PTRSFrequencyDensity','FrequencyDensity',~pusch.TransformPrecoding);
    ptrs = passign(chs,ptrs,'PTRSNumSamples','NumPTRSSamples',pusch.TransformPrecoding);
    ptrs = passign(chs,ptrs,'PTRSNumGroups','NumPTRSGroups',pusch.TransformPrecoding);
    ptrs = passign(chs,ptrs,'PTRSREOffset','REOffset');
    ptrs = passign(chs,ptrs,'PTRSPortSet');
    [ptrs,fieldPresent] = passign(chs,ptrs,'PTRSNID','NID',pusch.TransformPrecoding);
    if ~fieldPresent
        % Assign NID with DM-RS scrambling identity NRSID, when PTRSNID
        % field is not present or is empty
        ptrs.NID = pusch.DMRS.NRSID;
    end
    pusch.PTRS = ptrs;

    % Disable PT-RS
    if pusch.TransformPrecoding
        cond1 = isfield(chs,'PTRSNumSamples') && isempty(chs.PTRSNumSamples);
        cond2 = isfield(chs,'PTRSNumGroups') && isempty(chs.PTRSNumGroups);
        cond = cond1 | cond2;
    else
        cond = isfield(chs,'PTRSFrequencyDensity') && isempty(chs.PTRSFrequencyDensity);
    end
    if cond || (isfield(chs,'PTRSTimeDensity') && isempty(chs.PTRSTimeDensity))
        pusch.EnablePTRS = 0;
    end

    % Set the UCI related parameters
    pusch = passign(chs,pusch,'BetaOffsetACK');
    pusch = passign(chs,pusch,'BetaOffsetCSI1');
    pusch = passign(chs,pusch,'BetaOffsetCSI2');
    pusch = passign(chs,pusch,'ScalingFactor','UCIScaling');

end

function [o,cond] = passign(s,o,f,p,ac)

    cond = isfield(s,f) && ~isempty(s.(f));
    if nargin == 5
        cond = cond && ac;
    end

    if cond
        if nargin == 3
            o.(f) = s.(f);
        else
            o.(p) = s.(f);
        end
    end
end