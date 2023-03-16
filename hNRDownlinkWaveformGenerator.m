%hNRDownlinkWaveformGenerator Generate a 5G NR downlink carrier waveform
%   [WAVEFORM,GRIDSET,INFO] = hNRDownlinkWaveformGenerator(WAVECONFIG) generates
%   a 5G NR downlink waveform WAVEFORM given input WAVECONFIG parameters.
% 
%   WAVECONFIG is a structure containing the following fields:
%   NCellID               - Physical layer cell identity (0...1007)
%   ChannelBandwidth      - Channel bandwidth in MHz
%   FrequencyRange        - Frequency range ('FR1','FR2'). The frequency range
%                           and the channel bandwidth must be a valid ones as
%                           per TS 38.101-1 for FR1 and 38.101-2 for FR2
%   NumSubframes          - Number of 1ms subframes in generated waveform
%   Windowing             - Number of time-domain samples over which 
%                           windowing and overlapping of OFDM symbols is 
%                           applied
%   DisplayGrids          - Display the grids after signal generation (0,1)
%   Carriers - Structure array of SCS-specific carrier configurations, with fields:
%       SubcarrierSpacing - Subcarrier spacing configuration in kHz
%                           Standard configurations are 15,30,60,120,240
%                           for normal CP, and 60 for extended CP
%       NRB               - Number of resource blocks
%       RBStart           - Starting CRB index of SCS carrier, relative to 'Point A'
%   BWP - Structure array of bandwidth part configurations, with fields:
%       SubcarrierSpacing - Subcarrier spacing configuration in kHz. Note
%                           that a SCS carrier with the same SCS must be defined
%       CyclicPrefix      - Cyclic prefix ('normal','extended')
%       NRB               - Number of resource blocks in bandwidth part
%       RBOffset          - Starting index of BWP in the SCS carrier
%   SSBurst - Structure containing SS burst configuration, with fields:
%       Enable            - Enable SS Burst (0,1)
%       BlockPattern      - Burst type ('Case A','Case B','Case C','Case D','Case E')
%       SSBTransmitted    - Bitmap indicating blocks transmitted in a 5ms half-frame burst
%       SSBPeriodicity    - SS burst set periodicity in ms (5,10,20,40,80,160)
%       FrequencySSB      - Frequency offset of SS burst (Hz), relative to waveform center (multiples of 5kHz)
%       Power             - Power scaling (dB)
%   CSIRS - Structure array of CSI-RS configurations, with fields:
%       Enable            - Enable the set of CSI-RS resources
%       BWP               - Bandwidth part identity. It is indicated to
%                           which bandwidth part the set of CSI-RS
%                           resources are configured
%       Power             - Power scaling of CSI-RS resource(s) in dB
%       CSIRSType         - Type of CSI-RS resource(s) ('nzp','zp')
%       RowNumber         - Row number corresponds to CSI-RS resource(s)
%                           as defined in TS 38.211 Table 7.4.1.5.3-1 (1...18)
%       Density           - Frequency density of CSI-RS resource(s) ('one','three','dot5even','dot5odd')
%       SubcarrierLocations - Subcarrier locations of CSI-RS resource(s)
%                           within a resource block (RB)
%       NumRB             - Number of RBs allocated to CSI-RS resource(s) (1...275)
%       RBOffset          - Starting RB index of CSI-RS resource(s)
%                           allocation relative to carrier resource grid (0...2169)
%       SymbolLocations   - OFDM symbol locations of CSI-RS resource(s) within a slot
%       AllocatedSlots    - Slots allocated (0-based) to CSI-RS resource(s)
%                           within a period
%       AllocatedPeriod   - Periodicity of CSI-RS resource(s) allocation in slots
%       NID               - Scrambling identity corresponds to CSI-RS
%                           resource(s) (0...1023)
%   CORESET - Structure array of CORESET/search space sequence configurations, with fields:
%       AllocatedSymbols  - First symbol index to each CORESET opportunity in a slot
%       AllocatedSlots    - Slot indices of each CORESET opportunity in the sequence period
%       AllocatedPeriod   - Period of sequence, in terms of slots. Use empty for no repetition                          
%       Duration          - CORESET symbol duration (1,2,3)
%       AllocatedPRB      - PRB indices associated with REG in CORESET (in mapped BWP)
%       CCEREGMapping     - CCE-to-REG mapping, 'interleaved' or 'noninterleaved'
%       REGBundleSize     - REG bundle size, L, (2,6) or (3,6)
%       InterleaverSize   - Interleaver size, R, (2,3,6)
%       ShiftIndex        - Shift index, 0...274
%   PDCCH - Structure array of PDCCH sequence configurations, with fields:
%       Enable            - Enable PDCCH sequence
%       BWP               - Bandwidth part ID (1-based) for with PDCCH sequence.
%                           Note that the associated BWP configuration must 
%                           be defined
%       Power             - Power scaling (dB)
%       EnableCoding      - Enable DCI coding (0,1)
%       CORESET           - CORESET configuration ID associated this PDCCH sequence
%       AllocatedSearchSpaces - Indices of each CORESET opportunity, carrying PDCCH in the sequence period
%       AllocatedPeriod   - Period of sequence, in terms of slots. Use empty for no repetition          
%       NumCCE            - Number of CCE used by each PDCCH instance
%       StartCCE          - Starting CCE of PDCCH
%       RNTI              - RNTI
%       NID               - PDCCH and DM-RS scrambling NID
%       PowerDMRS         - Additional power boosting (dB)
%       DataBlkSize       - DCI payload size
%       DataSource        - Data source. Use a numerical array or one of 
%                           the following standard PN sequences: 'PN9-ITU',
%                           'PN9','PN11', 'PN15', 'PN23'. The generator seed for the
%                           can be specified using a cell array of the
%                           form {'PN9',seed}, otherwise the seed is all ones
%   PDSCH - Structure array of PDSCH sequence configurations, with fields:
%       Enable            - Enable the PDSCH sequence configuration (0,1)
%       BWP               - Bandwidth part ID (1-based) for with PDSCH sequence.
%                           Note that the associated BWP configuration must 
%                           be defined
%       Power             - Power scaling (dB)
%       EnableCoding      - Enable DL-SCH transport channel encoding (0,1)
%       DataSource        - Data source. Use a numerical array or one of 
%                           the following standard PN sequences: 'PN9-ITU',
%                           'PN9','PN11', 'PN15', 'PN23'. The generator seed for the
%                           can be specified using a cell array of the
%                           form {'PN9',seed}, otherwise the seed is all ones
%       TargetCodeRate    - Code rate used to calculate transport block sizes
%       Xoh_PDSCH         - Rate matching overhead parameter (nominally 0,6,12,18)
%       Modulation        - Modulation scheme ('QPSK','16QAM','64QAM','256QAM')
%       NLayers           - Number of layers (1...8)
%       RVSequence        - Redundancy version sequence, applied cyclically 
%                           all PDSCH instances in sequence
%       VRBToPRBInterleaving - Enable or disable the interleaving of the
%                           virtual to physical resource block mapping (0,1)
%       VRBBundleSize     - Bundle size for the interleaving, specified
%                           through the higher-layer parameter vrb-ToPRB-Interleaver (2,4)
%       AllocatedSymbols  - Symbol indices allocated to a PDSCH instance in a slot
%       AllocatedSlots    - Slot indices of each PDSCH instance in the sequence period
%       AllocatedPeriod   - Period of sequence, in terms of slots. Use empty for no repetition                          
%       AllocatedPRB      - PRB allocation
%       RNTI              - Radio network temporary identifier
%       NID               - Scrambling identity for data part
%       RateMatch         - Structure defining reserved PRB rate matching patterns 
%       PortSet           - DM-RS antenna ports associated with layers
%       DMRSTypeAPosition - DM-RS first symbol position for mapping type 'A' (2,3)
%       DMRSLength        - Number of front-loaded DM-RS symbols.
%                           Either 1 (single symbol) or 2 (double symbol)
%       DMRSAdditionalPosition  - Additional DM-RS symbols positions (0...3)
%       DMRSConfigurationType   - DM-RS configuration type (1,2)
%       NumCDMGroupsWithoutData - Number of DM-RS CDM groups without data (0...3)
%       NIDNSCID          - Scrambling identity for DM-RS
%       NSCID             - DM-RS scrambling initialization (0,1)
%       PowerDMRS         - Power offset for DM-RS (dB)
%       EnablePTRS        - Enable the PT-RS configuration (0,1)
%       PTRSFrequencyDensity - Frequency density of PT-RS (2,4)
%       PTRSTimeDensity   - Time density of PT-RS (1,2,4)
%       PTRSREOffset      - Resource element offset ('00','01','10','11')
%       PTRSPortSet       - Antenna ports of PT-RS
%       PowerPTRS         - Power boosting for PT-RS (dB)

%   Copyright 2018-2020 The MathWorks, Inc.

function [waveform,gridset,winfo] = hNRDownlinkWaveformGenerator(waveconfig)

    % Unbundle the channel specific parameter structures for easier access
    ssburst = waveconfig.SSBurst;
    carriers = waveconfig.Carriers;
    bwp = waveconfig.BWP;
    coreset = waveconfig.CORESET;
    pdcch = waveconfig.PDCCH;
    pdsch = waveconfig.PDSCH;

    % Defaulting for the grid plotting
    if ~isfield(waveconfig,'DisplayGrids')
        waveconfig.DisplayGrids = 0;
    end
    
    % Remove WindowCoeffs field if present, as the value of this field is
    % not necessarily valid for multiple carriers
    if (isfield(waveconfig,'WindowCoeffs'))
        waveconfig = rmfield(waveconfig,'WindowCoeffs');
    end
    
    % Define the instrumentation info variable 
    waveinfo = struct('PDCCH',[],'PDSCH',[]);
    
    % Cross-check the BWP and SCS carrier configs
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

    % Create BWP *PRB* resource grids
    ResourceGrids = arrayfun(@(bp)zeros(bp.NRB,waveconfig.NumSubframes*1*symbolsPerSlot(bp)*fix(bp.SubcarrierSpacing/15)),...
                                    bwp,'UniformOutput',false);

    % Create BWP subcarrier resource grids
    % Size ALL BWP RE grids by the max number of layers/ports in the enabled PDSCH
    maxlayers = max([pdsch(logical([pdsch(:).Enable])).NLayers]);
    if isempty(maxlayers) || maxlayers <= 0
        maxlayers = 1; % Ensure that there is at least one plane
    end
    numPorts = ones(1,numel(bwp))*maxlayers;

    % Create BWP *resource element grids*
    ResourceElementGrids = arrayfun(@(bp)zeros(bp.NRB*12,waveconfig.NumSubframes*1*symbolsPerSlot(bp)*fix(bp.SubcarrierSpacing/15),max(numPorts)),...
        bwp,'UniformOutput',false);

    % Process the set of CSI-RS transmissions
    if isfield(waveconfig,'CSIRS')
        csirs = waveconfig.CSIRS;
        % Unbundle CSI-RS specific parameter structure
        for nsig = 1:length(csirs)
            % Only process configuration if enabled
            if ~csirs(nsig).Enable
                continue;
            end

            % Check the referenced BWP index
            checkIndex('CSI-RS',csirs,nsig,'BWP',bwp);

            % Recreate carrier specific configuration objects using carrier
            % specific parameters
            carrierIdx = bwp(csirs(nsig).BWP).CarrierIdx;
            tempCarrier = nrCarrierConfig;
            tempCarrier.NSizeGrid = carriers(carrierIdx).NRB;
            tempCarrier.NStartGrid = carriers(carrierIdx).RBStart;
            tempCarrier.SubcarrierSpacing = carriers(carrierIdx).SubcarrierSpacing;
            tempCarrier.CyclicPrefix = bwp(csirs(nsig).BWP).CyclicPrefix;
            carrierCfgs(nsig) = tempCarrier; %#ok<AGROW>

            % Recreate CSI-RS specific configuration objects using the
            % parameters provided by CSIRS field of waveconfig
            tempCSIRS = nrCSIRSConfig;
            tempCSIRS.CSIRSType = csirs(nsig).CSIRSType;
            tempCSIRS.RowNumber = csirs(nsig).RowNumber;
            tempCSIRS.Density = csirs(nsig).Density;
            tempCSIRS.SymbolLocations = csirs(nsig).SymbolLocations;
            tempCSIRS.SubcarrierLocations = csirs(nsig).SubcarrierLocations;
            tempCSIRS.NumRB = csirs(nsig).NumRB;
            tempCSIRS.RBOffset = csirs(nsig).RBOffset;
            tempCSIRS.NID = csirs(nsig).NID;
            csirsCfgs(nsig) = tempCSIRS; %#ok<AGROW>
        end

        % Update the number of antenna ports corresponding to each BWP
        % based on CSI-RS configuration
        for bwpIdx = 1:numel(bwp)
            csirsBWPCheck = ([csirs.BWP] == bwpIdx);
            if any(csirsBWPCheck) && any([csirs(csirsBWPCheck).Enable]) % If a CSI-RS instance is configured for this BWP and if it is enabled
                tempCfg = csirsCfgs(csirsBWPCheck);
                numCSIRSPorts = tempCfg(logical([csirs(csirsBWPCheck).Enable])).NumCSIRSPorts;
                numPorts(bwpIdx) = max([numCSIRSPorts maxlayers]);
            end
        end

        % Update BWP resource element grids based on CSI-RS ports
        ResourceElementGrids = arrayfun(@(bp)zeros(bp.NRB*12,waveconfig.NumSubframes*1*symbolsPerSlot(bp)*fix(bp.SubcarrierSpacing/15),max(numPorts)),...
            bwp,'UniformOutput',false);

        numSubframes = waveconfig.NumSubframes;
        for nch = 1:length(pdsch)
            % Only process configuration if enabled
            if ~pdsch(nch).Enable
                continue;
            end

            % Check the referenced BWP index
            checkIndex('PDSCH',pdsch,nch,'BWP',bwp);

            scs = bwp(pdsch(nch).BWP).SubcarrierSpacing;
            numSlotsPerSubframe = scs/15;
            numSlots = numSlotsPerSubframe*numSubframes;
            pdsch(nch).CSIRSREs = cell(1,numSlots);
        end

        % Extract the power scalings of channel state information reference
        % signals
        PowerCSIRS = {csirs.Power};

        for nsig = 1:length(csirs)
            sig = csirs(nsig);
            % Only process configuration if enabled
            if ~sig.Enable
                continue;
            end

            carrierCfg = carrierCfgs(nsig);
            csirsCfg = csirsCfgs(nsig);
            % Validate CSI-RS resource allocation parameters to check
            % whether the CSI-RS resource is confined to the respective
            % BWP or not
            validateCSIRSAllocationParameters(csirsCfg,nsig,bwp(csirs(nsig).BWP),sig.BWP);

            % Get the number of CSI-RS resources configured
            if iscell(csirsCfg.CSIRSType)
                % Single or multiple resources
                numRes = numel(csirsCfg.CSIRSType);
            else
                % Single resource
                numRes = 1;
            end

            % Validate the number of values provided for CSI-RS power
            % scaling
            powerCSIRS = PowerCSIRS{nsig}(:);
            powerCSIRS = applyScalarExpansion(powerCSIRS,numRes,'Power',nsig);
            powerCSIRS = db2mag(powerCSIRS); % dB to magnitude conversion

            % Validate the number of values provided for AllocatedSlots
            allocSlots = sig.AllocatedSlots;
            if ~iscell(allocSlots)
                allocSlots = {allocSlots};
            end
            allocSlots = applyScalarExpansion(allocSlots,numRes,'AllocatedSlots',nsig);

            % Validate the number of values provided for AllocatedPeriod
            allocPeriod = sig.AllocatedPeriod;
            if ~iscell(allocPeriod)
                allocPeriod = {allocPeriod};
            end
            allocPeriod = applyScalarExpansion(allocPeriod,numRes,'AllocatedPeriod',nsig);

            nrb = bwp(csirs(nsig).BWP).NRB;
            symbperslot = carrierCfg.SymbolsPerSlot;
            bwpGridSize = [12*nrb symbperslot numPorts(sig.BWP)];

            for resIdx = 1:numRes
                % Establish slot allocation for CSI-RS resources
                allocatedSlotsCSIRS{resIdx} = expandbyperiod(allocSlots{resIdx},allocPeriod{resIdx},waveconfig.NumSubframes,carrierCfg.SubcarrierSpacing); %#ok<AGROW>
            end

            numSlots = carrierCfg.SlotsPerSubframe*numSubframes;
            % Loop over 0 to numSlots-1
            for slotIdx = 0:numSlots-1
                carrierCfg.NSlot = slotIdx;
                csirsPeriod = {};
                for resIdx = 1:numRes
                    if any(allocatedSlotsCSIRS{resIdx} == slotIdx)
                        csirsPeriod{resIdx} = 'on'; %#ok<AGROW>
                    else
                        csirsPeriod{resIdx} = 'off'; %#ok<AGROW>
                    end
                end
                csirsCfg.CSIRSPeriod = csirsPeriod;

                % Generate 1-based carrier oriented CSI-RS indices in
                % subscript form
                [csirsInd,~] = nrCSIRSIndices(carrierCfg,csirsCfg,'IndexStyle','subscript');

                % Generate CSI-RS symbols and apply power boosting
                [sym,info] = nrCSIRS(carrierCfg,csirsCfg,'OutputResourceFormat','cell');
                if slotIdx == 0 % This reordering should be done only once
                    powerCSIRS = powerCSIRS(info.ResourceOrder);
                end
                csirsSym = arrayfun(@(x) sym{x}*powerCSIRS(x),1:length(powerCSIRS),'UniformOutput',false);
                csirsSym = cell2mat(csirsSym(:));

                % Change the orientation of CSI-RS indices to BWP
                bwpRBOffset = bwp(sig.BWP).RBOffset;
                bwpCSIRSInd = [csirsInd(:,1)-bwpRBOffset*12 csirsInd(:,2) csirsInd(:,3)];
                % Linearize BWP oriented CSI-RS indices
                csirsLinInd = sub2ind(bwpGridSize,bwpCSIRSInd(:,1),bwpCSIRSInd(:,2),bwpCSIRSInd(:,3));

                % Create an empty slot grid
                slotgrid = zeros(bwpGridSize);

                % Write the CSI-RS symbols in the slot grid
                slotgrid(csirsLinInd) = csirsSym;

                % Combine CSI-RS instance with the rest of the BWP grid
                ResourceElementGrids{sig.BWP}(:,slotIdx*symbperslot+(1:symbperslot),1:bwpGridSize(3)) = ResourceElementGrids{sig.BWP}(:,slotIdx*symbperslot+(1:symbperslot),1:bwpGridSize(3)) + slotgrid;

                for nch = 1:length(pdsch)
                    dch = pdsch(nch);
                    if dch.Enable
                        if (dch.BWP == sig.BWP)
                            % 0-based CSI-RS REs
                            pdsch(nch).CSIRSREs{1,slotIdx+1} = [pdsch(nch).CSIRSREs{1,slotIdx+1}; csirsLinInd - 1];
                        end
                    end
                end
            end
        end
    end

    % Add an empty 'internal' placeholder field to the PDSCH configs to contain
    % any inter-channel reserved PRB which will be created locally below
    [pdsch.Reserved] = deal(struct('Name',{},'PRB',{},'Symbols',{},'Period',{}));

    % Define channel plotting ID markers
    chplevel.PDCCH = 1.3;
    chplevel.PDSCH = 0.8;
    chplevel.CORESET = 1.1;
    chplevel.SS_Burst = 0.6;

    % Mark the CORESET sequences for display in one of the BWP grids, for visualization purposes only
    for cs = coreset
        nbwp = 1;
        symbperslot = symbolsPerSlot(bwp(nbwp)); 

        % Establish slot allocation for CORESET 
        allocatedSlots = expandbyperiod(cs.AllocatedSlots,cs.AllocatedPeriod,waveconfig.NumSubframes,bwp(1).SubcarrierSpacing);    % Slots containing CORESET instances 

        % Establish PRB allocation for CORESET
        allocatedPRB = getCORESETPRB(carriers(bwp(nbwp).CarrierIdx),bwp(nbwp),cs);

        % Establish symbol allocation for CORESET 
        allocatedSymbols = cs.AllocatedSymbols(cs.AllocatedSymbols+cs.Duration <= symbperslot);

        for ns = allocatedSlots
            % Expand the PRB into REG
            ResourceGrids{nbwp}(1+allocatedPRB,1+symbperslot*ns+(expander(allocatedSymbols,cs.Duration))) = chplevel.CORESET;
        end
    end

    % Process the set of PDCCH transmission sequences
    for nch = 1:length(pdcch)

        % Get a copy of the current PDCCH channel parameters
        ch = pdcch(nch);

        % Only process configuration if enabled
        if ~ch.Enable 
            continue;
        end
        
        % Establish whether transport coding is enabled
        dcicoding = ~isfield(ch,'EnableCoding') || ch.EnableCoding;
        
        % Check the referenced BWP and CORESET indices
        checkIndex('PDCCH',pdcch,nch,'BWP',bwp,'CORESET',coreset);

        % Get a copy of the current CORESET for this PDCCH sequence
        cset = coreset(ch.CORESET);
        cset = validateCoreset(cset);

        % Get the number of symbols per slot for the associated BWP (CP dependent)
        symbperslot = symbolsPerSlot(bwp(ch.BWP));

        % Exclude any CORESET occasions that would fall outside a slot
        slotsymbs = cset.AllocatedSymbols(cset.AllocatedSymbols+cset.Duration <= symbperslot);
        if length(slotsymbs) ~= length(cset.AllocatedSymbols)
            warning('CORESET %d (%d symbol duration) in BWP %d includes positions which fall outside the slot (0...%d). Using only CORESET locations within a slot.',ch.CORESET,cset.Duration,ch.BWP,symbperslot-1);
            cset.AllocatedSymbols = slotsymbs;
        end

        % Calculate the initial symbol and slot numbers for the CORESET/search space
        % monitoring occasions by expanding by the period across the waveform length
        potentialslots = expandbyperiod(cset.AllocatedSlots,cset.AllocatedPeriod,waveconfig.NumSubframes,bwp(ch.BWP).SubcarrierSpacing);
        potentialsymbols = reshape(symbperslot*potentialslots + cset.AllocatedSymbols',1,[]);

        % Also need to expand the indices of the allocated monitoring locations
        % that will be used by the current PDCCH sequence
        % Expand by period so that it covers all the potential symbols 
        allocslotindices = expandbyperiod(ch.AllocatedSearchSpaces,ch.AllocatedPeriod,numel(potentialsymbols));    

        % Identify the absolute initial symbols associated with the CORESET instances carry the PDCCH 
        allocatedsymbols = potentialsymbols(1+allocslotindices);

        % Identify the locations of the PDCCH instances
        % - For displaying 
        % - For actual PDCCH instances mapping, to expand into resource elements (eventually data & DM-RS)
        % - For reserving in the associated PDSCH (connected by RNTI)
        % 
        % Initial format is array of (symbol,prb) in mapping order
        
        nrb = bwp(nbwp).NRB;
        if max(cset.AllocatedPRB) >= nrb
            error('The PRB allocation part of CORESET %d (set of 6 PRB/REG CCEs, max 0-based PRB index = %d) exceeds the number of available RB (%d) in BWP %d.',ch.CORESET,max(cs.AllocatedPRB),nrb,ch.BWP);
        end

        % Establish PRB allocation for CORESET
        [potentialPRB,ncces] = getCORESETPRB(carriers(bwp(nbwp).CarrierIdx),bwp(nbwp),cset);

        % Check the PDCCH allocation dimensionality, relative to the CORESET
        maxcce = cset.Duration*numel(ncces); 
        if ch.StartCCE + ch.NumCCE > maxcce    
            error('The CCE range ([%d,%d] inclusive) of PDCCH %d exceeds the maximum number of CCE (%d) in CORESET %d.',ch.StartCCE,ch.StartCCE+ch.NumCCE-1,nch, maxcce, ch.CORESET);
        end

        % Calculate the subscripts of the PDCCH PRB associated with CCE set 
        % Turn CCE into blocks of 6 REG/RB, in terms of symbols & PRB indices
        expcceindices = 6*ch.StartCCE + (0:6*ch.NumCCE-1);     

        % Find the CCE mapping
        if strcmp(cset.CCEREGMapping,'interleaved')
            numREGs = numel(potentialPRB)*cset.Duration;
            L = cset.REGBundleSize;
            R = cset.InterleaverSize;
            C = numREGs/(L*R);
            f = zeros(R*C,1);       % REG Bundles (interleaved)
            for cIdx = 0:C-1
                for rIdx = 0:R-1
                    x = cIdx*R + rIdx;
                    f(x+1) = mod(rIdx*C + cIdx + cset.ShiftIndex, R*C);
                end
            end

            % Slot based indices
            %   Supports only scalar value of AllocatedSymbols    
            if length(cset.AllocatedSymbols)>1
                error('AllocatedSymbols must be a scalar for interleaved CORESETs.');
            end
            rbGrid = zeros(bwp(ch.BWP).NRB,symbperslot);
            symNum = expander(cset.AllocatedSymbols,cset.Duration);
            rbGrid(potentialPRB+1,symNum+1) = 1;
            rbIdx = find(rbGrid);
            rbIdxTime = reshape(reshape(rbIdx,[],cset.Duration).',[],1);

            % Interleave
            regB = reshape(rbIdxTime,L,[]);
            rbIdxInt = regB(:,f+1);
            rbIdxTime = rbIdxInt(:);
            rbPDCCHIdx = rbIdxTime(1:ch.NumCCE*6)-1;

            % Build up the prblocations for downstream sharing
            tmp = sort(reshape(rbPDCCHIdx,cset.Duration,[])')' - ...
               bwp(ch.BWP).NRB.*(0:cset.Duration-1).' - ...
               bwp(ch.BWP).NRB.*(cset.AllocatedSymbols(1)*ones(cset.Duration,1)); %#ok

            prblocations = [mod(expcceindices,cset.Duration)', tmp(:)];
        else    % noninterleaved
            prblocations = [mod(expcceindices,cset.Duration)',potentialPRB(1+fix(expcceindices/cset.Duration))']; % 0-origin symbols
        end
        
        % Turn this into a cell array per symbol for reservation      
        celloc1d = cell(1,cset.Duration);
        for c = 1:cset.Duration
          celloc1d{c} = reshape(prblocations(prblocations(:,1)==(c-1),2),1,[]);  % Store each PRB set as a row
        end  
        % Repeat the 1 slot allocation cell array across all transmission occasions
        celloc = repmat(celloc1d,1,length(allocatedsymbols));

        % Expand symbols to match the CORESET duration
        symloc = expander(allocatedsymbols,cset.Duration);

        % Create a PRB reservation configuration
        reserved.Name = sprintf('Reserved for PDCCH transmission %d',nch);
        reserved.PRB = celloc;     % Reserved PRB (0-based indices, defined as a vector or cell array)
        reserved.Symbols = symloc; % OFDM symbols associated with reserved PRB (0-based indices, spanning one or more slots)
        reserved.Period = [];      % Total number of slots in the pattern period

        % If this PDCCH sequence is associated with any of the PDSCH (same RNTI and BWP)
        % then configure the PDSCH sequence with the associated PDCCH PRB reservation
        rnti = ch.RNTI;
        for i = 1:length(pdsch)
            dch = pdsch(i);
            if (dch.RNTI == rnti) && (dch.BWP == ch.BWP)
                pdsch(i).Reserved(end+1) = reserved;
            end
        end

        % Relative PRB linear indices in BWP (relative to 0-origin CORESET), for a single PDCCH instance 
        indblock = sum(prblocations .* [bwp(ch.BWP).NRB 1], 2);   % Column of 0-based linear indices

        % Mark the PRB used by this PDCCH for display in the associated BWP grid
        for s = allocatedsymbols
            ResourceGrids{ch.BWP}(1+(allocatedsymbols*bwp(ch.BWP).NRB)+indblock) = chplevel.PDCCH;  % Using column expansion to repeat the instances
        end

        % Turn subscripts into 0-based RE indices for a single PDCCH instance
        indpdcch = reSub2Ind(bwp(ch.BWP).NRB,celloc1d,0:cset.Duration-1);
        repdcch = expander(12*indpdcch,12,4,1,true)';  % Step by 4 RE with an offset of 1     
        redmrs = expander(12*indpdcch,12,4,1,false)';  % RE carrying DM-RS

        % PDCCH instance symbol/bit capacity
        Gd = length(repdcch);
        G = Gd*2;

        % Create a data source for this PDCCH sequence
        datasource = hVectorDataSource(ch.DataSource); 

        % Storage for PDCCH instance information
        controlstore = [];
          
        % Loop over all the PDCCH transmission occasions and write the encoded
        % DCI payloads into the resource elements of the associated PDCCH instances
        si = 0;
        for s = allocatedsymbols

            si = si+1;
            nID = ch.NID;   % Scrambling NID value for DM-RS (pdcch-DMRS-ScramblingID or NCellID, depending)

            if dcicoding
                % Get the DCI payload bits from the data source
                dcibits = datasource.getPacket(ch.DataBlkSize);

                % Encode the DCI payload to match the PDCCH bit capacity
                codeword = nrDCIEncode(dcibits,ch.RNTI,G);
            else
                % Get the PDCCH codeword directly from the data source
                codeword = datasource.getPacket(G);
                dcibits = [];
            end
            
            % Get the PDCCH QPSK symbols
            symbols = nrPDCCH(codeword,nID,ch.RNTI);

            % Combine with existing grid (shift indices to initial symbol for the search space instance)
            offset = 1+s*12*bwp(ch.BWP).NRB;
            ResourceElementGrids{ch.BWP}(offset+repdcch) = ResourceElementGrids{ch.BWP}(offset+repdcch) + symbols*db2mag(ch.Power); 

            % Construct and map the PDCCH DM-RS     
            nslot = mod(fix(s/symbperslot), bwp(ch.BWP).SubcarrierSpacing/15 * 10);  % Slot number, in a 10ms frame
            nsym = mod(s,symbperslot); % Symbol number in slot

            % Construct a single symbol column vector of DM-RS for the PDCCH
            dmrssym = [];   
            crboffset = carriers(bwp(ch.BWP).CarrierIdx).RBStart + bwp(ch.BWP).RBOffset; % First CRB associated with start of BWP
            for i = 1:cset.Duration

                prb = celloc1d{i};     % PRB indices associated with current symbol
                slen = max(prb)+1;     % Number of PRB that the continuous PRBS sequence needs to cover 

                % Construct PRBS for the transmission DM-RS, offsetting the 
                % sequence to account for the origin of the BWP
                cinit = mod(2^17*(symbperslot*nslot+nsym+1)*(2*nID+1)+2*nID,2^31);
                nsc = 6;   % 3 DM-RS symbols per RB and 2 PRBS bits needed per QPSK symbol
                cSeq = nrPRBS(cinit,nsc*[crboffset slen]);
                % Select onto the subset required to match the PRB
                cSeq = reshape(cSeq,nsc,[]);
                cSeq = cSeq(:,prb+1);
                % Create associated complex DM-RS symbols
                dmrssym = [dmrssym; nrSymbolModulate(cSeq(:),'QPSK')]; %#ok<AGROW>
                nsym = nsym + 1;  % Increment symbol number 
            end
            % Combine PDCCH with the grid
            ResourceElementGrids{ch.BWP}(offset+redmrs) = ResourceElementGrids{ch.BWP}(offset+redmrs) + dmrssym*db2mag(ch.Power+ch.PowerDMRS);

            % Indices offset for the current slot
            sindoffset = 1+mod(s,symbperslot)*12*bwp(ch.BWP).NRB;

            % Capture resource info for this PDCCH instance
            controlstore(si).NSlot = fix(s/symbperslot); %#ok<AGROW>
            controlstore(si).DCIBits = dcibits;  %#ok<AGROW>
            controlstore(si).Codeword = codeword;  %#ok<AGROW>
            controlstore(si).G = G; %#ok<AGROW>
            controlstore(si).Gd = Gd; %#ok<AGROW>
            controlstore(si).ChannelIndices = sindoffset+repdcch; %#ok<AGROW>
            controlstore(si).ChannelSymbols = symbols*db2mag(ch.Power); %#ok<AGROW>
            controlstore(si).DMRSIndices = sindoffset+redmrs; %#ok<AGROW>
            controlstore(si).DMRSSymbols = dmrssym*db2mag(ch.Power+ch.PowerDMRS); %#ok<AGROW>
            
        end
        
        % Capture all resources info for this PDCCH sequence
        name = "";
        if isfield(pdcch(nch),'Name')
            name = pdcch(nch).Name;
        end
        waveinfo.PDCCH(nch).Name = name;
        waveinfo.PDCCH(nch).CDMLengths = [1 1];
        waveinfo.PDCCH(nch).Resources = controlstore;
        
    % End of PDCCH sequence processing      
    end

    % Get the set of RB level resources, in each SCS carrier, that overlap 
    % with the SS burst, and apply modulation specific settings
    [ssbreserved,ssburst] = ssburstResources(ssburst,carriers);
    ofdmFields = {'WaveformType'; 'Windowing'; 'Alpha'; 'FilterLength'; 'ToneOffset'; 'CyclicExtension'};
    for i = 1:numel(ofdmFields)
        field = ofdmFields{i};
        if (isfield(waveconfig,field))
            ssburst.(field) = waveconfig.(field);
        end
    end

    % Process the set of PDSCH transmission sequences
    % Create a single DL-SCH channel processing object for use with 
    % all the PDSCH sequences
    dlsch = nrDLSCH('MultipleHARQProcesses',false);

    for nch = 1:length(pdsch)

        % Get a copy of the current PDSCH channel parameters
        ch = pdsch(nch);

        % Only process configuration if enabled
        if ~ch.Enable 
            continue;
        end

        % Establish whether transport coding is enabled
        trcoding = ~isfield(ch,'EnableCoding') || ch.EnableCoding;

        % Check the referenced BWP index
        checkIndex('PDSCH',pdsch,nch,'BWP',bwp);

        % Expand the allocated slot sequence by the repetition period, across
        % the length of the waveform
        allocatedSlots = expandbyperiod(ch.AllocatedSlots,ch.AllocatedPeriod,waveconfig.NumSubframes,bwp(ch.BWP).SubcarrierSpacing);
        if any(ch.AllocatedPRB >= bwp(ch.BWP).NRB)
            error('The allocated PRB indices (0-based, largest value = %d) for PDSCH %d exceed the NRB (%d) for BWP %d.',max(ch.AllocatedPRB),nch,bwp(ch.BWP).NRB,ch.BWP);
        end
        % Ensure that the allocated symbols for the slot are within a slot for the BWP CP
        symbperslot = symbolsPerSlot(bwp(ch.BWP));
        slotsymbs = ch.AllocatedSymbols(ch.AllocatedSymbols < symbperslot);
        if length(slotsymbs) ~= length(ch.AllocatedSymbols)
            warning('The slot-wise symbol allocation for PDSCH %d in BWP %d includes 0-based symbol indices which fall outside a slot (0...%d). Using only symbols within a slot.',nch,ch.BWP,symbperslot-1);
            ch.AllocatedSymbols = slotsymbs;
        end
    
        % Reserved PRB-level resources associated with SS burst
        rs = ssbreserved(bwp(ch.BWP).CarrierIdx);
        rs.PRB = rs.PRB - bwp(ch.BWP).RBOffset;
        ch.Reserved(end+1) = rs;         % Configure the channel with the pattern
        
        % Convert rate match pattern configurations into a format suitable 
        % for the hPDSCHResources function
        % 
        % Turn reserved CORESET indices into reserved patterns
        % with format, reserved = struct('Name',{},'PRB',{},'Symbols',{},'Period',{});
        for rmp = ch.RateMatch
            % Process CORESET ratematch pattern part
            if isfield(rmp,'CORESET')         
                for cidx = rmp.CORESET                
                    % Expand and project CORESET into the BWP 
                    % Pattern representation is single vector of PRB across all symbols  

                    % Check the CORESET index
                    if cidx < 1 || cidx > length(coreset)
                        error('For PDSCH %d, the ratematch CORESET index (%d) must be between 1 and the number of CORESET defined (%d)',...
                                         nch,cidx,length(coreset));
                    end

                    % Get a copy of the CORESET configuration
                    cs = coreset(cidx);

                    % Expand the allocated slots across the repetition period
                    rmallocatedSlots = expandbyperiod(cs.AllocatedSlots,cs.AllocatedPeriod,waveconfig.NumSubframes,bwp(ch.BWP).SubcarrierSpacing);
                    % All CORESET symbols in each allocated slot

                    % Expand to identify all symbols included in this CORESET sequence
                    slotsymbs = cs.AllocatedSymbols(cs.AllocatedSymbols+cs.Duration < symbperslot);
                    csetsymbols = expander(slotsymbs,cs.Duration);
                    rmallocatedsymbols = reshape(symbperslot*rmallocatedSlots+csetsymbols',1,[]);  % Use column expansion

                    % Turn the allocated PRB parameter into blocks of 6 RB/REG CCEs
                    allocatedPRB = expander(unique(6*fix(cs.AllocatedPRB/6)),6);

                    % Check that the associated PRB set fits within the associated BWP NRB
                    if max(allocatedPRB) >= bwp(ch.BWP).NRB
                        error('For PDSCH %d, the effective PRB allocation part of RateMatch CORESET %d (set of 6 PRB/REG CCEs, max PRB index = %d) exceeds the number of RB (%d) in BWP %d.',nch,cidx,max(allocatedPRB),bwp(ch.BWP).NRB,ch.BWP);
                    end

                    % Create reserved configuration structure and push it onto the copy of the PDSCH parameters
                    rs.Name = sprintf('Reserved for CORESET %d',cidx);
                    rs.PRB = allocatedPRB;           % Reserved PRB (0-based indices, defined as a vector or cell array)
                    rs.Symbols = rmallocatedsymbols; % OFDM symbols associated with reserved PRB (0-based indices, spanning one or more slots)
                    rs.Period = [];                  % Total number of slots in the pattern period (empty means don't cyclically repeat)
                    ch.Reserved(end+1) = rs;         % Configure the channel with the pattern
                end
            end
            % Process bitmap derived ratematch pattern part
            if isfield(rmp,'Pattern')
                for rmpat = rmp.Pattern
                    % Name this pattern for identification purposes
                    rs.Name = sprintf('Reserved for rate-matching pattern');

                    % Check that the associated PRB set fits within the associated BWP NRB
                    if max(rmpat.AllocatedPRB) >= bwp(ch.BWP).NRB
                        error('For PDSCH %d, the PRB allocation part of the RateMatch pattern exceeds the number of RB (%d) in BWP %d.',nch,bwp(ch.BWP).NRB,ch.BWP);
                    end
                    rs.PRB = rmpat.AllocatedPRB;
                    % Need to combine allocated symbols and allocated slots into a single list             
                    rmallocslots = expandbyperiod(rmpat.AllocatedSlots,rmpat.AllocatedPeriod,waveconfig.NumSubframes,bwp(ch.BWP).SubcarrierSpacing);
                    rmallocsymbols = reshape(symbperslot*rmallocslots + rmpat.AllocatedSymbols',1,[]);    
                    rs.Symbols = rmallocsymbols;    % OFDM symbols associated with reserved PRB (0-based indices, spanning one or more slots)

                    rs.Period = [];                % Total number of slots in the pattern period (empty means don't repeat)              
                    ch.Reserved(end+1) = rs;       % Configure the PDSCH channel with the pattern                 
                end       
            end
        end

        % Display related PRB level processing
        %
        % Calculate the *PRB* linear indices of all the PDSCH instances, primarily
        % for display purposes here.
        % This is performed by marking the allocated PRB in an empty PRB grid 
        % for the BWP in the entire waveform period, subtracting out the reserved
        % part then find the indices that have been used

        % Create an empty BWP spanning the length of the waveform
        rgrid = zeros(size(ResourceGrids{ch.BWP}));

        % Mark the PRB/symbols associated with all the PDSCH instances in this sequence
        for ns=allocatedSlots 
          rgrid(1+ch.AllocatedPRB,1+symbperslot*ns+ch.AllocatedSymbols) = 1;
        end 

        % Clear out all the 'reserved' PRB from the grid
        for rm = ch.Reserved            
            % Expand the symbols indices by the period
            symbols = expandbyperiod(rm.Symbols,rm.Period*symbperslot,(1+max(allocatedSlots))*symbperslot);
            ind = reSub2Ind(bwp(ch.BWP).NRB,rm.PRB,symbols);
            % Knock out the reserved PRB
            rgrid(1+ind)=0;
        end

        % Identify all the indices that remain 
        pdschindices = find(rgrid);

        % Mark the used PDSCH locations in the PRB grid
        ResourceGrids{ch.BWP}(pdschindices) = chplevel.PDSCH; %#ok<FNDSB>

        % Waveform generation RE level processing 
        %
        % The hPDSCHResources uses a slot-level set of parameters so map the
        % relevant parameter from the waveform level down to the slot level
        nrb = bwp(ch.BWP).NRB;
        ch.PRBSet = ch.AllocatedPRB;
        ch.SymbolSet = ch.AllocatedSymbols;     
        crboffset = carriers(bwp(ch.BWP).CarrierIdx).RBStart + bwp(ch.BWP).RBOffset; % First CRB associated with start of BWP
        ch.PRBRefPoint = crboffset; % First CRB associated with start of BWP

        % Create a data source for this PDSCH sequence
        datasource = hVectorDataSource(ch.DataSource);

        % Configure the DL-SCH processing object for this PDSCH sequence
        if trcoding
            dlsch.TargetCodeRate = ch.TargetCodeRate;
        end

        % Storage for PDSCH instance information
        datastore = [];

        % Initialize modinfo
        modinfo = struct('Gd',0,'G',0,'NREPerPRB',0,'DMRSSymbolSet',[],'CDMGroups',[],'PTRSSymbolSet',[],'CDMLengths',[]);

        % Loop over all the allocated slots
        for i = 1:length(allocatedSlots)

            % Get current slot number 
            s = allocatedSlots(i);
            ch.NSlot = s;

            % Create an empty slot grid to contain a single PDSCH instance
            slotgrid = zeros(12*nrb,symbperslot,ch.NLayers);

            ch.ReservedREs = [];
            if isfield(ch,'CSIRSREs')
                ch.ReservedREs = ch.CSIRSREs{s+1};
            end

            % Get the slot-oriented PDSCH indices, DM-RS indices and DM-RS symbol values
            [pdschREindices,dmrsREindices,dmrsSymbols,ptrsREindices,ptrsSymbols,modinfo] = ...
                hPDSCHResources(struct('NRB',bwp(ch.BWP).NRB,'CyclicPrefix',bwp(ch.BWP).CyclicPrefix,...
                'SubcarrierSpacing',bwp(ch.BWP).SubcarrierSpacing,'RBOffset',crboffset),ch);

            if trcoding
                % Get the RV value for this transmission instance
                rvidx = mod(i-1,length(ch.RVSequence))+1;
                rv = ch.RVSequence(rvidx);

                % For the first RV in a sequence, get a new transport block from
                % the data source and pass it to the DL-SCH processing
                if rvidx == 1
                    trblksize = nrTBS(ch.Modulation,ch.NLayers,length(ch.PRBSet),modinfo.NREPerPRB,ch.TargetCodeRate,ch.Xoh_PDSCH);
                    trblk = datasource.getPacket(trblksize);
                    setTransportBlock(dlsch,trblk);
                end

                % DL-SCH processing to create a codeword
                codeword = dlsch(ch.Modulation,ch.NLayers,modinfo.G,rv);
            else
                % If transport coding is not enabled then get the codeword
                % directly from the data source
                codeword = datasource.getPacket(modinfo.G);
                rv = [];
                trblk = [];
            end

            % PDSCH processing to create the PDSCH QAM symbols
            nID = ch.NID;
            symbols = nrPDSCH(codeword,ch.Modulation,ch.NLayers,nID,ch.RNTI);

            % Write the PDSCH, DM-RS and PT-RS symbols in the slot grid
            slotgrid(pdschREindices) = symbols*db2mag(ch.Power);
            slotgrid(dmrsREindices) = dmrsSymbols*db2mag(ch.Power+ch.PowerDMRS);
            ptpower = 1;
            if ~isempty(ptrsREindices)
                ptpower = db2mag(ch.Power+ch.PowerPTRS);
                slotgrid(ptrsREindices) = ptrsSymbols*db2mag(ch.Power+ch.PowerPTRS);
            end
            
            % Combine PDSCH instance with the rest of the BWP grid
            ResourceElementGrids{ch.BWP}(:,s*symbperslot+(1:symbperslot),1:ch.NLayers) = ResourceElementGrids{ch.BWP}(:,s*symbperslot+(1:symbperslot),1:ch.NLayers) + slotgrid; 

            % Capture resource info for this PDSCH instance
            datastore(i).NSlot = ch.NSlot; %#ok<AGROW>
            datastore(i).TransportBlockSize = length(trblk); %#ok<AGROW>
            datastore(i).TransportBlock = trblk;  %#ok<AGROW>
            datastore(i).RV = rv;  %#ok<AGROW>
            datastore(i).Codeword = codeword;  %#ok<AGROW>
            datastore(i).G = modinfo.G; %#ok<AGROW>
            datastore(i).Gd = modinfo.Gd; %#ok<AGROW>
            datastore(i).ChannelIndices = pdschREindices; %#ok<AGROW>
            datastore(i).ChannelSymbols = symbols*db2mag(ch.Power); %#ok<AGROW>
            datastore(i).DMRSSymbolSet = modinfo.DMRSSymbolSet;  %#ok<AGROW>           
            datastore(i).DMRSIndices = dmrsREindices; %#ok<AGROW>
            datastore(i).DMRSSymbols = dmrsSymbols*db2mag(ch.Power+ch.PowerDMRS); %#ok<AGROW>
            datastore(i).PTRSSymbolSet = modinfo.PTRSSymbolSet;  %#ok<AGROW>
            datastore(i).PTRSIndices = ptrsREindices; %#ok<AGROW>
            datastore(i).PTRSSymbols = ptrsSymbols*ptpower; %#ok<AGROW>
            
        end

        % Capture all resources info for this PDSCH sequence
        name = "";
        if isfield(pdsch(nch),'Name')
            name = pdsch(nch).Name;
        end
        waveinfo.PDSCH(nch).Name = name;
        waveinfo.PDSCH(nch).CDMLengths = modinfo.CDMLengths;
        waveinfo.PDSCH(nch).Resources = datastore;
        
    % End of PDSCH sequence processing         
    end
    
    % Create a new figure to display the plots
    % Display *PRB* resource grids
    % Map the BWPs into carrier sized PRB grids for display
    if waveconfig.DisplayGrids
        figure;
        for bp = 1:length(ResourceGrids)

            % Mark the unused RE in the overall BWP, relative to the carrier, so that
            % it is easier to see with respect to the complete carrier layout
            bgrid = ResourceGrids{bp};
            cgrid = zeros(carriers((bwp(bp).CarrierIdx)).NRB, size(bgrid,2));
            bgrid(bgrid==0) = 0.15;
     
            % Write the BWP into the grid representing the carrier
            cgrid(bwp(bp).RBOffset + (1:size(bgrid,1)),:) = bgrid;

            % Mark the SS blocks     
            nsymbolsperhalfframe = ssbreserved(bwp(bp).CarrierIdx).Period*symbolsPerSlot(bwp(bp));
            symbols = expandbyperiod(ssbreserved(bwp(bp).CarrierIdx).Symbols,nsymbolsperhalfframe,size(cgrid,2));                    
            cgrid(ssbreserved(bwp(bp).CarrierIdx).PRB+1,symbols+1) = chplevel.SS_Burst;
            
            % Plot the PRB BWP grid (relative to the carrier)
            cscaling = 40;
            subplot(length(ResourceGrids),1,bp)
            im = image(cscaling*cgrid); axis xy; title(sprintf('BWP %d in Carrier (SCS=%dkHz). PDSCH, PDCCH and CORESET Location',bp,bwp(bp).SubcarrierSpacing)); xlabel('Symbols'); ylabel('Carrier RB');
            cmap = parula(64);
            colormap(im.Parent,cmap);
            
            % Add a channel legend to the first BWP plot (applies to all)
            if bp == 1 
                % Extract channel names and color marker levels
                fnames = strrep(fieldnames(chplevel),'_',' ');
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
        if isfield(waveconfig,'Windowing') && ~isempty(waveconfig.Windowing)
            windowing = waveconfig.Windowing;
        else
            windowing = getfield(hOFDMInfo(struct('NDLRB',carrierObj.NSizeGrid,'SubcarrierSpacing',carrierObj.SubcarrierSpacing)),'Windowing');
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
    
    % Add SS burst sequence
    if ssburst.Enable
        % The hSSBurst function creates a 5ms half frames of data and the waveform 
        % is parameterized in terms of 1ms subframes so we can work out how
        % many complete instances are required, then generate and extract 
        % portion required in the output waveform
        ssburst.SampleRate = maxsr;
        ssburst.DisplayBurst = false;
        ssburst.NCellID = waveconfig.NCellID;

        % Number of complete half frames required to cover the waveform
        nhframes = ceil(waveconfig.NumSubframes/5);  

        % Burst waveform variable
        burstwaveform = [];
        for i=0:nhframes-1
            % Create the half frame sequences and concatenate
            ssburst.NHalfFrame = mod(i,2);
            ssburst.NFrame = fix(i/2);
            burstwaveform = [burstwaveform;hSSBurst(ssburst)]; %#ok<AGROW>
        end
        burstwaveform = db2mag(ssburst.Power)*burstwaveform;

        % Combine SS burst part with the rest of the waveform 
        waveform = waveform + burstwaveform(1:size(waveform,1),:);
    end
    
    winfo.WaveformResources = waveinfo;
    
% End of main function
end


% Expand 'd' by amount 'e', with optional non-unity strides and exclusion 
function expanded = expander(d,e,s,o,excl)
    if nargin < 5
        excl = 0;
    end
    if nargin < 4
        o = 0;
    end
    if nargin < 3
        s = 1;
    end
    eseq = (o:s:e-1)';
    if excl
        eseq = setdiff((0:e-1)',eseq);
    end
    expanded = reshape(reshape(d,1,[]) + eseq,1,[]);  % Use column expansion
end

% Expand 's' values with respect to period 'd', up to value 'nsf' (optionally accounting for the SCS)
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
        sp = reshape(s(s<p),[],1)+p*(0:ceil(ts/max(p,1))-1);
    end
    if ~isempty(ts)
        sp = reshape(sp(sp < ts),1,[]);            % Trim any excess
    else
        sp = ones(1,0);
    end
end

% Turn sets of PRB/symbol subscripts into linear indices, supporting PRB cell arrays
function ind = reSub2Ind(nrb,prb,symbols)

    % Use column expansion
    if iscell(prb)
        slen = min(length(prb),length(symbols));
        ind = cell2mat(cellfun(@(x,y)reshape(reshape(x,[],1)+nrb*reshape(y,1,[]),1,[]),prb(1:slen),num2cell(symbols(1:slen)),'UniformOutput',false));
    else
        ind = reshape(reshape(prb,[],1) + nrb*reshape(symbols,1,[]),1,[]);
    end
    
end

% Establish the number of symbols per slot from the cyclic prefix
function symbs = symbolsPerSlot(config)

    if isstruct(config)
        config = config.CyclicPrefix;
    end
    cpoptions = {'Normal','Extended'};
    symbs = sum(strcmpi(config,cpoptions) .* [14 12]);

end

% Map subcarriers/symbols from one numerology to another (source fs -> target ft)
function [mappedPRB,mappedSymbols] = mapNumerology(subcarriers,symbols,fs,ft) %#ok<DEFNU>

    % Project the subcarrier indices into target PRB
    mappedPRB = unique(fix(subcarriers*fs/(ft*12)),'stable');
    
    if (ft < fs)
        % If ft/fs < 1, reduction
        mappedSymbols = unique(fix(symbols*ft/fs),'stable');
    else
        % Else, repetition by ft/fs
        mappedSymbols = reshape((0:(ft/fs-1))' + symbols(:)'*ft/fs,1,[]);
    end
    
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

% Get the valid, usable PRB associated with a CORESET definition This
% function returns the PRB associated with the CORESET PRB which overlap
% the AllocatedPRB parameter. Accounts for a non-multiple-of-6 NStartBWP,
% by offsetting the CORESET to the next multiple-of-6 RB.
function [coresetprb,nblocks] = getCORESETPRB(carrier,bwp,coreset)

    nstartbwp = carrier.RBStart + bwp.RBOffset;     % Starting CRB of the BWP containing the CORESET
    offsetprb = mod(6-nstartbwp,6);                 % Offset from the BWP start to first complete block of 6 PRB

    % Turn each PRB into the associated block number that contains it

    % Block 0 should associate with BWP PRB of first complete block of 6, and
    % then expand the resulting set into the associated PRB
    nblocks = unique(floor(coreset.AllocatedPRB/6));                 % Indices of blocks associated by the CORESET AllocatedPRB parameter
    nblockmax = fix((bwp.NRB - offsetprb)/6);                        % Maximum number of complete block that can be contained in the BWP
    nblocks = nblocks( (nblocks >= 0) & (nblocks < nblockmax) );     % Complete block set contained within the BWP
    coresetprb = expander(6*nblocks,6) + offsetprb;                  % PRB associated with the resulting allocated CORESET

end

% Get the set of RB level resources in each SCS carrier that overlap 
% with the SS burst 
function [reserved, ssburst] = ssburstResources(ssburst,carriers)
   
    % Get dependent SS block/burst parameters from the burst config 
    modinfo = nrOFDMInfo(275,240);
    ssburst.SampleRate = modinfo.SampleRate;  % This parameter is required by the burst info function
    ssbinfo = hSSBurstInfo(ssburst);
    ssburst.SubcarrierSpacing = ssbinfo.SubcarrierSpacing; % Write the SCS back into the config
    
    % Create the output structure
    patternnames = arrayfun(@(x,y)sprintf('SS burst resources (SCS=%d kHz) in SCS carrier %d (NRB=%d, SCS=%d kHz)',...
                                    ssburst.SubcarrierSpacing, x, y.NRB, y.SubcarrierSpacing),1:length(carriers),carriers,...
                                    'UniformOutput',false);   
    reserved = struct('Name',patternnames,'PRB',[],'Symbols',[],'Period',[]);  
    
    % If the burst is not active then return early
    if ~ssburst.Enable
        return; 
    end
      
    % Establish the maximum carrier SCS configured
    % The center of the overall waveform will be on k0 of this SCS carrier 
    carrierscs = [carriers.SubcarrierSpacing]; 
    [maxcarrierscs,maxidx] = max(carrierscs); %#ok<ASGLU>
    
    % Relative center frequencies of first and last carriers of SS block
    bnds = [-120 120-1]*ssburst.SubcarrierSpacing*1e3 + ssburst.FrequencySSB;
    % Center frequency of DC subcarrier of combined SCS waveform, relative to 'point A' 
    f0 = ((carriers(maxidx).RBStart + fix(carriers(maxidx).NRB/2))*12 + 6*mod(carriers(maxidx).NRB,2))*carriers(maxidx).SubcarrierSpacing*1e3;
    % Relative center frequencies of SS block edge carriers in waveform (max SCS carrier)
    f = f0 + bnds;
    
    % Calculate SCS carrier oriented RB indices 
    for i = 1:length(carriers)
       crb = fix(f/(12*carriers(i).SubcarrierSpacing*1e3)) - carriers(i).RBStart;
       reserved(i).PRB = crb(1):crb(2);
    end
    
    % Burst symbols in active half frame
    symbols = reshape(ssbinfo.OccupiedSymbols',1,[])-1;
    fs = ssburst.SubcarrierSpacing;     % 'Source' SCS
    
    % Calculate SCS carrier oriented symbols indices 
    for i = 1:length(carriers)

        % Project the SS burst symbols into the SCS carrier
        ft = carriers(i).SubcarrierSpacing; % 'Target' SCS  
        if (ft < fs)
            % If ft/fs < 1, reduction
            mappedSymbols = unique(fix(symbols*ft/fs),'stable');
        else
            % Else, repetition by ft/fs
             mappedSymbols = reshape((0:(ft/fs-1))' + symbols(:)'*ft/fs,1,[]);
        end
        reserved(i).Symbols = mappedSymbols;
        
        % Set the half frame periodicity, relative to the SCS carrier
        period = ssburst.SSBPeriodicity*ft/15; % Number of slots in SS burst period
        reserved(i).Period = period;
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
        titleStr = [num2str(cbw) 'MHz Channel, '];
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

% TS 38.104 Section 5.3.3 
% Minimum guardband and transmission bandwidth configuration
function guardband = getGuardband(fr,cbw,scs)

    % Table 5.3.3-1: Minimum guardband [kHz] (FR1) (TS 38.104)
    cbwFR1    = [    5     10    15     20     25     30     40     50     60     70     80     90    100];
    guardsFR1 = [242.5  312.5 382.5  452.5  522.5  592.5  552.5  692.5    NaN    NaN    NaN    NaN    NaN; ...
                 505.0  665.0 645.0  805.0  785.0  945.0  905.0 1045.0  825.0  965.0  925.0  885.0  845.0; ...
                   NaN 1010.0 990.0 1330.0 1310.0 1290.0 1610.0 1570.0 1530.0 1490.0 1450.0 1410.0 1370.0];
    scsFR1 = [15 30 60].';

    % Table 5.3.3-2: Minimum guardband [kHz] (FR2) (TS 38.104)
    cbwFR2    = [  50  100  200  400];
    guardsFR2 = [1210 2450 4930  NaN; ...
                 1900 2420 4900 9860];
    scsFR2 = [60 120].';

    % return value in MHz
    if (strcmpi(fr,'FR1'))
        guardband = guardsFR1(scsFR1==scs,cbwFR1==cbw) / 1e3;
    else % FR2
        guardband = guardsFR2(scsFR2==scs,cbwFR2==cbw) / 1e3;
    end
    
end

function validateCSIRSAllocationParameters(csirs,nsig,bwp,bwpID)
    % Validate number of RBs allocated for CSI-RS
    CSIRSRBs = csirs.NumRB;
    flagRB = (CSIRSRBs > bwp.NRB);
    incorrectRB = CSIRSRBs(find(flagRB,1));
    if any(flagRB)
        error('For CSIRS %d, number of RB''s allocated to CSI-RS resource %d (%d) exceeds the number of RB''s (%d) in BWP %d.',nsig,find(flagRB,1),incorrectRB,bwp.NRB,bwpID);
    end

    % Validate CSI-RS RB offset relative to BWP
    BWPRBs = bwp.NRB;
    BWPRBOffset = bwp.RBOffset; % BWP offset relative to carrier
    CSIRSRBOffset = csirs.RBOffset; % CSI-RS offset relative to carrier
    minOffset = BWPRBOffset;
    maxOffset = BWPRBOffset + (BWPRBs - CSIRSRBs);
    flagOffset = (CSIRSRBOffset >= minOffset & CSIRSRBOffset <= maxOffset);
    incorrectOffset = CSIRSRBOffset(find(~flagOffset,1));
    if ~(all(flagOffset))
        error('For CSIRS %d, RB offset value of CSI-RS resource %d (%d) must be with in the range [%d %d] for BWP %d',nsig,find(~flagOffset,1),incorrectOffset,minOffset,maxOffset,bwpID);
    end
end

function val = applyScalarExpansion(val,n,propName,nsig)
    % Apply scalar expansion
    if ~((numel(val) == 1) || (numel(val) == n))
        error([propName ', field of CSIRS %d, must be either a scalar or a vector of length equals to the number of CSI-RS resources configured (%d)'],nsig,n);
    end
    if (numel(val) == 1)
        val = repmat(val,n,1);
    end
end

function coreset = validateCoreset(coreset)
% Validate CORESET parameters for interleaving.

    if ~isfield(coreset,'CCEREGMapping')
        % Add field, with default value
        coreset.CCEREGMapping = 'noninterleaved'; % default
    else
        if strcmp(coreset.CCEREGMapping,'interleaved')
            % Cross-check Duration and REGBundleSize(L)
            if coreset.Duration==3
                if coreset.REGBundleSize==2
                    error('With a CORESET duration of 3, REGBundleSize must be 3 or 6.');
                end
            else
                if coreset.REGBundleSize==3
                    error('With a CORESET duration of 1 or 2, REGBundleSize must be 2 or 6.');
                end
            end

            % Cross-check allocated REGs with REGBundlesize(L) and
            % InterleaverSize (R)
            numREGS = numel(coreset.AllocatedPRB)*6*coreset.Duration;
            C = numREGS/(coreset.REGBundleSize*coreset.InterleaverSize);
            if floor(C)~=C
                error(['Invalid interleaving configuration. Select a ',...
                       'different REGBundleSize or InterleaverSize.']);
            end
        end
    end
end