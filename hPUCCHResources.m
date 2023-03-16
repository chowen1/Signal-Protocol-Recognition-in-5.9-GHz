function [indices,dmrs] = hPUCCHResources(ue,pucch)
%hPUCCHResources 5G NR PUCCH and DM-RS resource element indices and DM-RS values
%   [INDICES,DMRS] = hPUCCHResources(UE,PUCCH) returns the PUCCH indices
%   information in a structure INDICES and the demodulation reference
%   signal DMRS based on the information provided in input structures UE
%   and PUCCH.
%
%   The UE-wide settings input, UE, must be a structure with the following 
%   fields:
%   CyclicPrefix      - Cyclic prefix ('Normal'(default),'Extended')
%   RBOffset          - Resource block offset. Default is 0
%   SubcarrierSpacing - Subcarrier Spacing (kHz) (15(default),30,60,120,240)
%
%   Note that the RBOffset parameter of input UE is used in the calculation
%   of DM-RS symbols for PUCCH format 2.
%
%   The input PUCCH must be a structure with the following fields:
%   NSLOT                - Absolute slot number of PUCCH transmission
%   PUCCHFormat          - PUCCH format number (0...4)
%   StartSymbol          - Starting symbol of PUCCH transmission
%   NrOfSymbols          - Number of OFDM symbols allocated for PUCCH
%                          transmission. It is either 1 or 2 for PUCCH
%                          formats 0 and 2. For formats 1, 3 and 4, it is
%                          in range 4 to 14.
%   StartPRB             - Index of first PRB prior to frequency hopping or
%                          for no frequency hopping
%   SecondHopPRB         - Index of first PRB after frequency hopping
%   NrOfRB               - Number of resource blocks allocated for PUCCH
%                          transmission for formats 2 and 3.
%   Modulation           - Modulation scheme for format 3 and 4. It must be
%                          one of the set {'pi/2-BPSK', 'QPSK'}.
%   HoppingId            - Scrambling identity. It is equal to the
%                          higher-layer parameter hoppingId (0...1023), if
%                          configured or else, it is equal to the physical
%                          layer cell identity, NCellID (0...1007). It is
%                          used in DM-RS sequence generation for formats
%                          1/3/4.
%   NIDDMRS              - Scrambling identity for DM-RS. It is used for
%                          format 2. The value is equal to higher-layer
%                          parameter scramblingID0 (0...65535), else, it is
%                          equal to the physical layer cell identity
%                          (0...1007).
%   GROUPHOPPING         - Group hopping configuration. It is one of the
%                          set {'neither','enable','disable'} provided by
%                          higher-layer parameter pucch-GroupHopping. It
%                          must be provided for all formats except format
%                          number 2.
%   IntraSlotFreqHopping - Intra-slot frequency hopping. It is one of the
%                          set {'enabled','disabled'} provided by
%                          higher-layer parameter
%                          intraSlotFrequencyHopping.
%   InitialCS            - Initial cyclic shift (m_0). It is in range 0 to
%                          11, provided by higher-layer parameter
%                          initialCyclicShift. It is provided for formats 0
%                          and 1.
%   OCCI                 - Orthogonal cover code index. It is in range 0 to
%                          6, for format 1 and the valid range depends on
%                          the number of OFDM symbols per hop which contain
%                          control information. For format 4, it must be
%                          less than the spreading factor.
%   SpreadingFactor      - Spreading factor for PUCCH format 4. It must be
%                          either 2 or 4.
%   AdditionalDMRS       - Additional DM-RS configuration flag for PUCCH
%                          format 3 and 4. It is either true or false,
%                          provided by higher-layer parameter
%                          additionalDMRS. The value set to true indicates,
%                          2 DM-RS symbols per hop, if both the hops are
%                          more than 4 symbols when IntraSlotFreqHopping is
%                          enabled and 4 DM-RS symbols with more than 9
%                          symbols when IntraSlotFreqHopping is disabled.
%
%   The fields of output structure INDICES are
%   AllocatedSymbols         - OFDM symbols indices in a slot allocated for
%                              PUCCH UCI
%   AllocatedSubcarriers     - Subcarrier indices allocated for PUCCH UCI
%   AllocatedDMRSSymbols     - OFDM symbol indices allocated for PUCCH DM-RS
%   AllocatedDMRSSubcarriers - Subcarrier indices allocated for PUCCH DM-RS
%   Gd                       - Symbol capacity
%   G                        - Bit capacity
%
%   Example:
%   % Get the DM-RS symbols for a specific PUCCH configuration
%
%   % Set the UE structure
%   ue = struct('CyclicPrefix','Normal','RBOffset',0,'SubcarrierSpacing',15);
%
%   % Specific PUCCH resource configuration
%   pucch = struct();
%   pucch.NSlot = 0;
%   pucch.PUCCHFormat = 3;
%   pucch.StartSymbol = 5;
%   pucch.NrOfSymbols = 8;
%   pucch.StartPRB = 0;
%   pucch.SecondHopPRB = 1;
%   pucch.NrOfRB = 15;
%   pucch.Modulation = 'pi/2-BPSK';
%   pucch.HoppingId = 0; % For formats 1/3/4
%   pucch.NIDDMRS = []; % Only for format 2
%   pucch.GroupHopping = 'enable';
%   pucch.IntraSlotFreqHopping = 'enabled';
%   pucch.InitialCS = 0;
%   pucch.OCCI = 0;
%   pucch.SpreadingFactor = 2;
%   pucch.AdditionalDMRS = 1;
%
%   % PUCCH indices and DM-RS sequence
%   [indices,dmrs] = hPUCCHResources(ue,pucch);
%
%   See also nrPUCCH0, nrPUCCH1, nrPUCCH2, nrPUCCH3, nrPUCCH4.

% Copyright 2019 The MathWorks, Inc.

    % Get the PUCCH indices
    indices = pucchIndices(pucch);

    % Get the number of symbols per slot
    if ~isfield(ue,'CyclicPrefix')
        ue.CyclicPrefix = 'Normal';
    end

    % Check the RB offset field
    if ~isfield(ue,'RBOffset')
        ue.RBOffset = 0;
    end

    % Capture the SCS and number of slots in a 10ms frame
    if isfield(ue,'SubcarrierSpacing')
        scs = ue.SubcarrierSpacing;
    else
        scs = 15;  % Default to 15kHz SCS
    end
    slotsperframe = 10*(scs/15);

    % Initialize parameters
    nID = pucch.HoppingId;
    nSlot = mod(pucch.NSlot,slotsperframe);

    % PUCCH DM-RS
    switch pucch.PUCCHFormat
        case 0
            % For format 0, there is no DM-RS
            dmrs = [];
        case 1
            dmrs = pucchDMRSF1([pucch.StartSymbol pucch.NrOfSymbols],...
                ue.CyclicPrefix,nSlot,nID,pucch.GroupHopping,pucch.InitialCS,pucch.IntraSlotFreqHopping,pucch.OCCI);
        case 2
            dmrs = pucchDMRSF2(indices.AllocatedDMRSSubcarriers,...
                indices.AllocatedDMRSSymbols,ue.CyclicPrefix,nSlot,pucch.NIDDMRS,ue.RBOffset);
        case 3
            dmrs = pucchDMRSF34(pucch.NrOfRB,[pucch.StartSymbol pucch.NrOfSymbols],...
                ue.CyclicPrefix,nSlot,nID,pucch.GroupHopping,pucch.IntraSlotFreqHopping,pucch.AdditionalDMRS);
        otherwise
            dmrs = pucchDMRSF34(1,[pucch.StartSymbol pucch.NrOfSymbols],ue.CyclicPrefix,nSlot,nID,...
                pucch.GroupHopping,pucch.IntraSlotFreqHopping,pucch.AdditionalDMRS,pucch.SpreadingFactor,pucch.OCCI);
    end

end

function indices = pucchIndices(resInfo)
%pucchIndices PUCCH indices information
%   INDICES = pucchIndices(RESINFO) provides the indices related
%   information (0-based) of PUCCH symbols and PUCCH DM-RS symbols in a
%   structure INDICES based on the input structure resource information
%   RESINFO. The RESINFO contains the following fields:
%   Modulation           - Modulation scheme. It must be one of the set
%                          {'pi/2-BPSK', 'QPSK'}.
%   PUCCHFormat          - PUCCH format (0...4)
%   StartSymbol          - Starting OFDM symbol index for PUCCH
%                          transmission
%   NrOfSymbols          - Number of OFDM symbols allocated for PUCCH
%                          transmission. For formats 1, 3 and 4, it is in
%                          range 4 to 14. For formats 0 and 2, it is either
%                          1 or 2.
%   StartPRB             - Index of first PRB prior to frequency hopping or
%                          for no frequency hopping
%   SecondHopPRB         - Index of first PRB after frequency hopping
%   IntraSlotFreqHopping - Intra-slot frequency hopping. It is one of the set
%                          {'enabled','disabled'} provided by higher-layer
%                          parameter intraSlotFrequencyHopping.
%   NrOfRB               - The number of resource blocks associated with
%                          the PUCCH format 2 or 3 transmission. Nominally
%                          the value is one of the set
%                          {1,2,3,4,5,6,8,9,10,12,15,16}.
%   SpreadingFactor      - Spreading factor for PUCCH format 4. It must be
%                          either 2 or 4
%   AdditionalDMRS       - Additional DM-RS flag for PUCCH format 3 and 4.
%                          It is either true or false, provided by
%                          higher-layer parameter additionalDMRS. The value
%                          set to true indicates, 2 DM-RS symbols per hop,
%                          if both the hops are more than 4 symbols when
%                          FREQHOPPING is enabled and 4 DM-RS symbols with
%                          more than 9 symbols when FREQHOPPING is
%                          disabled.
%
%   The output fields of INDICES are
%   AllocatedSymbols         - OFDM symbols indices in a slot allocated for
%                              PUCCH UCI
%   AllocatedSubcarriers     - Subcarrier indices allocated for PUCCH UCI
%   AllocatedDMRSSymbols     - OFDM symbol indices allocated for PUCCH DM-RS
%   AllocatedDMRSSubcarriers - Subcarrier indices allocated for PUCCH DM-RS
%   Gd                       - Symbol capacity
%   G                        - Bit capacity
%
%   Example:
%   % Get the indices information for a PUCCH specific configuration
%
%   pucch = struct();
%   pucch.Modulation = 'pi/2-BPSK';
%   pucch.PUCCHFormat = 3;
%   pucch.StartSymbol = 5;
%   pucch.NrOfSymbols = 8;
%   pucch.StartPRB = 0;
%   pucch.SecondHopPRB = 1;
%   pucch.IntraSlotFreqHopping = 'enabled';
%   pucch.NrOfRB = 15;
%   pucch.SpreadingFactor = 2;
%   pucch.AdditionalDMRS = 1;
%
%   indices = pucchIndices(pucch)

    % Check with intra-slot frequency hopping
    if (strcmpi(resInfo.IntraSlotFreqHopping,'enabled') && (isempty(resInfo.SecondHopPRB) || ~isfield(resInfo,'SecondHopPRB')))...
            || strcmpi(resInfo.IntraSlotFreqHopping,'disabled')
        resInfo.SecondHopPRB = resInfo.StartPRB;
    end

    % Initialize few parameters
    sf = 1;
    nrOfRB = 1;
    qm = 1;
    if resInfo.PUCCHFormat == 4
        sf = resInfo.SpreadingFactor;
    end
    if resInfo.PUCCHFormat == 2 || resInfo.PUCCHFormat == 3
        nrOfRB = resInfo.NrOfRB;
    end

    % Symbol and subcarrier indices of UCI and DM-RS on PUCCH
    switch resInfo.PUCCHFormat
        case 0
            % PUCCH indices
            l = resInfo.StartSymbol:resInfo.StartSymbol+resInfo.NrOfSymbols-1;
            if resInfo.NrOfSymbols == 2
                k = [resInfo.StartPRB resInfo.SecondHopPRB]*12 + (0:11)';
            else
                k = resInfo.StartPRB*12 + (0:11)';
            end

            % DM-RS indices
            ldmrs = [];
            kdmrs = [];

        case 1
            % PUCCH indices
            l = resInfo.StartSymbol+1:2:resInfo.StartSymbol+resInfo.NrOfSymbols-1;
            if strcmpi(resInfo.IntraSlotFreqHopping,'disabled')
                nSF0 = floor(resInfo.NrOfSymbols/2);
            else
                nSF0 = floor(resInfo.NrOfSymbols/4);
            end
            nSF1 = floor(resInfo.NrOfSymbols/2) - nSF0;
            k = [repmat(resInfo.StartPRB,1,nSF0) repmat(resInfo.SecondHopPRB,1,nSF1)]*12 + (0:11)';

            % DM-RS indices
            ldmrs = resInfo.StartSymbol:2:resInfo.StartSymbol+resInfo.NrOfSymbols-1;
            nSF = ceil(resInfo.NrOfSymbols/2);
            if strcmpi(resInfo.IntraSlotFreqHopping,'disabled')
                nSF0 = nSF;
            else
                if rem(resInfo.NrOfSymbols,2) == 1
                    nSF0 = floor(nSF/2);
                else
                    nSF0 = ceil(nSF/2);
                end
            end
            nSF1 = nSF - nSF0;
            kdmrs = [repmat(resInfo.StartPRB,1,nSF0) repmat(resInfo.SecondHopPRB,1,nSF1)]*12 + (0:11)';

        case 2
            qm = 2;
            if resInfo.NrOfSymbols == 2
                secondHopPRB = resInfo.SecondHopPRB;
            else
                secondHopPRB = [];
            end

            % DM-RS indices
            ldmrs = resInfo.StartSymbol:resInfo.StartSymbol+resInfo.NrOfSymbols-1;
            kdmrs = [resInfo.StartPRB secondHopPRB]*12+(1:3:nrOfRB*12)';

            % PUCCH indices
            l = ldmrs;
            k = setdiff([resInfo.StartPRB secondHopPRB]*12+(0:nrOfRB*12-1)',kdmrs,'rows');

        otherwise % 3 or 4
            qm = sum(strcmpi(resInfo.Modulation,{'pi/2-BPSK','QPSK'}).*[1 2]);

            % DM-RS indices
            ldmrs = getDMRSSymIndicesF34([resInfo.StartSymbol resInfo.NrOfSymbols],...
                resInfo.IntraSlotFreqHopping,resInfo.AdditionalDMRS);
            kdmrs1 = 0:nrOfRB*12-1;
            if strcmpi(resInfo.IntraSlotFreqHopping,'enabled')
                kdmrs = [repmat(resInfo.StartPRB,1,length(ldmrs)/2) repmat(resInfo.SecondHopPRB,1,length(ldmrs)/2)]*12 + kdmrs1';
            else
                kdmrs = repmat(resInfo.StartPRB,1,length(ldmrs))*12 + kdmrs1';
            end

            % PUCCH indices
            lset = resInfo.StartSymbol:resInfo.StartSymbol+resInfo.NrOfSymbols-1;
            l = setdiff(lset,ldmrs);
            if strcmpi(resInfo.IntraSlotFreqHopping,'enabled')
                k = [repmat(resInfo.StartPRB,1,floor(length(l)/2)) repmat(resInfo.SecondHopPRB,1,ceil(length(l)/2))]*12 + kdmrs1';
            else
                k = repmat(resInfo.StartPRB,1,length(l))*12 + kdmrs1';
            end
    end

    % Symbol/Bit capacity
    Gd = numel(k)/sf;
    G = Gd*qm;

    % Combine information into a structure
    indices.AllocatedSymbols = l;
    indices.AllocatedSubcarriers = k;
    indices.AllocatedDMRSSymbols = ldmrs;
    indices.AllocatedDMRSSubcarriers = kdmrs;
    indices.Gd = Gd;
    indices.G = G;

end


function sym = pucchDMRSF1(symAllocation,cp,nslot,nid,groupHopping,initialCS,freqHopping,occi)
%pucchDMRSF1 DM-RS symbols for PUCCH format 1
%   SYM = pucchDMRSF1(SYMALLOCATION,CP,NSLOT,NID,GROUPHOPPING,INITIALCS,FREQHOPPING,OCCI)
%   returns the PUCCH format 1 DM-RS symbols SYM, with the following inputs:
%   SYMALLOCATION - Symbol allocation for PUCCH transmission. It is a
%                   two-element vector, where first element is the symbol
%                   index corresponding to first OFDM symbol of the PUCCH
%                   transmission in the slot and second element is the
%                   number of OFDM symbols allocated for PUCCH
%                   transmission, which is in range 4 and 14.
%   CP            - Cyclic prefix ('normal','extended').
%   NSLOT         - Slot number in radio frame. It is in range 0 to 159 for
%                   normal cyclic prefix for different numerologies. For
%                   extended cyclic prefix, it is in range 0 to 39, as
%                   specified in TS 38.211 Section 4.3.2.
%   NID           - Scrambling identity. It is in range 0 to 1023 if
%                   higher-layer parameter hoppingId is provided, else, it
%                   is in range 0 to 1007, equal to the physical layer cell
%                   identity NCellID.
%   GROUPHOPPING  - Group hopping configuration. It is one of the set
%                   {'neither','enable','disable'} provided by higher-layer
%                   parameter pucch-GroupHopping.
%   INITIALCS     - Initial cyclic shift (m_0). It is in range 0 to 11,
%                   provided by higher-layer parameter initialCyclicShift.
%   FREQHOPPING   - Intra-slot frequency hopping. It is one of the set
%                   {'enabled','disabled'} provided by higher-layer
%                   parameter intraSlotFrequencyHopping.
%   OCCI          - Orthogonal cover code index. It is in range 0 to 6,
%                   provided by higher-layer parameter timeDomainOCC. It is
%                   less than the number of OFDM symbols containing DM-RS
%                   information in both the hops when FREQHOPPING is
%                   enabled and less than number of OFDM symbols containing
%                   DM-RS information when FREQHOPPING is disabled.
%
%   Example:
%   % Get the DM-RS symbols for a specific PUCCH 1 resource.
%   
%   cp = 'normal';
%   symAllocation = [0 4];
%   nslot = 0;
%   nid = 512;
%   groupHopping = 'enable';
%   initialCS = 5;
%   freqHopping = 'enabled';
%   occi = 0;
%   sym = pucchDMRSF1(symAllocation,cp,nslot,nid,groupHopping,initialCS,freqHopping,occi);

    symStart = symAllocation(1);
    nPUCCHSym = symAllocation(2);

    % Get the number of symbols allocated for DM-RS in each hop
    nSF = ceil(nPUCCHSym/2);
    if strcmpi(freqHopping,'disabled')
        nSF0 = nSF;
    else
        if rem(nPUCCHSym,2) == 1
            nSF0 = floor(nSF/2);
        else
            nSF0 = ceil(nSF/2);
        end
    end
    nSF1 = nSF - nSF0;

    seqCS = 0;
    info = nrPUCCHHoppingInfo(cp,nslot,nid,groupHopping,initialCS,seqCS);

    ind = symStart+1:2:14;
    nRBSC = 12;
    lps1 = nrLowPAPRS(info.U(1),info.V(1),info.Alpha(ind(1:nSF0)),nRBSC);
    if strcmpi(freqHopping,'enabled')
        lps2 = nrLowPAPRS(info.U(2),info.V(2),info.Alpha(ind(nSF0+1:nSF0+nSF1)),nRBSC);
        r = [lps1 lps2];
    else
        r = lps1;
    end

    % Get the orthogonal sequence from spreading factor and orthogonal
    % cover code index
    oSeq1 = nr5g.internal.PUCCH1Spreading(nSF0,occi);
    if strcmpi(freqHopping,'disabled')
        oSeq = oSeq1;
    else
        oSeq2 = nr5g.internal.PUCCH1Spreading(nSF1,occi);
        oSeq = [oSeq1 oSeq2];
    end

    % Get the PUCCH format 1 DM-RS sequence
    sym = r.*repmat(oSeq,nRBSC,1);

end

function sym = pucchDMRSF2(dmrsSubcarriers,nSym,cp,nSlot,nID,rbOffset)
%pucchDMRSF2 DM-RS symbols for PUCCH format 2
%   SYM = pucchDMRSF2(DMRSSUBCARRIERS,NSYM,CP,NSLOT,NID,RBOFFSET) returns
%   the DM-RS symbols for the PUCCH format 2, provided the following inputs:
%   DMRSSUBCARRIERS - Subcarrier indices of DM-RS PUCCH format 2 for each
%                     OFDM symbol. It must be a either a column vector or a
%                     matrix of two columns.
%   NSYM            - DM-RS symbol indices in a slot. It must be either a
%                     scalar or a two-element vector.
%   CP              - Cyclic prefix ('normal','extended').
%   NSLOT           - Slot number in radio frame. It is in range 0 to 159
%                     for normal cyclic prefix for different numerologies.
%                     For extended cyclic prefix, it is in range 0 to 39,
%                     as specified in TS 38.211 Section 4.3.2.
%   NID             - Scrambling identity. It is in range 0 to 65535 if
%                     higher-layer parameter scramblingID is provided,
%                     else, it is in range 0 to 1007, equal to the physical
%                     layer cell identity NCellID.
%   RBOFFSET        - Resource block offset.
%
%   Example:
%   % Get the DM-RS symbols for a specific PUCCH format 2 resource
%
%   dmrsSubcarriers = 3*(0:59)+1;
%   nsym = 0;
%   cp = 'normal';
%   nslot = 0;
%   nid = 0;
%   rbOffset = 0;
%   sym = pucchDMRSF2(dmrsSubcarriers,nsym,cp,nslot,nid,rbOffset);


    if strcmpi(cp,'normal')
        symbperslot = 14;
    else
        symbperslot = 12;
    end

    % Get the smallest and largest PRB (0-based) value
    [minRB, maxRB] = bounds(floor(dmrsSubcarriers/12));

    % Construct PRBS for the transmission DM-RS, offsetting the
    % sequence to account for the origin of the BWP
    cinit = mod(2^17*(symbperslot*nSlot+nSym+1)*(2*nID+1)+2*nID,2^31);
    cSeq1 = nrPRBS(cinit(1),2*4*[rbOffset+minRB(1) maxRB(1)-minRB(1)+1]);

    % Modulate the PRBS sequence and get the DM-RS
    dmrs1 = nrSymbolModulate(cSeq1,'QPSK');
    if length(nSym)==2
        cSeq2 = nrPRBS(cinit(2),2*4*[rbOffset+minRB(2) maxRB(2)-minRB(2)+1]);
        dmrs2 = nrSymbolModulate(cSeq2,'QPSK');
    else
        dmrs2 = [];
    end
    sym = [dmrs1 dmrs2];

end

function sym = pucchDMRSF34(Mrb,symAllocation,cp,nslot,nid,groupHopping,freqHopping,additionalDMRS,sf,occi)
%pucchDMRSF34 DM-RS symbols for PUCCH format 3 and 4
%   SYM = pucchDMRSF34(MRB,SYMALLOCATION,CP,NSLOT,NID,GROUPHOPPING,FREQHOPPING,ADDITIONALDMRS,SF,OCCI)
%   returns the DM-RS symbols for PUCCH format 3/4, provided the following
%   inputs:
%   MRB            - The number of resource blocks associated with the
%                    PUCCH format 3 transmission. Nominally the value of
%                    MRB will be one of the set
%                    {1,2,3,4,5,6,8,9,10,12,15,16}.
%   SYMALLOCATION  - Symbol allocation for PUCCH transmission. It is a
%                    two-element vector, where first element is the symbol
%                    index corresponding to first OFDM symbol of the PUCCH
%                    transmission in the slot and second element is the
%                    number of OFDM symbols allocated for PUCCH
%                    transmission, which is in range 4 and 14.
%   CP             - Cyclic prefix ('normal','extended')
%   NSLOT          - Slot number in radio frame. It is in range 0 to 159
%                    for normal cyclic prefix for different numerologies.
%                    For extended cyclic prefix, it is in range 0 to 39, as
%                    specified in TS 38.211 Section 4.3.2.
%   NID            - Scrambling identity. It is in range 0 to 1023 if
%                    higher-layer parameter hoppingId is provided, else, it
%                    is in range 0 to 1007, equal to the physical layer
%                    cell identity NCellID.
%   GROUPHOPPING   - Group hopping configuration. It is one of the set
%                    {'neither','enable','disable'} provided by
%                    higher-layer parameter pucch-GroupHopping.
%   FREQHOPPING    - Intra-slot frequency hopping. It is one of the set
%                    {'enabled','disabled'} provided by higher-layer
%                    parameter intraSlotFrequencyHopping.
%   ADDITIONALDMRS - Additional DM-RS flag. It is either true or false,
%                    provided by higher-layer parameter additionalDMRS. The
%                    value set to true indicates, 2 DM-RS symbols per
%                    hop, if both the hops are more than 4 symbols when
%                    FREQHOPPING is enabled and 4 DM-RS symbols with more
%                    than 9 symbols when FREQHOPPING is disabled.
%   SF             - Spreading factor for PUCCH format 4. It must be either
%                    2 or 4.
%   OCCI           - Orthogonal cover code sequence index. It must be
%                    greater than or equal to zero and less than SF.
%
%   Example:
%   % Get the DM-RS symbols for PUCCH format 3 for normal cyclic prefix with
%   % number of resource blocks as 5, starting symbol of PUCCH allocation
%   % as 2, number of symbols allocated as 11, slot number as 2, scrambling
%   % identity as 1007, group hopping as enable, intra slot frequency
%   % hopping enabled and with additional DM-RS.
%
%   cp = 'normal';
%   Mrb = 5;
%   symAllocation = [2 11];
%   nslot = 2;
%   nid = 1007;
%   groupHopping = 'enable';
%   freqHopping = 'enabled';
%   additionalDMRS = true;
%   sym = pucchDMRSF34(5,[2 11],'normal',2,1007,'enable','enabled',true)

    narginchk(8,10);

    cp = validatestring(cp,{'normal','extended'},'','CP');

    msc = Mrb*12;

    % Check freqHopping
    freqHopping = validatestring(freqHopping,{'enabled','disabled'},'','FREQHOPPING');

    % For format 4
    % Get the cyclic shift index
    if nargin == 10
        if sf == 2
            csTable = [0 6];
        else
            csTable = [0 6 3 9];
        end
        initialCS = csTable(1,occi+1);
    else
        initialCS = 0;
    end
    seqCS = 0;
    info = nrPUCCHHoppingInfo(cp,nslot,nid,groupHopping,initialCS,seqCS);

    l = getDMRSSymIndicesF34(symAllocation,freqHopping,additionalDMRS);

    if strcmpi(freqHopping,'enabled')
        nSF0 = floor(length(l)/2);
    else
        nSF0 = length(l);
    end

    lps1 = nrLowPAPRS(info.U(1),info.V(1),info.Alpha(l(1:nSF0)+1),msc);

    if strcmpi(freqHopping,'enabled')
        lps2 = nrLowPAPRS(info.U(2),info.V(2),info.Alpha(l(nSF0+1:end)+1),msc);
        sym = [lps1 lps2];
    else
        sym = lps1;
    end

end

function l = getDMRSSymIndicesF34(symAllocation,freqHopping,additionalDMRS)
%getDMRSSymIndicesF34 DM-RS OFDM symbol indices for PUCCH formats 3 and 4
%   L = getDMRSSymIndicesF34(SYMALLOCATION,FREQHOPPING,ADDITIONALDMRS)
%   provides the resource element indices (0-based) for the DM-RS symbols of
%   PUCCH formats 3 and 4, considering the following inputs:
%   SYMALLOCATION  - Symbol allocation for PUCCH transmission. It is a
%                    two-element vector, where first element is the symbol
%                    index corresponding to first OFDM symbol of the PUCCH
%                    transmission in the slot and second element is the
%                    number of OFDM symbols allocated for PUCCH
%                    transmission, which is in range 4 and 14.
%   FREQHOPPING    - Intra-slot frequency hopping. It is one of the set
%                    {'enabled','disabled'} provided by higher-layer
%                    parameter intraSlotFrequencyHopping.
%   ADDITIONALDMRS - Additional DM-RS flag. It is either true or false,
%                    provided by higher-layer parameter additionalDMRS. The
%                    value set to true indicates, 2 DM-RS symbols per
%                    hop, if both the hops are more than 4 symbols when
%                    FREQHOPPING is enabled and 4 DM-RS symbols with more
%                    than 9 symbols when FREQHOPPING is disabled.
%
%   Example:
%   % Get the OFDM symbols indices for DM-RS symbols given the starting
%   % symbol of PUCCH allocation as 0, number of symbols allocated as 10,
%   % frequency hopping disabled and with additional DM-RS.
%
%   symAllocation = [0 10];
%   freqHopping = 'enabled';
%   additionalDMRS = true;
%   l = getDMRSSymIndicesF34(symAllocation,freqHopping,additionalDMRS)

    nPUCCHSym = symAllocation(2);
    symIndex = symAllocation(1);

    if nPUCCHSym == 4 
        if strcmpi(freqHopping,'enabled')
            sym = [0 2];
        else
            sym = 1;
        end
    elseif nPUCCHSym <= 9
        indexTable = [0 3;...
                      1 4;...
                      1 4;...
                      1 5;...
                      1 6];
        sym = indexTable(nPUCCHSym-4,:);
    else
        if additionalDMRS
            indexTable = [1 3 6 8;...
                          1 3 6 9;...
                          1 4 7 10;...
                          1 4 7 11;...
                          1 5 8 12];
        else
            indexTable = [2 7;...
                          2 7;...
                          2 8;...
                          2 9;...
                          3 10];
        end
        sym = indexTable(nPUCCHSym-9,:);
    end

    l = symIndex+sym; % 0-based
end