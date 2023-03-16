function out = hMultiplexPUCCH(uciPart1,uciPart2,modulation,nPUCCHSym,dmrsSym,format,nRBOrSF)
%hMultiplexPUCCH UCI multiplexing for PUCCH formats 3 and 4
%   OUT = hMultiplexPUCCH(UCIPART1,UCIPART2,MODULATION,NPUCCHSYM,DMRSSYM,FORMAT,NRBORSF)
%   provides the multiplexed output OUT for UCIPART1 and UCIPART2 as per TS
%   38.212, Section 6.3.1.6, with the following inputs:
%   UCIPART1   - Coded UCI bits that correspond to CSI part1, with or
%                without HARQ-ACK and SR bits. Use empty ([]) to indicate
%                no UCI part 1 transmission.
%   UCIPART2   - Coded UCI bits that correspond to CSI part2. Use empty
%                ([]) to indicate no UCI part 2 transmission.
%   MODULATION - Modulation scheme. It is one of the set
%                {'pi/2-BPSK','QPSK'}.
%   NPUCCHSYM  - Number of OFDM symbols allocated for PUCCH transmission.
%                The value must be in range 4 to 14.
%   DMRSSYM    - PUCCH DM-RS symbol indices, as per Table 6.3.1.6-1, TS
%                38.212. The DM-RS symbol indices depends on the number of
%                OFDM symbols allocated for PUCCH transmission.
%   FORMAT     - Format number of PUCCH. It must be either 3 or 4.
%   NRBORSF    - Number of resource blocks or the spreading factor based on
%                the format. For format 3, the value is number of resource
%                blocks allocated for PUCCH and must be in range 1 and 16.
%                The nominal value for number of resource blocks is one of
%                the set {1,2,3,4,5,6,8,9,10,12,15,16}. For format 4, the
%                value is spreading factor and must be either 2 or 4.
%
%   The modulation order (qm) is 1 for pi/2-BPSK modulation and 2 for QPSK
%   modulation. The number of resource elements for each OFDM symbol
%   carrying UCI information (nSymbolUCI) for format 3 is equal to product
%   of number of subcarriers in a resource block (12) and the number of
%   resource blocks allocated for PUCCH (NRBORSF). Similarly for format 4,
%   it is equal to ratio of number of subcarriers in a resource block (12)
%   and the spreading factor (NRBORSF). The number of OFDM symbols carrying
%   UCI information (nPUCCHSymUCI) is difference of number of OFDM symbols
%   allocated for PUCCH (NPUCCHSYM) and the length of DMRSSYM.
%
%   The length of output OUT is the product of qm, nSymbolUCI and
%   nPUCCHSymUCI. Note that there will be zero entries in OUT when the
%   length of both the inputs is less than the length of OUT.
%
%   Example:
%   % Perform UCI multiplexing of UCI part 1 and part 2 for format 3 with
%   % rate matched lengths 874 and 1070 respectively, for QPSK modulation
%   % with number of OFDM symbols allocated for PUCCH as 11, number of
%   % resource blocks as 9 and the DM-RS symbol locations as 2 and 7.
%
%   uci1 = -1*ones(874,1);
%   uci2 = -2*ones(1070,1);
%   modulation = 'QPSK';
%   nPUCCHSym = 11;
%   dmrssym = [2 7];
%   format = 3;
%   nRB = 9;
%   cw = hMultiplexPUCCH(uci1,uci2,modulation,nPUCCHSym,dmrssym,format,nRB);
%
%   % Check the number of bits that are not uci1 or uci2
%   nnz(~((cw==-1)|(cw==-2)))
%
%   See also nrPUCCH3, nrPUCCH4, nrUCIEncode.

% Copyright 2019 The MathWorks, Inc.

%#codegen

    narginchk(7,7);

    % Validate inputs
    fcnName = 'hMultiplexPUCCH';
    validateattributes(uciPart1,{'numeric','logical'},{'2d'},fcnName,'UCIPART1');
    validateattributes(uciPart2,{'numeric','logical'},{'2d'},fcnName,'UCIPART2');
    modulation = validatestring(modulation,{'QPSK','pi/2-BPSK'},fcnName,'MODULATION');
    validateattributes(nPUCCHSym,{'numeric'},{'scalar','integer','>=',4,'<=',14},fcnName,'NPUCCHSYM');
    validateattributes(dmrsSym,{'numeric'},{'vector','integer'},fcnName,'DMRSSYM');

    if ~((format == 3) || (format == 4))
        error('The PUCCH format (%d) must be either 3 or 4.',format);
    end

    if (format == 3) && (isempty(nRBOrSF) || nRBOrSF < 1 || nRBOrSF > 16)
        error('The number of resource blocks (%d) must be in range 1 and 16.',nRBOrSF);
    end

    if (format == 4) && (isempty(nRBOrSF)|| ~(nRBOrSF == 2 || nRBOrSF == 4))
        error('The spreading factor (%d) must be either 2 or 4.',nRBOrSF);
    end

    % Get the input lengths
    g1 = length(reshape(uciPart1,[],1));
    g2 = length(reshape(uciPart2,[],1));

    % Check if both the inputs are zero and return empty output
    if g1==0 && g2 ==0
        out = zeros(0,1);
        return;
    end

    % Modulation order
    qm = nr5g.internal.getQm(modulation);

    % Get the set of symbol indices based on DM-RS configuration and the
    % number of symbols allocated for PUCCH transmission
    [nset,s1,s2,s3] = getSymbolSets(nPUCCHSym,dmrsSym);

    % Get the number of UCI OFDM symbols in each set
    nUCI = [numel(s1) numel(s2) numel(s3)];

    % Get the total number of OFDM symbols carrying UCI
    nPUCCHSymUCI = sum(nUCI);

    % Get the number of UCI symbols
    if format == 3
        nSymbolUCI = 12*nRBOrSF;
    else
        nSymbolUCI = 12/nRBOrSF;
    end

    % Cumulative sum of number of UCI OFDM symbols in each set
    cSumNUCI = cumsum(nUCI);

    % Find the set index that covers uciPart1
    j = 0;
    for i = 1:nset
        if cSumNUCI(i)*nSymbolUCI*qm >= g1
            j = i;
            break;
        end
    end

    % If the UCI part1 length is more than the total number of RE's
    % allocated
    if ~j
        j = nset;
    end

    % Find nBarSymbolUCI and M
    jMinus1 = j-1; % j minus 1
    if j == 1
        temp = 0;
    else
        temp = cSumNUCI(jMinus1)*nSymbolUCI*qm;
    end
    nBarSymbolUCI = floor((g1-temp)/(nUCI(jMinus1+1)*qm));
    M = mod((g1-temp)/qm,nUCI(jMinus1+1));

    % Get the sets to perform multiplexing
    uciset = sort([s1 s2 s3]);
    sets = {s1,s2,s3};
    setJminus1 = [];
    for i = 1:jMinus1
        setJminus1 = [setJminus1 sets{i}]; %#ok
    end

    % Initialize the intermediate value gBar
    if nBarSymbolUCI > nSymbolUCI
        gBar = zeros(nPUCCHSymUCI,nBarSymbolUCI+(M>0),qm,'like',uciPart1);
    else
        gBar = zeros(nPUCCHSymUCI,nSymbolUCI,qm,'like',uciPart1);
    end

    % Perform multiplexing
    n1 = 0;
    n2 = 0;
    for l = 0:nPUCCHSymUCI-1

        % Get the symbol index from the set of UCI symbols
        symIndex = uciset(l+1);

        % Check for the set which contains the symbol index and account
        % the uci.
        if ismember(symIndex,setJminus1)
            % If symbol index belongs to the union of sets up to j-1, uci
            % part 1 is accounted.
            for k = 0:nSymbolUCI-1
                for v = 0:qm-1
                    if n1 < g1 % To avoid out of bounds indexing
                        gBar(l+1,k+1,v+1) = uciPart1(n1+1);
                    end
                    n1 = n1+1;
                end
            end
        elseif ismember(symIndex,sets{j})
            % If symbol index belongs to set j, both uci part 1 and uci
            % part 2 are accounted.
            if M > 0
                gamma = 1;
            else
                gamma = 0;
            end
            M = M-1;
            for k = 0:nBarSymbolUCI+gamma-1
                for v = 0:qm-1
                    if n1 < g1 % To avoid out of bounds indexing
                        gBar(l+1,k+1,v+1) = uciPart1(n1+1);
                    end
                    n1 = n1+1;
                end
            end
            for k = nBarSymbolUCI+gamma:nSymbolUCI-1
                for v = 0:qm-1
                    if n2 < g2 % To avoid out of bounds indexing
                        gBar(l+1,k+1,v+1) = uciPart2(n2+1);
                    end
                    n2 = n2+1;
                end
            end
        else
            % If symbol index does not belong to union of sets j, uci part
            % 2 is accounted.
            for k = 0:nSymbolUCI-1
                for v = 0:qm-1
                    if n2 < g2 % To avoid out of bounds indexing
                        gBar(l+1,k+1,v+1) = uciPart2(n2+1);
                    end
                    n2 = n2+1;
                end
            end
        end
    end

    % Return multiplexed output
    out = reshape(permute(gBar(1:nPUCCHSymUCI,1:nSymbolUCI,1:qm),[3 2 1]),[],1);

end

function [nset,s1,s2,s3] = getSymbolSets(nPUCCHSym,dmrsSym)
%getSymbolSets provides the sets of symbols according to Table
%   6.3.1.6-1 TS 38.212, with the two inputs, number of symbols allocated
%   for PUCCH and the PUCCH DM-RS symbol indices. The valid combinations are
%   provided below:
%
%       PUCCH duration               PUCCH DM-RS symbol indices
%       --------------               -------------------------
%             4                            [1] or [0 2]
%             5                               [0 3]
%             6                               [1 4]
%             7                               [1 4]
%             8                               [1 5]
%             9                               [1 6]
%            10                          [2 7] or [1 3 6 8]
%            11                          [2 7] or [1 3 6 9]
%            12                          [2 8] or [1 4 7 10]
%            13                          [2 9] or [1 4 7 11]
%            14                         [3 10] or [1 5 8 12]
%
%   The DM-RS symbol indices provided in the table above are relative to the
%   start of PUCCH symbol.

    % Initialize the outputs
    nset = 0;
    s1 = [];
    s2 = [];
    s3 = [];

    % Make the dmrsSym as row
    dmrsSym = reshape(dmrsSym,1,[]);

    % Get the outputs based on valid combinations
    switch nPUCCHSym
        case 4
            if dmrsSym == 1
                nset = 2;
                s1 = [0 2];
                s2 = 3;
            elseif isequal(dmrsSym, [0 2])
                nset = 1;
                s1 = [1 3];
            end
        case 5
            if isequal(dmrsSym, [0 3])
                nset = 1;
                s1 = [1 2 4];
            end
        case 6
            if isequal(dmrsSym, [1 4])
                nset = 1;
                s1 = [0 2 3 5];
            end
        case 7
            if isequal(dmrsSym, [1 4])
                nset = 2;
                s1 = [0 2 3 5];
                s2 = 6;
            end
        case 8
            if isequal(dmrsSym, [1 5])
                nset = 2;
                s1 = [0 2 4 6];
                s2 = [3 7];
            end
        case 9
            if isequal(dmrsSym, [1 6])
                nset = 2;
                s1 = [0 2 5 7];
                s2 = [3 4 8];
            end
        case 10
            if isequal(dmrsSym, [2 7])
                nset = 2;
                s1 = [1 3 6 8];
                s2 = [0 4 5 9];
            elseif isequal(dmrsSym, [1 3 6 8])
                nset = 1;
                s1 = [0 2 4 5 7 9];
            end
        case 11
            if isequal(dmrsSym, [2 7])
                nset = 3;
                s1 = [1 3 6 8];
                s2 = [0 4 5 9];
                s3 = 10;
            elseif isequal(dmrsSym, [1 3 6 9])
                nset = 1;
                s1 = [0 2 4 5 7 8 10];
            end
        case 12
            if isequal(dmrsSym, [2 8])
                nset = 3;
                s1 = [1 3 7 9];
                s2 = [0 4 6 10];
                s3 = [5 11];
            elseif isequal(dmrsSym, [1 4 7 10])
                nset = 1;
                s1 = [0 2 3 5 6 8 9 11];
            end
        case 13
            if isequal(dmrsSym, [2 9])
                nset = 3;
                s1 = [1 3 8 10];
                s2 = [0 4 7 11];
                s3 = [5 6 12];
            elseif isequal(dmrsSym, [1 4 7 11])
                nset = 2;
                s1 = [0 2 3 5 6 8 10 12];
                s2 = 9;
            end
        otherwise % nPUCCHSym equal to 14
            if isequal(dmrsSym, [3 10])
                nset = 3;
                s1 = [2 4 9 11];
                s2 = [1 5 8 12];
                s3 = [0 6 7 13];
            elseif isequal(dmrsSym, [1 5 8 12])
                nset = 2;
                s1 = [0 2 4 6 7 9 11 13];
                s2 = [3 10];
            end
    end

    % Check if there is any set
    if ~nset
        error('Invalid combination of PUCCH duration and PUCCH DM-RS symbol indices specified.');
    end

end
