%hVectorDataSource Create data vector source of given input vector or PN source
%
%   VECSOURCE = hVectorDataSource(...) constructs a vector data
%   source object VECSOURCE, which then can be used to generate data.
% 
%   VECSOURCE = hVectorDataSource(DATASOURCE) constructs a vector data
%   source object VECSOURCE, from a scalar, vector, character vector or
%   cell array specified by DATASOURCE. If DATASOURCE is scalar or vector,
%   the values are looped around to generate data. DATASOURCE can also be
%   one of the pseudo-random sequences specified by a value in the set
%   ('PN9-ITU', 'PN9', 'PN11', 'PN15', 'PN23') or can be a cell array
%   containing the pseudo-random sequence source definition and its seed in
%   the format {PNSOURCE, SEED}, where PNSOURCE can be one of the above
%   mentioned values. If no seed is specified for the PN sources, the shift
%   register will be initialized with all ones. If seed is specified as 0,
%   the shift register will be random initialized.
%
%   VECSOURCE = hVectorDataSource(DATASOURCE,SEED) constructs a
%   pseudo-random sequence source, where DATASOURCE can take one of the
%   values ('PN9-ITU', 'PN9', 'PN11', 'PN15', 'PN23') and SEED is the
%   shift register initial state.
% 
%   VECTOR = VECSOURCE.getPacket(LENGTH) creates a data vector VECTOR by
%   appropriate looping of the vector specified in DATASOURCE or as per the
%   specified PN sequence for given LENGTH.
%
%   VECSOURCE.reset() resets the internal state to point to the start of
%   the input sequence.
%
%   Examples:
%   % Example 1: 
%   % Create a data source from the sequence [1 0 1]
%   source = hVectorDataSource([1 0 1]);
%   % Generate data of length 10
%   data = source.getPacket(10)' 
%   
%   % The above example returns:   
%   data = [1 0 1 1 0 1 1 0 1 1] 
% 
%   % Example 2: 
%   % Create a data source from standard PN sequence 'PN9-ITU' and seed 2
%   source = hVectorDataSource({'PN9-ITU',2});
%   % Generate data of length 10
%   data = source.getPacket(10)' 
%   
%   % The above example returns:   
%   data = [0 1 0 0 0 0 0 0 0 0] 

%   Copyright 2007-2018 The MathWorks, Inc.

function obj = hVectorDataSource(sourcespec,initialstate)
    
    defaultseed = 1; % Flag to indicate default PN seed
    if nargin == 0
        % Zeros as default
        sourcespec = 0;
    end
    
    if isempty(sourcespec)
       error('DataSource cannot be empty.') 
    end
    
    if nargin < 2
        % Zeros as default
        initialstate = 0;
    else
        % seed specified by user
        defaultseed = 0;
    end
        
    % If the source is PN, sourcespec can be a character vector or a cell
    % array {'PN9',seed} or {'PN9'} or 'PN9'
    if iscell(sourcespec)
        % {'PN9',seed} or {'PN9'} 
        if numel(sourcespec) == 2
            initialstate = sourcespec{2};
            defaultseed = 0;
        else
            % If {'PN9'}, use default seed
            defaultseed = 1;
        end
        sourcespec = sourcespec{1};
        % Now source can only be PN source
        if ~(ischar(sourcespec) || isstring(sourcespec))
            error('The DataSource cell array must be of the format {pnsource,seed} or {pnsource}');
        end
    end
    
    if ischar(sourcespec) || isstring(sourcespec)
        switch upper(sourcespec)
            case 'PN9-ITU'
                poly = [9 4 0];
                shiftreglength = 9;
            case 'PN9'
                poly = [9 5 0];
                shiftreglength = 9;
            case 'PN11'
                poly = [11 2 0];
                shiftreglength = 11;
            case 'PN15'
                poly = [15 14 0];
                shiftreglength = 15;
            case 'PN23'
                poly = [23 5 0];
                shiftreglength = 23;
            otherwise
                error('The allowed DataSource strings are ''PN9-ITU'',''PN9'',''PN11'',''PN15'',''PN23''');
        end
                
        % Set the initial seed to be all ones if not specified 
        if defaultseed
            initialstate = 2^shiftreglength-1;
        else
            % If user-specified seed is 0, randomize the seed
            if initialstate==0
                initialstate = randi([1, 2^shiftreglength-1],1);
            end
        end
        
        % Validate initial state - can only be positive integer or 0
        if ~(isnumeric(initialstate) && isscalar(initialstate) ...
                && ((initialstate >=0) && (initialstate <= 2^shiftreglength-1))) % Range validation
           error('For %s, the seed must be an integer between 0 to %d',sourcespec,2^shiftreglength-1);
        end
    end
    
    currentstate = initialstate;
  
    % Public definition of class methods    
    if ischar(sourcespec) || isstring(sourcespec)
        obj.getPacket = @getPacketPN;
        source = comm.PNSequence('Polynomial',poly,'InitialConditions',dec2bin(initialstate,shiftreglength)-'0','Mask',[zeros(1,shiftreglength-1) 1]);
    elseif isnumeric(sourcespec)
        if currentstate >= length(sourcespec)
            error('Invalid seed specified for vector input, must be 0 to %d', length(sourcespec)-1);
        end
        obj.getPacket = @getPacket;
    else
        error('Invalid DataSource, must be vector or one of the PN strings (''PN9-ITU'',''PN9'',''PN11'',''PN15'',''PN23'')');
    end
    obj.reset = @reset;
    
    % Class 'Method' implementations
    function bitsout = getPacket(psize)
        
          bitsout = zeros(psize,1);     
          for i = 1:psize
             bitsout(i) = sourcespec(currentstate+1);
             currentstate = mod(currentstate+1,length(sourcespec));
          end  
    end

    function bitsout = getPacketPN(psize)
        if isLocked(source)
            if source.SamplesPerFrame ~= psize
                release(source);
                source.SamplesPerFrame = psize;
            end
            bitsout = source();
        else
            source.SamplesPerFrame  = psize;
            bitsout = source();
        end
    end

    function reset()
        currentstate = initialstate;
    end

end
