%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

classdef MS_S2D_InputParser < inputParser
    properties
        defaultVerbose   = 0;
        defaultDispOn    = 0;
        defaultDispOn2   = 0;
        defaultPadding   = 0;
        defaultCloseAll  = 1;
        defaultNumSubX   = 1;
        defaultNumSubY   = 1;
        defaultDX        = 50;
        defaultDY        = 50;
        defaultSubX      = 0;
        defaultSubY      = 0;
        defaultMaxSize   = Inf; 
        defaultOutBW     = ''; % name for output black/white image file
        defaultOutSeg    = ''; % name for output segmentation file
        defaultOutRGB    = ''; % name for output colored labels file
        defaultHistogram = 0;
    end

    methods
        function obj = MS_S2D_InputParser % class constructor 
            obj  = obj@inputParser;     
        end
        function parse(obj, fracBlack, fracBlack2, varargin)
            compiled_code = 0;
            % If compiled code
            if ~isnan(str2double(fracBlack))
                compiled_code = 1;
                fracBlack  = str2double(fracBlack);
                fracBlack2 = str2double(fracBlack2);
            end

            % Required inputs
            obj.addRequired('fracBlack',  @isnumeric);
            obj.addRequired('fracBlack2', @isnumeric);
            % If varargin{1} is a string convertable to number (i.e. this is compiled code)
            if size(varargin,1) > 0 && ~isnan(str2double(varargin{1}))
                varargin{1} = int32(str2double(varargin{1}));
            end

            % Optional parameters
            if compiled_code
                obj.defaultVerbose   = num2str(obj.defaultVerbose);
                obj.defaultHistogram = num2str(obj.defaultHistogram);
                obj.defaultDispOn    = num2str(obj.defaultDispOn);        
                obj.defaultDispOn2   = num2str(obj.defaultDispOn2);
                obj.defaultCloseAll  = num2str(obj.defaultCloseAll);
                obj.defaultPadding   = num2str(obj.defaultPadding);
                obj.defaultNumSubX   = num2str(obj.defaultNumSubX);     
                obj.defaultNumSubY   = num2str(obj.defaultNumSubY);
                obj.defaultSubX      = num2str(obj.defaultSubX );        
                obj.defaultSubY      = num2str(obj.defaultSubY );        
                obj.defaultSectionX  = num2str(obj.defaultSectionX);
                obj.defaultSectionY  = num2str(obj.defaultSectionY);
                obj.defaultMaxSize   = num2str(obj.defaultMaxSize);
            end

            % Optional parameters
            obj.addParamValue('verbose', obj.defaultVerbose);
            obj.addParamValue('hist',    obj.defaultHistogram);
            obj.addParamValue('dispOn',  obj.defaultDispOn); 
            obj.addParamValue('dispOn2', obj.defaultDispOn2); 
            obj.addParamValue('closeAll',obj.defaultCloseAll);
            obj.addParamValue('padding', obj.defaultPadding);
            obj.addParamValue('nx',      obj.defaultNumSubX);
            obj.addParamValue('ny',      obj.defaultNumSubX);
            obj.addParamValue('dx',      obj.defaultDX);          
            obj.addParamValue('dy',      obj.defaultDY);            
            obj.addParamValue('sx',      obj.defaultSubX);
            obj.addParamValue('sy',      obj.defaultSubY);
            obj.addParamValue('maxSize', obj.defaultMaxSize);
            obj.addParamValue('outBW',   obj.defaultOutBW')
            obj.addParamValue('outSeg',  obj.defaultOutSeg);
            obj.addParamValue('outRGB',  obj.defaultOutRGB);
            obj.KeepUnmatched = true;
            parse@inputParser(obj, fracBlack, fracBlack2, varargin{:});
        end
    end
end
