%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

classdef MS_S2D_InputParser < inputParser
    properties
        defaultFracBlack = 0.;
        defaultVerbose   = 0;
        defaultDispOn    = 0;
        defaultDispOn2   = 0;
        defaultCloseAll  = 1;
        defaultNumTilesX = 1;
        defaultNumTilesY = 1;
        defaultTileX     = 1;
        defaultTileY     = 1;
        defaultOutputName= '';
    end

    methods
        function obj = MS_S2D_InputParser % class constructor 
            obj  = obj@inputParser;     
        end
        function parse(obj, fracBlack, varargin)
             compiled_code = 0;
            % If compiled code
            if ~isnan(str2double(fracBlack))
                compiled_code = 1;
                fracBlack = str2double(fracBlack);
            end

            % Required inputs
            obj.addRequired('fracBlack', @isnumeric);
            % If varargin{1} is a string convertable to number (i.e. this is compiled code)
            if size(varargin,1) > 0 && ~isnan(str2double(varargin{1}))
                varargin{1} = int32(str2double(varargin{1}));
            end

            % Optional parameters
            if compiled_code
                obj.defaultVerbose   = num2str(obj.defaultVerbose);
                obj.defaultDispOn    = num2str(obj.defaultDispOn);        
                obj.defaultDispOn2   = num2str(obj.defaultDispOn2);
                obj.defaultCloseAll  = num2str(obj.defaultCloseAll);
                obj.defaultNumTilesX = num2str(obj.defaultNumTilesX);
                obj.defaultNumTilesY = num2str(obj.defaultNumTilesY);
                obj.defaultTileX     = num2str(obj.defaultTileX);        
                obj.defaultTileY     = num2str(obj.defaultTileY);        
            end

            % Optional parameters
            obj.addParamValue('verbose', obj.defaultVerbose);
            obj.addParamValue('dispOn',  obj.defaultDispOn); 
            obj.addParamValue('dispOn2', obj.defaultDispOn2); 
            obj.addParamValue('closeAll',obj.defaultCloseAll);
            obj.addParamValue('nx',      obj.defaultNumTilesX);
            obj.addParamValue('ny',      obj.defaultNumTilesY);
            obj.addParamValue('ix',      obj.defaultTileX);          
            obj.addParamValue('iy',      obj.defaultTileY);
            obj.addParamValue('outName', obj.defaultOutputName);
            obj.KeepUnmatched = true;

            parse@inputParser(obj, fracBlack, varargin{:});
        end
    end
end
