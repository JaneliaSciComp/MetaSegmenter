%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

classdef MS_S2D_InputParser < inputParser
    properties
        defaultVerbose      = 0;
        defaultDispOn       = 0;
        defaultDispOn2      = 0;
        defaultPadding      = 0;
        defaultCloseAll     = 1;
        defaultNoDark       = 0; % don't detect dark structures
        defaultNumSubX      = 1;
        defaultNumSubY      = 1;
        defaultDX           = 50;
        defaultDY           = 50;
        defaultSubX         = 0;   % if > 0, process data only in this subsection
        defaultSubY         = 0;   % if > 0, process data only in this subsection
        defaultThresh       = 0;
        defaultThresh2      = 0;
        defaultOutBW        = '';  % name for output black/white image file
        defaultOutSeg       = '';  % name for output segmentation file
        defaultOutThr       = '';  % name for output intensity thresholds file
        defaultRGB          = 0;   % show RGB labels                       
        defaultHistogram    = 0;
        defaultScale        = 1;
        defaultMembProb     = ''; % path to the membranes probability file
        defaultMitoProb     = ''; % path to the mitochondria probability file
        defaultMitoMembProb = ''; % path to the mitochondria probability file
        defaultUseMembPr    = 0;  % use membrane probabilities instead of grayscale signals
        defaultVesicles     = 0;
        defaultDistType     = 1; % may be 1 or 2
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
                obj.defaultNoDark    = num2str(obj.defaultNoDark);   
                obj.defaultPadding   = num2str(obj.defaultPadding);
                obj.defaultNumSubX   = num2str(obj.defaultNumSubX);     
                obj.defaultNumSubY   = num2str(obj.defaultNumSubY);
                obj.defaultSubX      = num2str(obj.defaultSubX );        
                obj.defaultSubY      = num2str(obj.defaultSubY );        
                obj.defaultThresh    = num2str(obj.defaultThresh);
                obj.defaultThresh2   = num2str(obj.defaultThresh2);
                obj.defaultDX        = num2str(obj.defaultDX);          
                obj.defaultDY        = num2str(obj.defaultDY);                
                obj.defaultScale     = num2str(obj.defaultScale);
                obj.defaultDistType  = num2str(obj.DistType);
            end

            % Optional parameters
            obj.addParamValue('verbose',   obj.defaultVerbose);
            obj.addParamValue('hist',      obj.defaultHistogram);
            obj.addParamValue('dispOn',    obj.defaultDispOn); 
            obj.addParamValue('dispOn2',   obj.defaultDispOn2); 
            obj.addParamValue('closeAll',  obj.defaultCloseAll);
            obj.addParamValue('noDark',    obj.defaultNoDark);
            obj.addParamValue('padding',   obj.defaultPadding);
            obj.addParamValue('nx',        obj.defaultNumSubX);
            obj.addParamValue('ny',        obj.defaultNumSubX);
            obj.addParamValue('dx',        obj.defaultDX);          
            obj.addParamValue('dy',        obj.defaultDY);            
            obj.addParamValue('sx',        obj.defaultSubX);
            obj.addParamValue('sy',        obj.defaultSubY);
            obj.addParamValue('thr',       obj.defaultThresh);
            obj.addParamValue('thr2',      obj.defaultThresh2);
            obj.addParamValue('resize',    obj.defaultScale);  
            obj.addParamValue('outBW',     obj.defaultOutBW')
            obj.addParamValue('outSeg',    obj.defaultOutSeg);
            obj.addParamValue('outThr',    obj.defaultOutThr);
            obj.addParamValue('RGB',       obj.defaultRGB);
            obj.addParamValue('useMembPr', obj.defaultUseMembPr);
            obj.addParamValue('vesicles',  obj.defaultVesicles);
            obj.addParamValue('membPr',    obj.defaultMembProb);
            obj.addParamValue('mitoPr',    obj.defaultMitoProb);
            obj.addParamValue('mitoMembPr',obj.defaultMitoMembProb);
            obj.addParamValue('distType',  obj.defaultDistType);
            obj.KeepUnmatched = true;
            parse@inputParser(obj, fracBlack, fracBlack2, varargin{:});
        end
    end
end
