%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

classdef MS_S2D_InputParser < inputParser
    properties
        defaultVerbose      = 0;
        defaultHeatMap      = 0;
        defaultDebug        = 0;
        defaultDispOn       = 0;
        defaultDispOn2      = 0;
        defaultPadding      = 0;
        defaultCloseAll     = 1;
        defaultDarkStr      = 0; % don't detect dark structures
        defaultNumSubX      = 0;
        defaultNumSubY      = 0;
        defaultDX           = 50;
        defaultDY           = 50;
        defaultMaxRegSize   = 1200; % 700 needed by 00001; needs to be larger to handle glial cells 
        defaultMinRegSize   = 200; % mininal size of a segmentation region (# pixels). Smaller regions will be discarded
        defaultMaxWinSize   = 250;
        defaultMinWinSize   = 64;  % mininal size of a window for recursive thresholding; iterate while dx,dy >= minWin
        defaultMinMitoFrac  = 0.1; % region below this frac is non-mito
        defaultMaxMitoFrac  = 0.4; % region above this frac is the true mito
        defaultNumRefPasses = 1;   % number of passes / cycles of the refinement procedure
        defaultIncRef       = 0;   % perform only incrementation refinement; skip decrementation refinement
        defaultSubX         = 0;   % if > 0, process data only in this subsection
        defaultSubY         = 0;   % if > 0, process data only in this subsection
        defaultSubwindow    = '';  % if not '', process only subwindow 'xmin,ymin,xmax,ymax'
        defaultThrMembPr    = 0.2;
        defaultThresh       = 0;
        defaultThresh2      = 0;
        defaultOutBW        = '';  % name for the output black/white image file
        defaultOutSeg       = '';  % name for the output segmentation file
        defaultOutRGB       = '';  % name for the output RGB file
        defaultOutThr       = '';  % name for the output intensity thresholds file
        defaultRGB          = 0;   % display the RGB labels                       
        defaultHistogram    = 0;
        defaultShowIter     = 0;
        defaultScale        = 1;
        defaultMembProb     = ''; % path to the membranes probability file
        defaultMitoProb     = ''; % path to the mitochondria probability file
        defaultMitoMembProb = ''; % path to the mitochondria probability file
        defaultUseMembPr    = 0;  % use membrane probabilities instead of grayscale signals
        defaultVesicles     = 0;
        defaultTsallisq     = 0;
        defaultRenyia       = 0;
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
                obj.defaultVerbose     = num2str(obj.defaultVerbose);
                obj.defaultHeatMap     = num2str(obj.defaultHeatMap);
                obj.defaultDebug       = num2str(obj.defaultDebug);
                obj.defaultHistogram   = num2str(obj.defaultHistogram);
                obj.defaultShowIter    = num2str(obj.defaultShowIter);
                obj.defaultDispOn      = num2str(obj.defaultDispOn);        
                obj.defaultDispOn2     = num2str(obj.defaultDispOn2);
                obj.defaultCloseAll    = num2str(obj.defaultCloseAll);
                obj.defaultDarkStr     = num2str(obj.defaultDarkStr);   
                obj.defaultPadding     = num2str(obj.defaultPadding);
                obj.defaultNumSubX     = num2str(obj.defaultNumSubX);     
                obj.defaultNumSubY     = num2str(obj.defaultNumSubY);
                obj.defaultSubX        = num2str(obj.defaultSubX );        
                obj.defaultSubY        = num2str(obj.defaultSubY );        
                obj.defaultThrMembPr   = num2str(obj.defaultThrMembPr);
                obj.defaultThresh      = num2str(obj.defaultThresh);
                obj.defaultThresh2     = num2str(obj.defaultThresh2);
                obj.defaultDX          = num2str(obj.defaultDX);          
                obj.defaultDY          = num2str(obj.defaultDY);         
                obj.defaultMaxRegSize  = num2str(obj.defaultMaxRegSize);       
                obj.defaultMinRegSize  = num2str(obj.defaultMinRegSize);
                obj.defaultMaxWinSize  = num2str(obj.defaultMaxWinSize);
                obj.defaultMinWinSize  = num2str(obj.defaultMinWinSize);
                obj.defaultNumRefPasses= num2str(obj.defaultNumRefPasses);
                obj.defaultIncRef      = num2str(obj.defaultIncRef);
                obj.defaultMinMitoFrac = num2str(obj.defaultMinMitoFrac);
                obj.defaultMaxMitoFrac = num2str(obj.defaultMaxMitoFrac);
                obj.defaultScale       = num2str(obj.defaultScale);
                obj.defaultTsallisq    = num2str(obj.defaultTsallisq);
                obj.defaultRenyia      = num2str(defaultRenyia);
            end

            % Optional parameters
            obj.addParamValue('verbose',   obj.defaultVerbose);
            obj.addParamValue('heatmap',   obj.defaultHeatMap);
            obj.addParamValue('debug',     obj.defaultDebug);
            obj.addParamValue('hist',      obj.defaultHistogram);
            obj.addParamValue('showIter',  obj.defaultShowIter);
            obj.addParamValue('dispOn',    obj.defaultDispOn); 
            obj.addParamValue('dispOn2',   obj.defaultDispOn2); 
            obj.addParamValue('closeAll',  obj.defaultCloseAll);
            obj.addParamValue('darkStr',   obj.defaultDarkStr);
            obj.addParamValue('padding',   obj.defaultPadding);
            obj.addParamValue('nx',        obj.defaultNumSubX);
            obj.addParamValue('ny',        obj.defaultNumSubY);
            obj.addParamValue('dx',        obj.defaultDX);          
            obj.addParamValue('dy',        obj.defaultDY);       
            obj.addParamValue('maxSize',   obj.defaultMaxRegSize); 
            obj.addParamValue('minSize',   obj.defaultMinRegSize); 
            obj.addParamValue('maxWin',    obj.defaultMaxWinSize);
            obj.addParamValue('minWin',    obj.defaultMinWinSize);   
            obj.addParamValue('minMitoFr', obj.defaultMinMitoFrac);
            obj.addParamValue('maxMitoFr', obj.defaultMaxMitoFrac);
            obj.addParamValue('numRef',    obj.defaultNumRefPasses);
            obj.addParamValue('incRef',    obj.defaultIncRef);
            obj.addParamValue('sx',        obj.defaultSubX);
            obj.addParamValue('sy',        obj.defaultSubY);
            obj.addParamValue('thrMembPr', obj.defaultThrMembPr);
            obj.addParamValue('thr',       obj.defaultThresh);
            obj.addParamValue('thr2',      obj.defaultThresh2);
            obj.addParamValue('resize',    obj.defaultScale);  
            obj.addParamValue('outBW',     obj.defaultOutBW)
            obj.addParamValue('outSeg',    obj.defaultOutSeg);
            obj.addParamValue('outRGB',    obj.defaultOutRGB);
            obj.addParamValue('outThr',    obj.defaultOutThr);
            obj.addParamValue('RGB',       obj.defaultRGB);
            obj.addParamValue('useMembPr', obj.defaultUseMembPr);
            obj.addParamValue('vesicles',  obj.defaultVesicles);
            obj.addParamValue('membPr',    obj.defaultMembProb);
            obj.addParamValue('subWin',    obj.defaultSubwindow);
            obj.addParamValue('mitoPr',    obj.defaultMitoProb);
            obj.addParamValue('mitoMembPr',obj.defaultMitoMembProb);
            obj.addParamValue('q',         obj.defaultTsallisq);
            obj.addParamValue('a',         obj.defaultRenyia);
            obj.KeepUnmatched = true;
            parse@inputParser(obj, fracBlack, fracBlack2, varargin{:});
        end
    end
end
