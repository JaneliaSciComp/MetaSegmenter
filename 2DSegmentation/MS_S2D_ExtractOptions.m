function options = MS_S2D_ExtractOptions(p, fracBlack)
    options.outBW  = p.Results.outBW;
    options.outSeg = p.Results.outSeg;
    options.outRGB = p.Results.outRGB;
    compiled_code = 0;
    if isnumeric(fracBlack)     % original Matlab code
        options.fracBlack =         p.Results.fracBlack;
        options.verbose   = logical(p.Results.verbose);
        options.dispOn    = logical(p.Results.dispOn);
        options.dispOn2   = logical(p.Results.dispOn2);
        options.closeAll  = logical(p.Results.closeAll);
        options.nx        =   int32(p.Results.nx);
        options.ny        =   int32(p.Results.ny);
        options.ix        =   int32(p.Results.ix);
        options.iy        =   int32(p.Results.iy);
        options.maxSize   =   int32(p.Results.maxSize);
    else
        compiled_code     = 1;    % compiled code
        options.fracBlack =         str2double(p.Results.fracBlack);
        options.verbose   = logical(str2double(p.Results.verbose));
        options.dispOn    = logical(str2double(p.Results.dispOn));
        options.dispOn2   = logical(str2double(p.Results.dispOn2));
        options.closeAll  = logical(str2double(p.Results.closeAll));
        options.nx        =   int32(str2double(p.Results.nx));
        options.ny        =   int32(str2double(p.Results.ny));
        options.ix        =   int32(str2double(p.Results.ix));
        options.iy        =   int32(str2double(p.Results.iy));
        options.maxSize   =   int32(str2double(p.Results.maxSize));
    end
    return

