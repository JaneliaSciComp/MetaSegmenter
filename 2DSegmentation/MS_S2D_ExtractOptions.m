function options = MS_S2D_ExtractOptions(p, fracBlack, fracBlack2)
    options.outBW  = p.Results.outBW;
    options.outSeg = p.Results.outSeg;
    options.outRGB = p.Results.outRGB;
    options.compiled_code = 0;
    if isnumeric(fracBlack)     % original Matlab code
        options.fracBlack = fracBlack;
        options.fracBlack2= fracBlack2;
        options.verbose   = p.Results.verbose;
        options.dispOn    = p.Results.dispOn;
        options.dispOn2   = p.Results.dispOn2;
        options.closeAll  = p.Results.closeAll;
        options.nx        = p.Results.nx;
        options.ny        = p.Results.ny;
        options.ix        = p.Results.ix;
        options.iy        = p.Results.iy;
        options.maxSize   = p.Results.maxSize;
    else
        options.compiled_code = 1;    % compiled code
        options.fracBlack     =       str2double(fracBlack);
        options.fracBlack2    =       str2double(fracBlack2);
        options.verbose       = int32(str2double(p.Results.verbose));
        options.dispOn        = int32(str2double(p.Results.dispOn));
        options.dispOn2       = int32(str2double(p.Results.dispOn2));
        options.closeAll      = int32(str2double(p.Results.closeAll));
        options.nx            = int32(str2double(p.Results.nx));
        options.ny            = int32(str2double(p.Results.ny));
        options.ix            = int32(str2double(p.Results.ix));
        options.iy            = int32(str2double(p.Results.iy));
        options.maxSize       = int32(str2double(p.Results.maxSize));
    end
    return

