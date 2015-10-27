function options = MS_S2D_ExtractOptions(p, fracBlack, fracBlack2)
    options.outBW  = p.Results.outBW;
    options.outSeg = p.Results.outSeg;
    options.outRGB     = p.Results.outRGB;
    options.membPr     = p.Results.membPr;
    options.mitoPr     = p.Results.mitoPr;
    options.mitomembPr = p.Results.mitoMembPr;
    options.compiled_code = 0;
    if isnumeric(fracBlack)     % original Matlab code
        options.fracBlack = fracBlack;
        options.fracBlack2= fracBlack2;
        options.verbose   = p.Results.verbose;
        options.hist      = p.Results.hist;
        options.dispOn    = p.Results.dispOn;
        options.dispOn2   = p.Results.dispOn2;
        options.closeAll  = p.Results.closeAll;
        options.noDark    = p.Results.noDark;
        options.padding   = p.Results.padding;
        options.nx        = p.Results.nx;
        options.ny        = p.Results.ny;
        options.dx        = p.Results.dx;
        options.dy        = p.Results.dy;
        options.sx        = p.Results.sx;
        options.sy        = p.Results.sy;
        options.maxSize   = p.Results.maxSize;
        options.resize    = p.Results.resize;
        options.thr       = p.Results.thr;
        options.thr2      = p.Results.thr2;
    else
        options.compiled_code = 1;    % compiled code
        options.fracBlack     =       str2double(fracBlack);
        options.fracBlack2    =       str2double(fracBlack2);
        options.verbose       = int32(str2double(p.Results.verbose));
        options.hist          = int32(str2double(p.Results.hist));
        options.dispOn        = int32(str2double(p.Results.dispOn));
        options.dispOn2       = int32(str2double(p.Results.dispOn2));
        options.closeAll      = int32(str2double(p.Results.closeAll));
        options.noDark        = int32(str2double(p.Results.noDark));
        options.padding       = int32(str2double(p.Results.padding));
        options.nx            = int32(str2double(p.Results.nx));
        options.ny            = int32(str2double(p.Results.ny));
        options.dx            = int32(str2double(p.Results.dx));
        options.dy            = int32(str2double(p.Results.dy));
        options.sx            = int32(str2double(p.Results.sx));
        options.sy            = int32(str2double(p.Results.sy));
        options.maxSize       = int32(str2double(p.Results.maxSize));
        options.resize        = int32(str2double(p.Results.resize)); 
        options.thr           =       str2double(p.Results.thr);
        options.thr2          =       str2double(p.Results.thr2);
    end
    return

