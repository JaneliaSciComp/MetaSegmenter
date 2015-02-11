function subsections = MS_S2D_DivideImageIntoSubsections(im_size, options)
    subsections(1).ypixels = [1:im_size(1)];
    subsections(1).xpixels = [1:im_size(2)];

    % Determine the max # of sections based on image dimensions and options
    MS = options.maxSize;
    ny_MS = 1+(im_size(1) - mod(im_size(1),MS))/MS;
    nx_MS = 1+(im_size(2) - mod(im_size(2),MS))/MS;
    ny_best = max(ny_MS, options.ny);
    nx_best = max(nx_MS, options.nx);

    % Determine section size
    y_size = round(im_size(1)/ny_best);
    x_size = round(im_size(2)/nx_best);

    % Populate subsections structure with pixels
    k = 0;
    for iy=1:ny_best
        if options.iy > 0 && iy ~= options.iy
            continue;
        end
        for ix=1:nx_best
            if options.ix > 0 && ix ~= options.ix
                continue;
            end
            k = k+1;
            ymin =     y_size*(iy-1) + 1;
            ymax = min(y_size* iy    , im_size(1));
            subsections(k).ypixels = [ymin:ymax];
            xmin =     x_size*(ix-1) + 1;
            xmax = min(x_size* ix    , im_size(2));
            subsections(k).xpixels = [xmin:xmax];
            if options.verbose
                disp(['k=', num2str(k) ' x_size=' num2str(xmax-xmin+1)  ...
                                       ' y_size=' num2str(ymax-ymin+1)]);
            end
        end
    end

