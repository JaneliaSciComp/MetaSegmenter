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
    if options.verbose
        disp(['im_size=' num2str(im_size)]);
    end
    k = 0;
    for iy=1:ny_best
        if options.sy > 0 && iy ~= options.sy
            continue;
        end
        for ix=1:nx_best
            if options.sx > 0 && ix ~= options.sx
                continue;
            end
            k = k+1;
            ymin =     y_size*(iy-1) + 1;
            if iy < ny_best
                ymax = y_size* iy;                         
            else
                ymax = im_size(1);
            end
            subsections(k).ypixels = [ymin:ymax];
            xmin =     x_size*(ix-1) + 1;
            if ix < nx_best
                xmax = x_size* ix;
            else
                xmax = im_size(2);
            end
            subsections(k).xpixels = [xmin:xmax];
            if options.verbose
                disp(['subsection=', num2str(k) ' x_size=' num2str(xmax-xmin+1)  ...
                                                ' y_size=' num2str(ymax-ymin+1)]);
            end
        end
    end

