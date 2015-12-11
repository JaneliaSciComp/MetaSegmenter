function subsections = MS_S2D_DivideImageIntoSubsections(im_size, options)
    subsections(1).ypixels = [1:im_size(1)];
    subsections(1).xpixels = [1:im_size(2)];

    % Determine section size
    y_size = round(im_size(1)/options.ny);
    x_size = round(im_size(2)/options.nx);

    % Populate subsections structure with pixels
    if options.verbose
        disp(['im_size=' num2str(im_size)]);
    end
    k = 0; % subsection index
    for iy=1:options.ny
        if options.sy > 0 && iy ~= options.sy
            continue;
        end
        for ix=1:options.nx
            if options.sx > 0 && ix ~= options.sx
                continue;
            end
            k = k+1;
            ymin =     y_size*(iy-1) + 1;
            if iy < options.ny
                ymax = y_size* iy;                         
            else
                ymax = im_size(1);
            end
            subsections(k).ypixels = [ymin:ymax];
            xmin =     x_size*(ix-1) + 1;
            if ix < options.nx
                xmax = x_size* ix;
            else
                xmax = im_size(2);
            end
            subsections(k).xpixels = [xmin:xmax];
            disp(['subswection k=' num2str(k) ' ypixels_min=' num2str(min(subsections(k).ypixels)) ' ypixels_max=' num2str(max(subsections(k).ypixels))]);
        end
    end

