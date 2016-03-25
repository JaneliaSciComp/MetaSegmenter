function MS_S2D_ShowImage(I, image_title, options)
    figure
    if ~options.heatmap
        imshow(I);
        title(image_title);
        impixelinfo;
        colorbar('off');
    else
        colormap('hot');
        imagesc(I);
        colorbar;
    end
    waitforbuttonpress;
    if options.closeAll
        close all;
    end
