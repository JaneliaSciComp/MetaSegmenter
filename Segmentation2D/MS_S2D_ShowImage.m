function MS_S2D_ShowImage(I, image_title, options)
    figure
    imshow(I);
    title(image_title);

    
    
    waitforbuttonpress;
    if options.closeAll
        close all;
    end
