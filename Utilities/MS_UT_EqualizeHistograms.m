function match_trn_hist(trndir,inputdir,outputdir)

    allfiles = dir(fullfile(trndir,'*.png'));

    minval=50;
    maxval=200;
    
    nfiles = length(allfiles)
    
    rowcount=1;
    colcount=1;
    bigim=[];
    for i=1:nfiles
        imname = allfiles(i).name
        im=imread(fullfile(trndir,imname));
        
        im(find(im<minval)) = minval;
        im(find(im>maxval)) = maxval;
        
        bigim = [bigim im];
        
    end
    
    bigim = double(bigim-minval);
    bigim = bigim./(maxval-minval);
    bigim = bigim.*255;
    uim = uint8(bigim);
    
    allfiles_in = dir(fullfile(inputdir,'*.png'));

    for i=1:length(allfiles_in)
        imname = allfiles_in(i).name
        im=imread(fullfile(inputdir,imname));
        
        im(find(im<minval)) = minval;
        im(find(im>maxval)) = maxval;
        
        eqim = imhistmatch(im, uim,256);
        
        savename = fullfile(fullfile(outputdir,imname));
        imwrite(eqim, savename,'png');
    end
