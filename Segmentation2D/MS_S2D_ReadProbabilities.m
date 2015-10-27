function probs = MS_S2D_ReadProbabilities(infile)
    % Read probabilities from either HDF5 or image file
    probs = [];
    if numel(findstr('.h5', infile)) > 0 
        probs  = h5read(infile, '/probabilities')';
    else
        probs = im2double(imread(infile))';
    end
    disp(['size(probs=' num2str(size(probs))]);
