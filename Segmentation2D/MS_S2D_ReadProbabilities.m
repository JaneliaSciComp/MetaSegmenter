function probs = MS_S2D_ReadProbabilities(infile, id)
    % Read probabilities from either HDF5 or image file
    probs = [];
    if numel(findstr('.h5', infile)) > 0 
        try
            probs  = h5read(infile, '/probabilities')';
        catch  
            % NeuroProof
            probs4 = h5read(infile, '/stack');
            probs  = squeeze(probs4(id,:,:))';
            disp(['size(probs)=' num2str(size(probs))]);
        end
    else
        probs = im2double(imread(infile))';
    end
    disp(['size(probs=' num2str(size(probs))]);
