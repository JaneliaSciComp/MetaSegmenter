function show_h5stack(h5file, id)
    h5data   = h5read(h5file, '/main');
    class(h5data);
    
    disp(['size(h5data)=' num2str(size(h5data))]);
    stack_size = size(h5data, 3);
    len = length(h5file);
    if numel(size(h5data)) == 3
        for i=1:size(h5data,1)
            h5data_uint16(i,:,id)=typecast(h5data(i,:,id), 'uint16');
        end
        
        imshow(mat2gray(squeeze(h5data(:,:,id))), []); % scale automatically
%       imagesc(squeeze(h5data(:,:,id)));
    else
        layer = double(squeeze(h5data(:,:,id,:)));
        disp(['size(layer)=' num2str(size(layer))]);
        imshow(layer);
%       imagesc(squeeze(h5data(:,:,id,3)));
    end
