function show_h5stack_labels(h5file, z)
    h5data   = h5read(h5file, '/main');
    u8 = uint8(z); 
    disp(['size(h5data)=' num2str(size(h5data))]);
    disp('. = none');
    disp('+ = all ');
    disp('4 = 1 + 2');
    disp('5 = 1 + 3');
    disp('6 = 2 + 3');

    stack_size = size(h5data, 3);
    len = length(h5file);
%   disp(['size(h5data1)=' num2str(size(h5data(:,:,1))) ' size(squeeze(h5data1))=' num2str(size(squeeze(h5data(:,:,1))))]);
    labels = zeros(size(h5data, 1), size(h5data, 2));
    labels(:,:) = '.';
    labels(h5data(:,:,u8,1)>0 & h5data(:,:,u8,2) >0 & h5data(:,:,u8,3) >0) = '+';
    labels(h5data(:,:,u8,1)>0 & h5data(:,:,u8,2)==0 & h5data(:,:,u8,3)==0) = '1';
    labels(h5data(:,:,u8,2)>0 & h5data(:,:,u8,1)==0 & h5data(:,:,u8,3)==0) = '2';
    labels(h5data(:,:,u8,3)>0 & h5data(:,:,u8,1)==0 & h5data(:,:,u8,2)==0) = '3';
    labels(h5data(:,:,u8,1)>0 & h5data(:,:,u8,2)> 0 & h5data(:,:,u8,3)==0) = '4';
    labels(h5data(:,:,u8,1)>0 & h5data(:,:,u8,2)==0 & h5data(:,:,u8,3)> 0) = '5';
    labels(h5data(:,:,u8,3)>0 & h5data(:,:,u8,1)==0 & h5data(:,:,u8,2)> 0) = '6';
    char(labels)  
%   size(labels,1)
%   for i=1:size(labels,1)
%       class(uint8(i))
%       char(labels(h5data(uint8(i),:)))
%   end
