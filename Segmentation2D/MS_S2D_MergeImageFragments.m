%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%
% dir = direction of merge ('x' or 'y')
% nf  = total number of fragments in the direction of merge
% df  = pixels overlap in the direction of merge
function [] = MS_S2D_MergeImageFragments(dir, nf0, df0, verbose0, ...
                                         Ysize0, Xsize0, output_name, varargin)
    disp(['nargin=' num2str(nargin)]);
    num_input_files = nargin - 7;
    disp(['class(df0)=' class(df0) ' class(nf0)='  class(nf0)]);
    verbose = int32(str2double(verbose0));
    df      = int32(str2double(df0));
    nf      = int32(str2double(nf0));
    disp(['dir=' dir ' df=' df ' num_input_files=' num2str(num_input_files) ...
          ' output_name=' output_name]);
    for i=1:num_input_files
        disp(['i=' num2str(i) ' input_file=' varargin{i}]);
    end
    len = length(output_name);
    xsize = 0; % x dimension of one fragment
    ysize = 0; % y dimension of one fragment
    Xsize = 0; % x dimension of the final output image
    Ysize = 0; % y dimension of the final output image
    mx    = 0;
    my    = 0;
    M0    = [];
    M1    = [];
    for i=1:num_input_files
        disp(['...loading '   varargin{i}]);
        ind_y = get_index(varargin{i}, '_y');
        ind_z = get_index(varargin{i}, '_z');
        if     strcmp(dir, 'x')
            ind_x = get_index(varargin{i}, '_x');
            if verbose
                disp(['ind_x=' num2str(ind_x) ' ind_y=' num2str(ind_y) ...
                     ' ind_z=' num2str(ind_z)]);
            end
        else
            if verbose
                disp(['ind_y=' num2str(ind_y) ' ind_z=' num2str(ind_z)]); 
            end
        end
        im = imread(varargin{i});
        disp(['i=' num2str(i) ' size(im)=' num2str(size(im))]);
        if i == 1
            if strcmp(dir, 'x')
                dir
                nf
                size(im, 2)
                df
                Xsize = int32(str2double(Xsize0));
                ysize = size(im, 1);
                Ysize = ysize;
                disp(['Xsize=' num2str(Xsize) ' Ysize=' num2str(Ysize)]);
                output_image = uint8(zeros(Ysize, Xsize));
                disp(['size(output_image)=' num2str(size(output_image))]);
                if num_input_files == 1
                    output_image = im
                else % num_input_files > 1
                    xsize = size(im, 2) - df;
                    xsize
                    mx    = 2*df;
                    my    = ysize;
                    output_image(1:ysize, 1:(xsize+df)) = im;
                    R     = uint8(zeros(my, mx));
                    M1    = uint8(zeros(my, mx));
                    Ones  = uint8(ones( my, mx));
                    for j=1:mx
                        R(:,j) = randi(mx, my, 1); % column of size my containing random integers 1:mx
                    end
                    for j=2:mx
                        Rj = R(:,j);
                        Mj = M1(:,j);
                        Mj(Rj >= j) = 1;
                        M1(:,j) = Mj;                 
                    end
                    M0 = Ones - M1;
                end
            else % dir == y
                dir
                nf
                size(im, 1)
                df
                xsize = size(im, 2);
                Xsize = xsize;
                Ysize = int32(str2double(Ysize0));
                output_image = uint8(zeros(Ysize, Xsize));
                if num_input_files == 1
                    output_image = im;
                else
                    ysize = size(im, 1) - df;
                    ysize
                    output_image(1:(ysize+df), :) = im;
                    mx    = xsize;
                    my    = 2*df;
                    R     = uint8(zeros(my, mx));
                    M1    = uint8(zeros(my, mx));
                    Ones  = uint8(ones( my, mx));
                    for j=1:my
                        R(j,:) = randi(my, mx, 1)'; % row of size mx containing random integers 1:my
                    end
                    for j=2:my
                        Rj = R(j,:);
                        Mj = M1(j,:);
                        Mj(Rj >= j) = 1;
                        M1(j,:) = Mj;
                    end
                    M0 = Ones - M1;
                end
            end
        else % Merge current fragment with the output image
            if strcmp(dir, 'x')
                xdim = size(output_image, 2);
                xmid = (i-1)*xsize;
                xmin =       xmid  - df + 1; % left  boundary of the overlap zone
                xmax =       xmid  + df;     % right boundary of the overlap zone
                disp(['i=' num2str(i) ' xsize=' num2str(xsize) ' xmin=' num2str(xmin) ' xmax=' num2str(xmax) ' dx=' num2str(xmax-xmin)]);
                disp(['im_min=' num2str(2*df+1) ' im_max=' num2str(xsize+  df) ' d_im=' num2str(xsize+  df - 2*df -1)]);
                if i < num_input_files
                    output_image(:,(xmax + 1):(i*xsize+df)) = im(:, (2*df+1):(xsize+2*df));
                else
                    im_xsize = size(im, 2);
                    output_image(:,(xmax + 1):Xsize       ) = im(:, (2*df+1):im_xsize);    
                end
                output_image(:, xmin:xmax) = ...
                    round(M0 .* output_image(:, xmin:xmax)...
                        + M1 .* uint8(im(:, 1:(2*df))));
            else % dir == 'y'
                ydim = size(output_image, 1);
                ymid = (i-1)*ysize;
                ymin =       ymid  - df + 1;
                ymax =       ymid  + df;
                if i < num_input_files
                    disp(['i=' num2str(i) ' size(output_image)=' num2str(size(output_image)) ' i*ysize+df=' num2str(i*ysize+df)]);
                    disp(['i=' num2str(i) ' size(im)='  num2str(size(im)) ' ysize+2*df=' num2str(ysize+2*df)]);
                    output_image((ymax + 1):(i*ysize+df),:) = im((2*df+1):(ysize+2*df),:);
                else
                    im_ysize = size(im, 1);
                    output_image((ymax + 1):Ysize       ,:) = im((2*df+1):im_ysize,:);
                end
                disp(['i=' num2str(i) ' ymin=' num2str(ymin) ' ymax=' num2str(i*ysize+df)]);
                output_image( ymin:ymax,:) = ...
                    round(M0 .* output_image(ymin:ymax,:)...
                        + M1 .* uint8(im(1:(2*df),:)));
            end
        end
    end
    % Output the merged BW image (*_BW.png)
    imwrite(logical(output_image), output_name);
    output_image_bw = im2bw(mat2gray(output_image), 0.5);

    if strcmp(dir, 'y') 
        % Produce and output labels image
        disp(['size(output_image_bw)=' num2str(size(output_image_bw))]);
        seg_output_name = [ output_name(1:numel(output_name)-6) 'Seg.h5'];
        disp(' ');
        disp(['...generating file ' seg_output_name ]);
        L = MS_S2D_GenerateLabelsMatrix(output_image_bw, verbose);                  
        disp(['size(L)=' num2str(size(L)) ' class(L)='  num2str(class(L(1,1))) ' max(L)=' num2str(max(L(:)))]);
        hdf5write(seg_output_name,'/main',uint64(L));

        produce_RGB = 0;
        if produce_RGB > 0
            % Produce and output RGB segmentattion image
            rgb_output_name = [ output_name(1:numel(output_name)-6) 'RGB.png'];
            disp(' ');
            disp(['...generating file ' rgb_output_name ]);
            Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
            disp('');
            disp(['size(Lrgb)=' num2str(size(Lrgb)) ' class(Lrgb)='  class(Lrgb(1,1,1))]);
            imwrite(Lrgb, rgb_output_name);
        end
    end

% ------------------------------------------------------------------------------

function [ ind ] = get_index(mystring, pattern)
    if ~isempty(strfind(mystring, pattern))
        split1 = strsplit(mystring, pattern);
        disp(['pattern=' pattern ' split1=' split1]);
        split2 = strsplit(char(split1(2)), '_');
        ind    = int32(str2double(split2(1)));
    else
        ind = -1;
    end
    disp(['ind=' num2str(ind)]);
