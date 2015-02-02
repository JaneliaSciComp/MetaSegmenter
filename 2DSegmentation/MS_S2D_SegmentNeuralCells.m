%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

function MS_S2D_SegmentNeurons2D(varargin)   
warning off;
%%
clc;
try
    image_file = varargin{1};
catch
    output_usage_message();
    return
end

frac_black = 0.45;
if nargin > 1
    frac_black = double(varargin{2});
end 

% Original image    
I = imread(image_file); I = mat2gray(I(:,:,1));
disp(['min(I)=' num2str(min(min(I))) ' max(I)=' num2str(max(max(I))) ' mean2(I)=' num2str(mean2(I)) ' std2(I)=' num2str(std2(I))]);
In = round(double(I)*255.0/max(max(double(I))));
disp(['min(In)=' num2str(min(min(In))) ' max(In)=' num2str(max(max(In))) ' mean2(In)=' num2str(mean2(In)) ' std2(In)=' num2str(std2(In))]);
imshow(I);
title('Original image (I)');
waitforbuttonpress;
close all;

% Black-white image
num_bin_locs = 1001;
[N,bin_locs] = imhist(I,num_bin_locs);
Nsum = sum(N);
disp(['numel(N)=' num2str(numel(N)) ' Nsum=' num2str(Nsum) ' numel(I)=' num2str(numel(I)) ' numel(bin_locs)=' num2str(numel(bin_locs))]);
cum_hist = get_cumulative_histogram(N, Nsum);
%figure
%plot(cum_hist);
%drawnow;
%waitforbuttonpress;
%close all;

threshold = get_threshold_intensity(N, frac_black, bin_locs);
disp(['final threshold_intensity=' num2str(threshold)]);
disp(['size(I)=' num2str(size(I)) ' class(I)=' class(I(1,1))]);
Ibw = my_im2bw(I, threshold);  
disp(['size(Ibw)=' num2str(size(Ibw)) ' class(Ibw)=' class(Ibw(1,1))]);
Ibw = bwareaopen(Ibw,20); 
disp(['my frac_black=' num2str(sum(sum(Ibw == 0))/numel(Ibw))]);
figure
imshow(Ibw), title('Black-white image (Ibw)')
drawnow;
waitforbuttonpress;
close all;

% Fill holes and open gaps
Ibwf = imfill(Ibw, 'holes');
Ibwf = bwareaopen(Ibwf,20);
%Ibwfr = imcomplement(Ibwf);
%Ibwfr = imclose(Ibwfr, strel('disk',1));
%Ibwf = imcomplement(Ibwfr);
%Ibwfo = imopen(Ibwf, strel('disk',2));
%Ibwfo = bwareaopen(Ibwf,20);
imshow(Ibwf)
drawnow;
waitforbuttonpress;
close all;

% Label components in the inverse BW image
L = generate_labels_matrix(Ibwf);

% Color components
Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
figure
imshow(Lrgb)
title('Colored watershed label matrix (Lrgb)')
drawnow;
waitforbuttonpress;
close all;

function cum_hist = get_cumulative_histogram(N, Nsum)
    cum_hist = zeros(numel(N));
    my_sum = 0;
    for i=1:numel(N)
        my_sum = my_sum + N(i);
        cum_hist(i) = my_sum/Nsum; 
    end

% -----------------------------------------------------------------------------

function threshold_intensity = get_threshold_intensity(N, frac_black, bin_locs)
    Ncurr = 0;
    Nsum = sum(N);
    threshold_count = double(Nsum) * frac_black;
    imax = numel(N);
    disp(['N(1)=' num2str(N(1)) ' N(last)=' num2str(N(imax)) ' N(pre-last)=' num2str(N(imax-1)) ' Nsum=' num2str(Nsum) ' threshold_count=' num2str(threshold_count)]);
    threshold_intensity = 0;
    num_bins = numel(N);
    i = 0;
    while i<= num_bins 
        i = i + 1;
        Ncurr = Ncurr + N(i);
        disp(['i=' num2str(i) ' frac_black=' num2str(frac_black) ' Ncurr/Nsum=' num2str(Ncurr/Nsum) ' (Ncurr + N(i))/Nsum=' num2str((Ncurr + N(i))/Nsum)]);
        if Ncurr + N(i) >= threshold_count
            deltaI = double(bin_locs(i+1) - bin_locs(i));
            deltaN = double(N(i));
            threshold_intensity = bin_locs(i) + deltaI*(double(threshold_count - Ncurr)/deltaN);
            disp(['bin_locs(i)=' num2str(bin_locs(i)) ' deltaI=' num2str(bin_locs(i+1) - bin_locs(i)) ' threshold_intensity=' num2str(threshold_intensity)]); 
            break
        end
        Ncurr = Ncurr + N(i);
    end

function output_usage_message()
    disp('Usage: MS_S2D_SegmentNeurons2D(image_file [ ,fraction_of_black ])');
    return;

function L = generate_labels_matrix(Ibwf)
    CC = bwconncomp(Ibwf); % NOTE: this fucntion does not label holes
    Nc = CC.NumObjects;
    imsize = CC.ImageSize;

    disp(['imsize=' num2str(imsize) ' num_comp=' num2str(Nc) ]);
    L = zeros(imsize);
    for k = 1:Nc
        for m=1:numel(CC.PixelIdxList{k})
            lin_ind = CC.PixelIdxList{k}(m);
            [r, c] = ind2sub(CC.ImageSize, lin_ind);
%           disp(['lin_ind=' num2str(lin_ind) ' r=' num2str(r) ' c=' num2str(c)]);
            L(r,c) = k;
        end
    end

% -----------------------------------------------------------------------------

function Ibw = my_im2bw(I, threshold_intensity)
   Ibw = logical(zeros(size(I)));
   Ibw(I > threshold_intensity) = logical(1);
   disp(['my fraction of black=' num2str((numel(Ibw)-sum(sum(Ibw)))/numel(Ibw))]);
