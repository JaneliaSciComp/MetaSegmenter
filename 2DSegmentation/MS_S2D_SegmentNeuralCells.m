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

frac_black = 0.34;
if nargin > 1
    frac_black = double(varargin{2});
end 

% Original image    
I = imread(image_file); 
I = mat2gray(I(:,:,1));
disp(['min(I)=' num2str(min(min(I))) ' max(I)=' num2str(max(max(I))) ...
      ' mean2(I)=' num2str(mean2(I)) ' std2(I)=' num2str(std2(I))]);
In = round(double(I)*255.0/max(max(double(I))));
disp(['min(In)=' num2str(min(min(In))) ' max(In)=' num2str(max(max(In))) ...
      ' mean2(In)=' num2str(mean2(In)) ' std2(In)=' num2str(std2(In))]);
imshow(I);
title('Original image (I)');
waitforbuttonpress;
close all;

% Determine intensity threshold for conversion to black-white image
% NOTE: imhist assumes that I is a grayscale image, with range of values [0, 1]
% NOTE2: edges are computed incorrectly by Matlab's imhist function,
%        so a custom function has been implemented for compuring the edges

num_edges = 1001;
N     = imhist(I,num_edges);      % counts
edges = compute_edges(num_edges); % using a custom function for edges
Nsum  = sum(N);
disp(['numel(N)=' num2str(numel(N)) ' Nsum=' num2str(Nsum) ...
      ' numel(I)=' num2str(numel(I)) ' numel(edges)=' num2str(numel(edges))]);
cum_hist = get_cumulative_histogram(N, Nsum);

% Black-white image
threshold = get_threshold_intensity(cum_hist, frac_black, edges);
disp(['final threshold_intensity=' num2str(threshold)]);
disp(['size(I)=' num2str(size(I)) ' class(I)=' class(I(1,1))]);
Ibw = im2bw(I, threshold);  
disp(['size(Ibw)=' num2str(size(Ibw)) ' class(Ibw)=' class(Ibw(1,1))]);
Ibw = bwareaopen(Ibw,20); 
disp(['my frac_black=' num2str(sum(sum(Ibw == 0))/numel(Ibw))]);
figure
imshow(Ibw), title('Black-white image (Ibw)')
drawnow;
waitforbuttonpress;
close all;

% Fill holes and open gaps
Ibw = adjust_edges(Ibw);
Ibwf = imfill(Ibw, 'holes');
Ibwf = bwareaopen(Ibwf,20);
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

% -----------------------------------------------------------------------------

function cum_hist = get_cumulative_histogram(N, Nsum)
    cum_hist = zeros(numel(N));
    my_sum = 0;
    for i=2:numel(N)
        my_sum = my_sum + N(i-1);
        cum_hist(i) = my_sum/Nsum; 
    end
    % Adjusting the last bin
    my_sum = my_sum + N(i);
    cum_hist(numel(N)) = my_sum/Nsum;

% -----------------------------------------------------------------------------

function threshold_intensity = get_threshold_intensity(cum_hist, frac_black, edges)
    threshold_intensity = 0;
    num_bins = numel(cum_hist);
    i = 0;
    while i<= num_bins-1
        i = i + 1;
        if cum_hist(i) >= frac_black      
            deltaI = double(edges(i) - edges(i-1));
            deltaH = cum_hist(i) - cum_hist(i-1);
            threshold_intensity = edges(i-1) + deltaI*(frac_black - cum_hist(i-1))/deltaH;           
            disp(['edges(i)=' num2str(edges(i)) ' cum_hist(i)=' num2str(cum_hist(i))...
                  ' threshold_intensity=' num2str(threshold_intensity)]);
            break
        end
    end

% -----------------------------------------------------------------------------

function output_usage_message()
    disp('Usage: MS_S2D_SegmentNeurons2D(image_file [ ,fraction_of_black ])');
    return;

% -----------------------------------------------------------------------------

function L = generate_labels_matrix(Ibwf)
    CC = bwconncomp(Ibwf); % NOTE: this function does not label holes
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
   Ibw(I > round(threshold_intensity)) = logical(1);
   disp(['threshold_intensity=' num2str(threshold_intensity) ...
         ' my fraction of black=' num2str((double(numel(Ibw)-sum(sum(Ibw))))/double(numel(Ibw)))]);

% -----------------------------------------------------------------------------

function edges = compute_edges(num_edges)
    edges = zeros(1, num_edges);
    step = 1/double(num_edges-1);
    edges(1) = 0.;
    edges(2) = step/2.;
    edges(num_edges  ) = 1.;
    edges(num_edges-1) = 1. - step/2.;
    for i=3:(num_edges-2)
        edges(i) = edges(i-1) + step;
    end

% -----------------------------------------------------------------------------

% Add zero padding at the adges, 
% to ensure thgat that entire image != one BW component
function Ibw = adjust_edges(Ibw)
    size1 = size(Ibw, 1);
    size2 = size(Ibw, 2);

    i = 1;
    while sum(Ibw(i,:)) >= size2*0.9
        Ibw(i,:) = logical(0);
        i = i+1;
    end

    i = 1;
    while sum(Ibw(:,i)) >= size1*0.9
        Ibw(:,i) = logical(0);
        i = i+1;
    end

    i = 0;
    while sum(Ibw(size1 - i,:)) >= size2*0.9
        Ibw(size1 - i,:) = logical(0);
        i = i+1;
    end

    i = 0;
    while sum(Ibw(:,size2 - i)) >= size1*0.9
        Ibw(:,size2 - i) = logical(0);
        i = i+1;
    end

