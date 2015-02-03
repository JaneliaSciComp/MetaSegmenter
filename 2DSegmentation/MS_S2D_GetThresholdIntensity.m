%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

%
% Purpose: detect the location of "dark" structures, 
% such as mitochondria, T-,E-bars, and glial cells 
% Determine intensity threshold for conversion to black-white image
% NOTE: imhist assumes that I is a grayscale image, with range of values [0, 1]
% NOTE2: edges are computed incorrectly by Matlab's imhist function,
%        so a custom function has been implemented for compuring the edges

function threshold = MS_S2D_GetThresholdIntensity(Igr, varargin)
frac_black = 0.5;
num_edges  = 1001;
if nargin == 0 || nargin > 3
    output_usage_message();
    return
elseif nargin == 2
    frac_black = varargin{1};
elseif nargin == 3
    frac_black = varargin{1};
    num_edges  = varargin{2};
end

N     = imhist(Igr,num_edges);      % counts
edges = compute_edges(num_edges); % using a custom function for edges
Nsum  = sum(N);
disp(['numel(N)=' num2str(numel(N)) ' Nsum=' num2str(Nsum) ...
      ' numel(Igr)=' num2str(numel(Igr)) ' numel(edges)=' num2str(numel(edges))]);
cum_hist = get_cumulative_histogram(N, Nsum);

% Black-white image
threshold = get_threshold_intensity(cum_hist, frac_black, edges);

% -----------------------------------------------------------------------------

function output_usage_message()
    disp('Usage: MS_S2D_GetThresholdIntensity(grayScaleImage [,desired_fraction_of_black [,num_edges]])');
    return

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

