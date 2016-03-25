% Add zero padding at the adges,
% to ensure thgat that entire image != one BW component
function Ip = MS_S2D_AddBoundaryPadding(I, value)
    size1 = size(I, 1);
    size2 = size(I, 2);
    Ip = I;
    Ip(1:2  ,:)           = logical(value);
    Ip((size1-1):size1,:) = logical(value);
    Ip(:,            1:2) = logical(value);
    Ip(:,(size2-1):size2) = logical(value);

