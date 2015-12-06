% Add zero padding at the adges,
% to ensure thgat that entire image != one BW component
function Ibwp = MS_S2D_AddBoundaryPadding(Ibw, value)
    size1 = size(Ibw, 1);
    size2 = size(Ibw, 2);
    Ibwp = Ibw;
    Ibwp(1:2  ,:)           = logical(value);
    Ibwp((size1-1):size1,:) = logical(value);
    Ibwp(:,            1:2) = logical(value);
    Ibwp(:,(size2-1):size2) = logical(value);

