function [u] = demosaicing_Dummy(g,lambda)
% input: g: single bayer-filtered image
%        lambda: parameter
% output: u: demosaiced image

% dummy output
u = repmat( 0*g, [1 1 3] );

end

