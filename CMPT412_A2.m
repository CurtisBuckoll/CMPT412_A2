% -------------------------------------------------------------------------
% Read in a source image
img = im2double(imread('./imgs/OneBallLetteringVerticalLarge.jpg'));  %OK
%img = im2double(imread('./imgs/OneBallVerticalLarge.jpg'));           %OK
%img = im2double(imread('./imgs/OneBallCornerLarge.jpg'));             %OK
%img = im2double(imread('./imgs/TwoBallsVerticalLarge.jpg'));          %OK
%img = im2double(imread('./imgs/ThreeBallsNetLarge.jpg'));             %OK
%img = im2double(imread('./imgs/OneBallLarge.jpg'));                   %OK
%img = im2double(imread('./imgs/TwoBallsTouchingVerticalLarge.jpg'));  %OK
%img = im2double(imread('./imgs/NewBallsLarge.jpg'));                  %OK
%img = im2double(imread('./imgs/ThreeBallsCloseUpTouching.jpg'));      %OK
%img = im2double(imread('./imgs/TwoBallsShadowLarge.jpg'));            %OK
%img = im2double(imread('./imgs/ThreeBallsShadowLarge.jpg'));          %OK
%img = im2double(imread('./imgs/OneBallShadowsLarge.jpg'));            %OK
%img = im2double(imread('./imgs/extra/4.jpg'));
% -------------------------------------------------------------------------

imshow(img);

% -------------------------------------------------------------------------
% The main idea here is to interpret each RGB pixel as a 3 component vector 
% with some geometric orientation. We can then sample a number of different 
% images to determine pixel values of the tennis balls in different
% lighting environments (in the sun, shade, etc). These vectors are all
% normalized to have unit length, and then dot products can be taken for
% each pixel to compare them to the similarity between source to the target
% pixels. The pixels we are interested are those with the score closest to
% 1.
% -------------------------------------------------------------------------

img_cpy = img;

% Is the top half of the screen overall brighter than the bottom? If so,
% dim the top half. This as an attempt to remove false positives.
[H, W, XX] = size(img);
scale = (0:1.5/H:1)';
scale = (scale .* 0.9) + 0.1;
scale_col = ones([1,H]);
scale_col(1:size(scale,1)) = scale;
dim_mat = repmat(scale_col,W,1)';
sum_top = sum(sum(sum(img(1:floor(H/2),:,:))));
sum_bot = sum(sum(sum(img(ceil(H/2):H,:,:))));
vert_intensity_ratio = sum_top / sum_bot;

% Set a hard ratio threshold of 1.5 If above, we dim.
if ( vert_intensity_ratio > 1.5 )
    img_cpy = img_cpy .* dim_mat;
end

imshow(img_cpy);

% Remove any pixels with intensity below 0.5 It seems we can be pretty
% generous here.
mask = rgb2gray(img_cpy) > 0.5;

% Create the 'gradient' map as a result of taking per pixel dot products,
% and threshold to take only those pixels that score above 0.992.
[gradient_map, sz] = tb_gradient_map(img_cpy);
gradient_map = gradient_map .* mask;
filt = gradient_map > (0.992 * sz);

imshow(filt);

% Try to clean up any leftover noise.
filt = imerode(filt, strel('disk', 1));
filt = bwareaopen(filt, 150);

% Expand and fill the regions to create closed areas.
filt = imdilate(filt, strel('diamond', 15));
filt = imfill(filt, 'holes');

imshow(filt);

% Extract single pixels, ideally representing the TB centers.
[Y,X] = find(bwmorph(filt,'shrink',Inf) > 0);

% This loop sequence is to try to break one closed area that actually
% contains two TB's into two seperate areas. We do this by looking at
% the width and height of each region and examining the ratios. This
% will probably fail if the balls are not horizontally or vertically
% aligned or if more than two TB's are in a line.
for i=1:size(Y,1)
    seg_w = 1;
    seg_h = 1;
    start_x = X(i);
    start_y = Y(i);
    while true
        if filt(start_y, start_x - 1) == 0
            break;
        end
        seg_w = seg_w + 1;
        start_x = start_x - 1;
    end
    start_x = X(i);
    while true
        if filt(start_y, start_x + 1) == 0
            break;
        end
        seg_w = seg_w + 1;
        start_x = start_x + 1;
    end
    start_x = X(i);
    while true
        if filt(start_y - 1, start_x) == 0
            break;
        end
        seg_h = seg_h + 1;
        start_y = start_y - 1;
    end
    start_y = Y(i);
    while true
        if filt(start_y + 1, start_x) == 0
            break;
        end
        seg_h = seg_h + 1;
        start_y = start_y + 1;
    end

    h_to_w = seg_h / seg_w;
    w_to_h = seg_w / seg_h;

    % If we find that width is at least twice the height or vice-versa,
    % then split the region into two.
    if ( h_to_w > 2.02 )
        filt(Y(i)-1:Y(i)+1, X(i)-seg_w:X(i)+seg_w) = 0;
    elseif ( w_to_h > 2.02 )
        filt(Y(i)-seg_h:Y(i)+seg_h, X(i)-1:X(i)+1) = 0;
    end
end

imshow(filt);

% Extract single pixels once more in case we have seperated the regions in
% the loop above. Plot the final points on the source image.
filt_shrink = bwmorph(filt,'shrink',Inf);
[y,x] = find(filt_shrink > 0);
imshow(img);
hold on;
plot(x, y, 'r+', 'MarkerSize', 30, 'LineWidth', 1.5);
hold off;

% Output the radii.
radii = find_radius(filt, x, y)





