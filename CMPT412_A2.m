%
%img = im2double(imread('./imgs/OneBallLetteringVerticalLarge.jpg'));  %OK
%img = im2double(imread('./imgs/OneBallVerticalLarge.jpg'));           %OK
%img = im2double(imread('./imgs/OneBallCornerLarge.jpg'));             %OK
%img = im2double(imread('./imgs/TwoBallsVerticalLarge.jpg'));          %OK
%img = im2double(imread('./imgs/ThreeBallsNetLarge.jpg'));             %OK
%img = im2double(imread('./imgs/OneBallLarge.jpg'));                   %OK
%img = im2double(imread('./imgs/TwoBallsTouchingVerticalLarge.jpg'));  %F2A1
%img = im2double(imread('./imgs/NewBallsLarge.jpg'));                  %OK
%img = im2double(imread('./imgs/ThreeBallsCloseUpTouching.jpg'));      %F2A1
%img = im2double(imread('./imgs/TwoBallsShadowLarge.jpg'));            %FWO
%img = im2double(imread('./imgs/ThreeBallsShadowLarge.jpg'));          %F2A1

img = im2double(imread('./imgs/OneBallShadowsLarge.jpg'));            %OK

imshow(img);

[H, W, XX] = size(img);
scale = (0:2/H:1)';
scale = (scale .* 0.9) + 0.1;
scale_col = ones([1,H]);
scale_col(1:size(scale,1)) = scale;
dim_mat = repmat(scale_col,W,1)';


testim = ones(H,W) .* 2;
testim = testim .* dim_mat;

%Trying to normalize colour vectors and taking a dot product.
if ( true )
    % could probably fine-tune this per channel.
    img_cpy = img .* dim_mat;
    %img_cpy = imgaussfilt(img_cpy, 4.5);
    
    imshow(img_cpy);
    
    mask = rgb2gray(img_cpy) > 0.5;
    
    divisor = sqrt(img_cpy(:,:,1).^2 + img_cpy(:,:,2).^2 + img_cpy(:,:,3).^2);
    %should check for zero entries.
    im_unit_vecs = (img_cpy ./ divisor) .* mask;
    imshow(im_unit_vecs);
    
    v1 = [0.5565 0.7138 0.4251]; % OneBallLetteringVerticalLarge / OneBallVerticalLarge / OneBallCornerLarge

    v2 = [0.6622 0.6622 0.3506]; % TwoBallsVerticalLarge
    v3 = [0.5871 0.7205 0.3690]; % ThreeBallsNetLarge / OneBallLarge
    v4 = [0.5996 0.6957 0.3956];
    
    %v = [0.6187 0.6876 0.3799];
    v5 = [0.5478 0.7259 0.4160];
    
    v6 = [0.5320 0.7435 0.4052]; % Slight adjust for OneBallCornerLarge
    v7 = [0.6637 0.6675 0.3376];
    
    %v = [0.6494 0.6494 0.3957];
    %v = [0.6105 0.6335 0.4754];
    %v = [0.6475 0.6525 0.3937];
    %v = v ./ norm(v)
    
    vec_mat = [v1; v2; v3; v4; v5; v6; v7];
    
    gradient_map = zeros(H, W);
    
    sz = size(vec_mat, 1);
    for i=1:sz
        v = vec_mat(i,:);
        to_dot = zeros(H, W, 3);
        to_dot(:,:,1) = v(1);
        to_dot(:,:,2) = v(2);
        to_dot(:,:,3) = v(3);
    
        gradient_map = gradient_map + dot(im_unit_vecs, to_dot, 3);
    end
    
    %imshow( gradient_map ./ sz );
    
    
    filt = gradient_map > (0.992 * sz);
    
    imshow(filt);
    
    filt = bwareaopen(filt, 150); % was 75 before.
    
    imshow(filt);

    filt = imdilate(filt, strel('disk', 15));
    
    
    
    %filt = imdilate(filt, strel('disk', 5));
    
    %imshow(filt);
    
    %filt = imerode(filt, strel('disk', 10));
    
    
    
    imshow(filt);
    
    filt = imfill(filt, 'holes');
    
    imshow(filt);
    
    filt=bwmorph(filt,'shrink',Inf);
   
    
    [y,x] = find(filt > 0);
    imshow(img);
    hold on;
    plot(x, y, 'r+', 'MarkerSize', 30, 'LineWidth', 2);
    hold off;
    
    return;
end


% trying  the built in matlab find circles
if ( false )
    d = imdistline;
    delete(d)
    filt = rgb2gray(img);
    imshow(filt)
    
    [centers,radii] = imfindcircles(img, [20 60], 'ObjectPolarity', 'bright', 'Sensitivity', 0.9);
    [centers2,radii2] = imfindcircles(img, [200 250], 'ObjectPolarity', 'bright', 'Sensitivity', 0.95);
       
    imshow(img)
    h = viscircles(centers,radii);
    h = viscircles(centers2,radii2);
    
    return
end


% Trying the convolution with a circle thing
if ( false )
    w = 75; r = 31.5
    [x, y] = meshgrid(1:w, 1:w);
    f_circ = ((x - (w/2)).^2 + (y - (w/2)).^2 <= r^2);
    f_circ = double(f_circ);
    f_circ(f_circ==0) = -1;
    figure; imshow(f_circ)

    title('Plot of the circle mask')
    img_blur = img; %imgaussfilt(img, 20);
    img_blur = imgaussfilt(img_blur, 10);

%     img_blur = imgaussfilt(img, 20);
%     img_blur = rgb2hsv(img_blur);
%     img_blur(:,:,3) = histeq(img_blur(:,:,3));
%     img_blur = hsv2rgb(img_blur);

    imshow(img_blur);
    
    mask = rgb2gray(img_blur) > 0.6;
    %mask = mask & rgb2gray(img_blur) < 0.9;

    %filt(filt<0.7) = 0;
    %imshow(mask);
    
    divisor = sqrt(img_blur(:,:,1).^2 + img_blur(:,:,2).^2 + img_blur(:,:,3).^2);
    %should check for zero entries.
    filt = double((rgb2gray(img_blur) .* mask ) > 0);
    imshow(filt);
    
    filt = normalize(conv2(filt, f_circ) > 100);
    
    imshow(filt)
    v = [0.5565 0.7138 0.4251];
    %v = [0.6475 0.6525 0.3937];
    %v = v ./ norm(v)
    to_dot = img_blur;
    to_dot(:,:,1) = v(1);
    to_dot(:,:,2) = v(2);
    to_dot(:,:,3) = v(3);

    filt = dot(filt, to_dot, 3) > 0.995;

    imshow(filt);
    
    filt=bwmorph(filt,'shrink',Inf);
   
    
    [y,x] = find(filt > 0);
    imshow(img);
    hold on;
    plot(x, y, 'r+', 'MarkerSize', 30, 'LineWidth', 2);
    hold off;
    
    return;
end





%Trying to normalize colour vectors and taking a dot product.
if ( true )
    img_blur = img; %imgaussfilt(img, 20);
    img_blur = imgaussfilt(img_blur, 20);

%     img_blur = imgaussfilt(img, 20);
%     img_blur = rgb2hsv(img_blur);
%     img_blur(:,:,3) = histeq(img_blur(:,:,3));
%     img_blur = hsv2rgb(img_blur);

    imshow(img_blur);
    
    mask = rgb2gray(img_blur) > 0.4;
    %mask = mask & rgb2gray(img_blur) < 0.9;

    %filt(filt<0.7) = 0;
    %imshow(mask);
    
    divisor = sqrt(img_blur(:,:,1).^2 + img_blur(:,:,2).^2 + img_blur(:,:,3).^2);
    %should check for zero entries.
    filt = (img_blur ./ divisor) .* mask;
    imshow(filt);
    
    v = [0.5565 0.7138 0.4251];
    %v = [0.6475 0.6525 0.3937];
    %v = v ./ norm(v)
    to_dot = img_blur;
    to_dot(:,:,1) = v(1);
    to_dot(:,:,2) = v(2);
    to_dot(:,:,3) = v(3);

    filt = dot(filt, to_dot, 3) > 0.995;

    imshow(filt);
    
    filt=bwmorph(filt,'shrink',Inf);
   
    
    [y,x] = find(filt > 0);
    imshow(img);
    hold on;
    plot(x, y, 'r+', 'MarkerSize', 30, 'LineWidth', 2);
    hold off;
    
    return;
end

% comparing only the two colour channels
if (false)
    %filt = imgaussfilt(img, 7);
    %filt = 1 + -1 * normalize(abs(filt(:,:,2) - 2*filt(:,:,3)));
    %imshow(filt);
end


F = fspecial('laplacian',1)
L = [0 1 0; 1 -4 1; 0 1 0];
filt = rgb2gray(normalize(img));
%filt = imgaussfilt(filt, 1);
imshow(filt);

%filt(filt<0.6) = 0;

imshow(filt);

%filt = filt.^(1/2);

%filt = img(:,:,2);
filt = normalize((conv2(filt, L)));

imshow(filt);

%filt = filt > 0.55;




imshow(double(filt));

function res = normalize(img)
    res = img - min(img(:)); 
    MAX = max(res(:));
    if( MAX > 0 )
       res = res ./ MAX; 
    end
end