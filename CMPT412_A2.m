%
%img = im2double(imread('./imgs/OneBallLetteringVerticalLarge.jpg'));
%img = im2double(imread('./imgs/OneBallVerticalLarge.jpg'));
img = im2double(imread('./imgs/OneBallCornerLarge.jpg'));
%img = im2double(imread('./imgs/TwoBallsVerticalLarge.jpg'));
%img = im2double(imread('./imgs/ThreeBallsNetLarge.jpg'));
%img = im2double(imread('./imgs/OneBallLarge.jpg'));
%img = im2double(imread('./imgs/TwoBallsTouchingVerticalLarge.jpg'));
%img = im2double(imread('./imgs/NewBallsLarge.jpg'));
%img = im2double(imread('./imgs/ThreeBallsCloseUpTouching.jpg'));
%img = im2double(imread('./imgs/TwoBallsShadowLarge.jpg'));

imshow(img);

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