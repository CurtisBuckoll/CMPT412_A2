%
%img = im2double(imread('./imgs/OneBallLetteringVerticalLarge.jpg'));  %OK
%img = im2double(imread('./imgs/OneBallVerticalLarge.jpg'));           %OK
%img = im2double(imread('./imgs/OneBallCornerLarge.jpg'));             %OK
%img = im2double(imread('./imgs/TwoBallsVerticalLarge.jpg'));          %OK
%img = im2double(imread('./imgs/ThreeBallsNetLarge.jpg'));             %OK
%img = im2double(imread('./imgs/OneBallLarge.jpg'));                   %OK
img = im2double(imread('./imgs/TwoBallsTouchingVerticalLarge.jpg'));  %F2A1
%img = im2double(imread('./imgs/NewBallsLarge.jpg'));                  %OK
%img = im2double(imread('./imgs/ThreeBallsCloseUpTouching.jpg'));      %F2A1
%img = im2double(imread('./imgs/TwoBallsShadowLarge.jpg'));            %FWO
%img = im2double(imread('./imgs/ThreeBallsShadowLarge.jpg'));          %F2A1
%img = im2double(imread('./imgs/OneBallShadowsLarge.jpg'));            %OK

%img = im2double(imread('./imgs/extra/4.jpg'));

imshow(img);

[H, W, XX] = size(img);
scale = (0:1.5/H:1)'; % can make this scale adaptively by taking into consideration the top/bottom intensity ratio.
scale = (scale .* 0.9) + 0.1;
scale_col = ones([1,H]);
scale_col(1:size(scale,1)) = scale;
dim_mat = repmat(scale_col,W,1)';


%Trying to normalize colour vectors and taking a dot product.
if ( true )
    
    img_cpy = img;
    
    % is the top half of the screen darker?
    teess = img(1:floor(H/2),:,:);
    sum_top = sum(sum(sum(img(1:floor(H/2),:,:))))
    sum_bot = sum(sum(sum(img(ceil(H/2):H,:,:))))
    vert_intensity_ratio = sum_top / sum_bot
    if ( vert_intensity_ratio > 1.5 )
        % could probably fine-tune this per channel.
        img_cpy = img_cpy .* dim_mat;
    end
    
    rough_mask = rgb2gray(img_cpy) > 0.4;
    img_blurred = imgaussfilt(img_cpy, 10);
    [rough_gradient, sz] = tb_gradient_map(img_blurred);
    rough_gradient = rough_gradient .* rough_mask;
    rough_area = rough_gradient > (0.990 * sz);
    
    imshow(rough_area);
    rough_area = imerode(rough_area, strel('disk', 1));
    imshow(rough_area);
    rough_area = bwareaopen(rough_area, 100);
    imshow(rough_area);
    rough_area = imdilate(rough_area, strel('disk', 50));
    imshow(rough_area);
    rough_area = imfill(rough_area, 'holes');
    imshow(rough_area);
    
    %img_cpy = imgaussfilt(img_cpy, 4.5);
    %img_cpy(:,:,1) = imadjust(img_cpy(:,:,1));
    %img_cpy(:,:,2) = imadjust(img_cpy(:,:,2));
    %img_cpy(:,:,3) = imadjust(img_cpy(:,:,3));
    
    imshow(img_cpy);
    
    mask = rgb2gray(img_cpy) > 0.5;
    
%     divisor = sqrt(img_cpy(:,:,1).^2 + img_cpy(:,:,2).^2 + img_cpy(:,:,3).^2);
%     %should check for zero entries.
%     im_unit_vecs = (img_cpy ./ divisor) .* mask;
%     imshow(im_unit_vecs);
%     
%     v1 = [0.5565 0.7138 0.4251]; % OneBallLetteringVerticalLarge / OneBallVerticalLarge / OneBallCornerLarge
% 
%     v2 = [0.6622 0.6622 0.3506]; % TwoBallsVerticalLarge
%     v3 = [0.5871 0.7205 0.3690]; % ThreeBallsNetLarge / OneBallLarge
%     v4 = [0.5996 0.6957 0.3956];
%     
%     %v = [0.6187 0.6876 0.3799];
%     v5 = [0.5478 0.7259 0.4160];
%     
%     v6 = [0.5320 0.7435 0.4052]; % Slight adjust for OneBallCornerLarge
%     v7 = [0.6637 0.6675 0.3376];
%     
%     %v = [0.6494 0.6494 0.3957];
%     %v = [0.6105 0.6335 0.4754];
%     %v = [0.6475 0.6525 0.3937];
%     %v = v ./ norm(v)
%     
%     vec_mat = [v1; v2; v3; v4; v5; v6; v7];
%     
%     gradient_map = zeros(H, W);
%     
%     sz = size(vec_mat, 1);
%     for i=1:sz
%         v = vec_mat(i,:);
%         to_dot = zeros(H, W, 3);
%         to_dot(:,:,1) = v(1);
%         to_dot(:,:,2) = v(2);
%         to_dot(:,:,3) = v(3);
%     
%         gradient_map = gradient_map + dot(im_unit_vecs, to_dot, 3);
%     end
    
    %imshow( gradient_map ./ sz );
    
    [gradient_map, sz] = tb_gradient_map(img_cpy);
    gradient_map = gradient_map .* mask;
    filt = gradient_map > (0.992 * sz);
    
    imshow(filt);
    
    filt = imerode(filt, strel('disk', 1)); % this was added lol - Skews the result for NewBallsLarge/OneBallShadowsLarge, but still works.
    
    imshow(filt);
    
    filt = bwareaopen(filt, 150); % was 75 before.
    
    imshow(filt);

    filt = imdilate(filt, strel('diamond', 15)); % maybe try a square?
    
    
    
    %filt = imdilate(filt, strel('disk', 5));
    
    %imshow(filt);
    
    %filt = imerode(filt, strel('disk', 10));
    
    
    
    imshow(filt);
    
    filt = imfill(filt, 'holes');
    
    imshow(filt);
    
    % --- New test code
    
    [Y,X] = find(bwmorph(filt,'shrink',Inf) > 0);
    
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
        
        h_to_w = seg_h / seg_w
        w_to_h = seg_w / seg_h
        
        if ( h_to_w > 2.02 )
            filt(Y(i)-1:Y(i)+1, X(i)-seg_w:X(i)+seg_w) = 0;
            rough_area(Y(i)-1:Y(i)+1, X(i)-seg_w:X(i)+seg_w) = 0;
        elseif ( w_to_h > 2.02 )
            filt(Y(i)-seg_h:Y(i)+seg_h, X(i)-1:X(i)+1) = 0;
            rough_area(Y(i)-seg_h*2:Y(i)+seg_h*2, X(i)-1:X(i)+1) = 0;
        end
        
    end
    
    imshow(filt);
    imshow(rough_area);
    
    % --- end
    
    filt_shrink = bwmorph(filt,'shrink',Inf);
  
    
    [y,x] = find(filt_shrink > 0);
    imshow(img);
    hold on;
    plot(x, y, 'r+', 'MarkerSize', 30, 'LineWidth', 1.5);
    hold off;
    
    radii = find_radius(filt, x, y)
    
    % probably not bother with this.. improves half the time but makes
    % worse the other half. just roll with the initial findings.
    [X,Y] = find_centers(rough_area, x, y);
    imshow(img);
    hold on;
    plot(X, Y, 'r+', 'MarkerSize', 30, 'LineWidth', 1.5);
    hold off;
    
    radii = find_radius(filt, X', Y')
    
    %     % watershed thing ***********************************************
%     %use watershed analysis
%     filt2 = -bwdist(~filt,'euclidean');
%     %imshow(filt2);
%     filt2(~filt) = -inf;  %set background to be infinitely far away
%     filt2 = watershed(filt2);
%     %imshow(filt2)
%     
%     %imshow(img)
%     %get the region properties
%     props = regionprops(filt2,'centroid','equivdiameter','perimeter');
%     %plot the results
%     figure;
%     imshow(img);
%     hold on;
%     %loop over all detected objects
%     for ai = 1:length(props)
%         %check if it is a circle
%         %if abs(props(ai).Perimeter-pi*props(ai).EquivDiameter)/props(ai).EquivDiameter > 13
%             %plot the detected centroid
%             plot(props(ai).Centroid(1), props(ai).Centroid(2), 'r+', 'MarkerSize', 30, 'LineWidth', 2);
%         %end
%     end
%     hold off;
%     
%     return
%     % //watershed thing *********************************************
    
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