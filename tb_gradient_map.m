function [gmap, num_vecs] = tb_gradient_map(img)
    [H, W, XX] = size(img);
    % probably should check for zero entries.
    divisor = sqrt(img(:,:,1).^2 + img(:,:,2).^2 + img(:,:,3).^2);
    
    im_unit_vecs = (img ./ divisor); % .* mask;
    imshow(im_unit_vecs);

    v1 = [0.5565 0.7138 0.4251]; % OneBallLetteringVerticalLarge / OneBallVerticalLarge / OneBallCornerLarge
    v2 = [0.6622 0.6622 0.3506]; % TwoBallsVerticalLarge
    v3 = [0.5871 0.7205 0.3690]; % ThreeBallsNetLarge / OneBallLarge
    v4 = [0.5996 0.6957 0.3956];
    v5 = [0.5478 0.7259 0.4160];
    v6 = [0.5320 0.7435 0.4052]; % Slight adjust for OneBallCornerLarge
    v7 = [0.6637 0.6675 0.3376];

    vec_mat = [v1; v2; v3; v4; v5; v6; v7];
    gmap    = zeros(H, W);
    sz      = size(vec_mat, 1);
    
    for i=1:sz
        v = vec_mat(i,:);
        to_dot = zeros(H, W, 3);
        to_dot(:,:,1) = v(1);
        to_dot(:,:,2) = v(2);
        to_dot(:,:,3) = v(3);
    
        gmap = gmap + dot(im_unit_vecs, to_dot, 3);
    end
    
    num_vecs = sz;
end