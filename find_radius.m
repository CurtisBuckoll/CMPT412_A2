function r = find_radius(img, X, Y)

    for i=1:size(Y,1)
        seg_w = 1;
        seg_h = 1;
        start_x = X(i);
        start_y = Y(i);
        while true
            if img(start_y, start_x - 1) == 0
                break;
            end
            seg_w = seg_w + 1;
            start_x = start_x - 1;
        end
        start_x = X(i);
        while true
            if img(start_y, start_x + 1) == 0
                break;
            end
            seg_w = seg_w + 1;
            start_x = start_x + 1;
        end
        start_x = X(i);
        while true
            if img(start_y - 1, start_x) == 0
                break;
            end
            seg_h = seg_h + 1;
            start_y = start_y - 1;
        end
        start_y = Y(i);
        while true
            if img(start_y + 1, start_x) == 0
                break;
            end
            seg_h = seg_h + 1;
            start_y = start_y + 1;
        end
        
        if( seg_w > seg_h )
            r(i) = seg_w / 2;
        else
            r(i) = seg_h / 2;
        end   
    end
end
