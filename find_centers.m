function [x, y] = find_centers(img, X, Y)

    for i=1:size(Y,1)
        xlo = X(i);
        xhi = X(i);
        ylo = Y(i);
        yhi = Y(i);
        
        start_x = X(i);
        start_y = Y(i);
        while true
            if img(start_y, start_x - 1) == 0
                xlo = start_x;
                break;
            end
            start_x = start_x - 1;
        end
        start_x = X(i);
        while true
            if img(start_y, start_x + 1) == 0
                xhi = start_x;
                break;
            end
            start_x = start_x + 1;
        end
        start_x = X(i);
        while true
            if img(start_y - 1, start_x) == 0
                ylo = start_y;
                break;
            end
            start_y = start_y - 1;
        end
        start_y = Y(i);
        while true
            if img(start_y + 1, start_x) == 0
                yhi = start_y;
                break;
            end
            start_y = start_y + 1;
        end
        
        x(i) = floor((xlo + xhi) / 2);
        y(i) = floor((ylo + yhi) / 2);
        
    end
end