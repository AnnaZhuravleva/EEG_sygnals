%% 1

% 1.1
A = rand(5, 6);
for i = 1:5
    for j = 1:6
        if A(i, j) > 0.5;
            disp([num2str(A(i, j)) ' is greater than 0.5']);
        else
            disp([num2str(A(i, j))  ' is less than 0.5']);
        end
    end
end

% 1.2
A_bin = false(5, 6);
for i = 1:5
    for j = 1:6
        if A(i, j) > 0.5;
            A_bin(i, j) = 1;
        end
    end
end

A
A_bin

clear A A_bin
% 1.3
A = randn(5, 6);
A_bin = matrix_gt_05(A)

A
A_bin

%% 2
% 2.1 
Moscow = imread('Moscow.bmp');

figure
imagesc(Moscow);
hold on
edgecolors = 'rgb';
for cur_color = 1:3
    cur_color_matrix = Moscow(:, :, cur_color);
    [max_val, max_idx] = max(cur_color_matrix);
    [max_max_val, max_max_idx] = max(max_val);
    idx_x = max_idx(max_max_idx);
    POINTS(cur_color).x = idx_x;
    POINTS(cur_color).y = max_max_idx;
    POINTS(cur_color).z = Moscow(idx_x, max_max_idx, cur_color);
    plot(idx_x, max_max_idx, 'p','linewidth',10,'markerEdgeColor',edgecolors(cur_color));
end
POINTS(1)
POINTS(2)
POINTS(3)
