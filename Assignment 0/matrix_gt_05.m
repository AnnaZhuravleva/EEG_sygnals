function A_bin = matrix_gt_05(A)
    A_bin = false(5, 6);
    for i = 1:5
        for j = 1:6
            if A(i, j) > 0.5;
                A_bin(i, j) = 1;
            end
        end
    end
end