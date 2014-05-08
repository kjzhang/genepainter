function [MASK, IND] = BTTC_filter(img, epsilon, verbose)
%{
    Our implementation of the B-tree triangular coding (BTTC) algorithm for
    image compression. The algorithm recursively decomposes the image
    domain into right-angled triangles in a binary tree structure, such
    that every point in each triangle is below the error threshold epsilon.

    The algorithm requires square images of side lengths 2^n + 1, so we pad
    the images to the next highest valid edge length.

    NOTE: BTTC_filter has a dependency on java.util.Stack()

    INPUTS:
        img, [m n c] input image, supports grayscale and RGB
        epsilon, error threshold

    OUTPUTS:
        MASK, [m n] mask of pixels in the BTTC decomposition
        IND, indices of pixels in the BTTC decomposition
%}

    if( nargin < 3 )
        verbose = 0;
    end

    [m, n, ~] = size(img);
    
    d = max([m n]);
    pad = (2 ^ nextpow2(d - 1) + 1) - d;
    d = d + pad;
    img = padarray(img, [pad pad], 'post');
    
    L_mask = zeros(d, d);
    
    % initialize stack
    stack = java.util.Stack();
    stack.push([[1, 1]; [1, d]; [d, 1]]);
    stack.push([[d, d]; [d, 1]; [1, d]]);

    while stack.empty() ~= 1
        % get the next triangle
        T = stack.pop();
        P1 = T(1, :);
        P2 = T(2, :);
        P3 = T(3, :);

        % find indices of points within the triangle
        min_x = min(T(:, 1));
        max_x = max(T(:, 1));
        min_y = min(T(:, 2));
        max_y = max(T(:, 2));
        [X, Y] = meshgrid(min_x:max_x, min_y:max_y);
        X = reshape(X, [numel(X) 1]);
        Y = reshape(Y, [numel(Y) 1]);
        [IN, ~] = inpolygon(X, Y, T(:, 1), T(:, 2));

        success = 1;
        num_points = numel(IN);
        for i = 1 : num_points
            if IN(i) == 1
                % check if triangle approximates P within epsilon
                P = [X(i) Y(i)];
                G = BTTC_G(img, P, P1, P2, P3);
                F = img(P(1), P(2), :);
                error = norm(F(:) - G(:), 2);
                if error > epsilon
                    success = 0;
                    break;
                end
            end
        end

        if success == 0
            % bad approximation, subdivide triangles
            PM = (P2 + P3) / 2;
            T1 = [PM; P1; P2];
            T2 = [PM; P3; P1];
            stack.push(T1);
            stack.push(T2);
        else
            % good approximation, add finished triangle to mask
            L_mask(P1(1), P1(2)) = 1;
            L_mask(P2(1), P2(2)) = 1;
            L_mask(P3(1), P3(2)) = 1;
        end
        
        if verbose == 1
            % print stack size
            stack.size()
        end
    end
    
    spy(L_mask);
    
    L_mask = L_mask(1 : m, 1 : n);
    
    % calculate image mask and pixel indices
    MASK = (L_mask == 1);
    IND = find(MASK == 1);
end