function alpha = Fun_SmoothingMatting(image, pack)

%% Parameters
omega    = 1e-1;
lambda   = 100;
epsilon  = [];
win_size = [];

%% Normalize inputs
image      = double(image);
alpha_hat  = double(pack(:,:,1));
confidence = double(pack(:,:,2));
trimap     = pack(:,:,3);

if max(image(:)) > 1
    image = image / 255;
end

if max(confidence(:)) > 1
    confidence = confidence / 255;
end

if max(alpha_hat(:)) > 1
    alpha_hat = alpha_hat / 255;
end

consts_map = ~(trimap(:,:,1) > 0 & trimap(:,:,1) < max(trimap(:)));

[h,w,c] = size(image);
img_size = w * h;

%% Generate matting Laplacian matrix
tlaplacian = start_timer('Generating matting Laplacian... ');
L = getLaplacian1(image, consts_map, epsilon, win_size);
end_timer(tlaplacian);

D = spdiags(consts_map(:), 0, img_size, img_size);
P = spdiags(confidence(:) .* ~consts_map(:), 0, img_size, img_size);

K = lambda * D + omega * P;

%% Solve for alpha
talpha = start_timer('Solving for alpha... ');
x = (L + K) \ (K * alpha_hat(:));
end_timer(talpha);

alpha = max(min(reshape(x,h,w),1),0);


function t = start_timer(msg)
	fprintf(1, msg); t = tic;

function end_timer(timer)
	t = toc(timer); fprintf(1, '%.8f sec\n', t);
