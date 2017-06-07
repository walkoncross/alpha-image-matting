function falpha = decideFB(alpha);

falpha = alpha;

ones_mat = ones(size(alpha));
border_size = sum(ones_mat(:))-sum(sum(ones_mat(2:end-1,2:end-1)));
border_alpha = sum(alpha(:))-sum(sum(alpha(2:end-1,2:end-1)));

flipped = (border_alpha/border_size) > 0.5;


if (flipped)
  falpha = 1-falpha;
end

