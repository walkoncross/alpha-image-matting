
load my_images
load scribbles
%%
% if this variable is set to 0, only the final alpha matting is saved.
% Otherwise, partial results (such as the matting components) are saved.
save_partial_results = 0;

% try to group the matting components for forground/background matting
do_grouping = 1;

% for all examples here we use 50 eigen vectors...
eigs_num=50;
%%

if (1)
    nclust =10;
    [alpha_comps,alpha] = SpectralMatting(my_images.woman, [], 'woman_', eigs_num, nclust, ...
        [], save_partial_results);
end

%%

if (1)
  nclust =10;
  [alpha_comps,alpha] = SpectralMatting(my_images.face, [], 'face_', eigs_num, nclust, ...
      [], save_partial_results);
end

%%

if (1)
  nclust = 8;
  [alpha_comps,alpha] = SpectralMatting(my_images.kim, scribbles.kim, 'kim_', eigs_num, nclust, ...
      [], save_partial_results);
end

%%

if (1)
  nclust = 40;
  apply_final_enhancement = 1;
  [alpha_comps,alpha] = SpectralMatting(my_images.wind, scribbles.wind, 'wind_', eigs_num, nclust, ...
      apply_final_enhancement, save_partial_results);
end

%%

if (1)
  nclust = 20;
  [alpha_comps,alpha] = SpectralMatting(my_images.kid, scribbles.kid, 'kid_', eigs_num, nclust, ...
      [], save_partial_results);
end

