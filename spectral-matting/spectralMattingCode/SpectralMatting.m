function [alpha_comps,alpha] = SpectralMatting(I, scribbles, bname, eigs_num, nclust,...
                                              apply_final_enhancement, save_partial_results)

if (~exist('eigs_num','var') | isempty(eigs_num))
  eigs_num=50
end
if (~exist('nclust','var')||isempty(nclust))
  nclust = 10;
end
if (~exist('apply_final_enhancement','var') | isempty(apply_final_enhancement))
  apply_final_enhancement = 0;
end
if (~exist('save_partial_results','var') | isempty(save_partial_results))
  save_partial_results = 0;
end

% compute laplacian matrix.
disp('computing Laplacian matrix');
[n,m,c]=size(I);
sig=0.1^5;
L=getLaplacian1(I,zeros(n,m),sig);

% compute matting components.
disp('computing matting components');
[alpha_comps,eig_vectors,eig_values] = calcMattingComponents(I,L,bname,save_partial_results,...
                                                            nclust,eigs_num);

% group together matting components to get forground-background matte.           
if (nargout > 1)
    disp('grouping matting components');
    if (isempty(scribbles)) % unsupervised grouping
        alpha=unsp_GroupMattingComponents(L,alpha_comps,I,bname, save_partial_results);
    else % supervised grouping
        alpha=sup_GroupMattingComponents(L,alpha_comps,I,scribbles, bname, save_partial_results);
    end
    if (apply_final_enhancement)
        alpha = finalEnhancment(alpha,eigs_num,eig_vectors,eig_values);
    end
    imwrite(alpha,[bname,'alpha.tif']);
end
