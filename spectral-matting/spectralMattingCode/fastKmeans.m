function [idx, C, sumd, D]=fastKmeans(X,K);
% X: points in the N-by-P data matrix
% idx - an N-by-1 vector containing the cluster indices of each point
%
% c - the K cluster centroid locations in a K-by-P matrix.
% sumd - the within-cluster sums of point-to-centroid distances in a 1-by-K vector.
% distances from each point to every centroid in a N-by-K.

ttime = cputime;
startK = 5;
startK = min(K,startK);
maxIters = 100;  % Defualt of matlab is 100    
FinalMaxItr = 10;  % Maximal Number of Iteration for the last stage.

[idx, C, sumd, D]=kmeans(X,startK,'EmptyAction','singleton','Start','cluster', ...
    'Maxiter', maxIters, 'Replicates', 7);
disp('done initial kmeans');

valid_vec = zeros(1,startK);
scr_vec = zeros(1,startK)-1;

for nclust = startK+1:K
    % create a new cluster by splitting each cluster to two...
     max_scr=-1;
     clear min_C;
     for cl = 1:nclust-1
       cl_mask = idx == cl;
       cl_idxs = find(cl_mask);
       clX = X(cl_idxs,:);
       if (size(clX,1)> 2*size(clX,2))
         if (valid_vec(cl) == 0)
           [tmp_idx, tmp_C, tmp_sumd, tmp_D]=kmeans(clX,2,'EmptyAction','singleton','Start','cluster', 'Maxiter', maxIters);
           % chk how much the partition helps ...
           scr=sum(min(D(cl_idxs,:),[],2))-sum(min(tmp_D,[],2));
           scr_vec(cl) = scr;
         else % we already saved it...
           scr = scr_vec(cl);
         end         
       else
         scr=-2;
         scr_vec(cl) = scr;
       end  
       if (scr > max_scr)
         if (valid_vec(cl)==1) % not for the scr. Just for the idxs.
           [tmp_idx, tmp_C, tmp_sumd, tmp_D]=kmeans(clX,2,'EmptyAction','singleton','Start','cluster', 'Maxiter', maxIters);
         end
         
         max_scr = scr;
         bestC = [C;tmp_C(2,:)];  
         bestC(cl,:) = tmp_C(1,:);
         best_cl = cl;         

         best_idx = idx;
         best_idx(cl_idxs) = (tmp_idx == 1)*best_cl + (tmp_idx == 2)*nclust;
       end
       valid_vec(cl) = 1;
     end        
     C = bestC;
     idx = best_idx;
     
     valid_vec = [valid_vec,0];   % the two new clusers are new, so their
     valid_vec(best_cl) = 0;      % score have not been computed yet.
     scr_vec = [scr_vec,-1];
     scr_vec(best_cl) = -1;

     if (nclust < 13)
       [idx, C, sumd, D]=kmeans(X,nclust,'EmptyAction','singleton','Start',C, ...
                                'Maxiter', maxIters);       
       valid_vec = zeros(1,nclust);
     end
     
     cur_err = sum(min(D,[],2));
end
%[idx, C, sumd, D]=kmeans(X,K,'EmptyAction','singleton','Start',C, ...
%                         'Maxiter',FinalMaxItr);
