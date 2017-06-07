function alpha=unsp_GroupMattingComponents(L,alpha_comps,I,bname,save_partial_results)
% groups components together.                          

SizeThresh = 0.35;
nclust = size(alpha_comps,2);
[n,m,c]=size(I);
N=n*m;

% compute matting cost between pairs
perrM = zeros(nclust);
for i=1:nclust
    clust_size(i)=sum(sum(alpha_comps(:,i)));
    for j=1:i
      terrp=alpha_comps(:,i)'*L*alpha_comps(:,j);
      perrM(i,j)=terrp; 
      perrM(j,i)=terrp; 
    end
  end  
free_clust_len = nclust;

% find grouping for which the sum of costs is small, and 
% whose wize is within SizeThresh to (1-SizeThresh) if the image.
scr=ones(1,2^(free_clust_len-1));
sizes = ones(1,2^(free_clust_len-1));
for subg=1:2^(free_clust_len-1)
  bits = ind2bin(subg,free_clust_len)';
  scr(subg)=bits'*perrM*bits; 
  sum_alpha = sum(clust_size*bits);
  sizes(subg) = sum_alpha;
  if ((sum_alpha < SizeThresh*N) ||  (sum_alpha > (1-SizeThresh)*N))
      scr(subg) = 10^7;
  end
end  

nbest = 10;
[mv,mi]=sort(scr);

best_x = zeros(N,nbest);
for l = 1:nbest,
    subg=mi(l);
    bits = ind2bin(subg,free_clust_len)'; 
    alpha=reshape(alpha_comps*bits,n,m);
    alpha = decideFB(alpha);
    best_x(:,l)=alpha(:);
end
if (save_partial_results)
    imwrite(flat3DArray(reshape(best_x,n,m,nbest),2),[bname,'alpha_all_cands.tif']);
end

% now work on the best one...
best_idx = 1;
subg=mi(best_idx);
bits = ind2bin(subg,free_clust_len)'; 
alpha=reshape(alpha_comps*bits,n,m);
alpha = decideFB(alpha);

if (save_partial_results)
    imwrite(alpha,[bname,'best_alpha.tif']);
end

