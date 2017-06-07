function alpha=sup_GroupMattingComponents(L,alpha_comps,I,scribbles,bname,save_partial_results)
% groups components together.                          

% This function uses the mincut code from 
%   http://www.adastral.ucl.ac.uk/~vladkolm/software.html
% also described in the paper: 
% An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Computer Vision.
% Yuri Boykov and Vladimir Kolmogorov.
% In IEEE Transactions on Pattern Analysis and Machine Intelligence, September 2004

SizeThresh = 0.35;
nclust = size(alpha_comps,2);
[n,m,c]=size(I);
N=n*m;

mask1=(sum(abs(I-scribbles),3)>0.01).*(min(scribbles,[],3)>0.99);
mask0=(sum(abs(I-scribbles),3)>0.01).*(min(scribbles,[],3)<0.01);

W0 = zeros(nclust,1);
W1 = zeros(nclust,1);    
mm0=zeros(n*m,1); mm1=zeros(n*m,1);
for i=1:nclust  
  s0=sum(alpha_comps(:,i).*(alpha_comps(:,i)>0.2).*mask0(:));
  s1=sum(alpha_comps(:,i).*(alpha_comps(:,i)>0.2).*mask1(:));
  W0(i) = (s0>s1)*(s0>1)*10000;
  W1(i) = (s1>s0)*(s1>1)*10000;
  mm0=mm0+alpha_comps(:,i).*double((s0>s1)*(s0>1));
  mm1=mm1+alpha_comps(:,i).*double((s1>s0)*(s1>1));
end

free_inds=find(max(W1,W0)==0);
free_clust_len=length(free_inds);

% the number of matting comp' is small 
% enough for exhostive (2^K) enumeration.
if (free_clust_len<=20) 
  perrM = zeros(nclust);
  for i=1:nclust
    clust_size(i)=sum(sum(alpha_comps(:,i)));
    for j=1:i
      terrp=alpha_comps(:,i)'*L*alpha_comps(:,j);
      perrM(i,j)=terrp; 
      perrM(j,i)=terrp; 
    end
  end  
  bits=W1>0;
  scr=ones(1,2^free_clust_len);
  for subg=1:2^free_clust_len
    tbits = ind2bin(subg,free_clust_len );
    bits(free_inds)=tbits';
    scr(subg)=bits'*perrM*bits;   
  end  

  [mv,mi]=min(scr);
  subg=mi;
    
  tbits = ind2bin(subg, free_clust_len);
  bits(free_inds)=tbits';

else

    % the number of matting components is big: use mincut to find
    % the best grouping
    perrM = zeros(nclust);
    for i=1:nclust
        for j=1:i-1
            terrp=max(0,-alpha_comps(:,i)'*L*alpha_comps(:,j));
            perrM(i,j)=terrp;
            perrM(j,i)=terrp;
        end
    end
    disp('using min-cut !!!');
    [bits,e] = minCutPull(perrM, W0, W1);
end


alpha=reshape(alpha_comps*bits,n,m);

if (save_partial_results)
    imwrite(alpha,[bname,'alpha_before_final_enhancment.tif']);
end

