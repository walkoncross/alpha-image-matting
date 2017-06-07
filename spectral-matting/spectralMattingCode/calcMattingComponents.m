function [alpha_comps,eig_vectors,eig_values] = calcMattingComponents(I,L,bname,save_partial_results,...
                                                nclust,eigs_num,maxItr,exp_a)
if (~exist('maxItr','var')|| isempty(maxItr))
  maxItr=20;
end
if (~exist('exp_a','var'))
    exp_a=0.9;
end                                         
                                         
[n,m,c]=size(I);
N=n*m;
OPT.disp = 0;
[u00,d00]=eigs(L,eigs_num+1,'sm', OPT);
disp('done computing eigen vectors');
if (d00(1,1) > d00(end,end))
    d00=d00(end:-1:1,end:-1:1); u00=u00(:,end:-1:1);
end
u=u00; d=d00;

% we use only 20 first eigen vectors for the kmeans initialization.
s_eigs_num = 20;
eigs_weights =  diag(1./diag(d00(2:s_eigs_num+1,2:s_eigs_num+1).^0.5));
eigss=u00(:,2:s_eigs_num+1)*eigs_weights;

[tkmind,tctrs,tsumd,tdis]=fastKmeans(eigss,nclust);
if (save_partial_results)
    imwrite(reshape((tkmind-1)/(max(tkmind(:))-1), [n,m]),[bname,'unsp_kmeans.tif']);
end

for i=1:nclust
  x(:,i)=(tkmind==i);
end  
thr_e=0.1^10;
w1=0.3; w0=0.3; wa=1;
e1=w1.^exp_a*max(abs(x-1),thr_e).^(exp_a-2);
e0=w0.^exp_a*max(abs(x),thr_e).^(exp_a-2);

scld=1;
eig_vectors=u(:,1:eigs_num); eig_values=d(1:eigs_num,1:eigs_num);

omaxItr=maxItr;
maxItr=ceil(omaxItr/4);

% in each iteration: use only "active" components (not zero).
for nzitr=1:3
  % comute matting component with sparsity prior.
  iterateProjComps
  
  % remove matting components which are close to zero.
  nzii=find(max(abs(x))>0.1);
  nclust=length(nzii);
  x=x(:,nzii);
  e1=w1.^exp_a*max(abs(x-1),thr_e).^(exp_a-2);
  e0=w0.^exp_a*max(abs(x),thr_e).^(exp_a-2);
end
iterateProjComps

alpha_comps=x;
for l = 1:nclust
  compname = sprintf('%scomponents%.2d.tif', bname, l);
  if (save_partial_results)
      imwrite(reshape(alpha_comps(:,l),n,m),compname);
  end
end

