% iterateProjComp - Finding sparse linear combination of the eigenvectors
% using Newton's method.

if (~exist('exp_a','var'))
exp_a=0.9;
end
thr_e=0.1^10;
w1=0.3; w0=0.3; wa=1;

for itr=1:maxItr
  tA=zeros((nclust-1)*eigs_num);
  tb=zeros((nclust-1)*eigs_num,1);
  eig_vectors=u(:,1:eigs_num); 
  for k=1:nclust-1
    weighted_eigs=repmat(e1(:,k)+e0(:,k),1,eigs_num).*eig_vectors;
    tA((k-1)*eigs_num+1:k*eigs_num,(k-1)*eigs_num+1:k*eigs_num)= ...
          eig_vectors'*weighted_eigs+scld*eig_values;
    tb((k-1)*eigs_num+1:k*eigs_num)=eig_vectors'*e1(:,k);
  end 
  k=nclust;
  weighted_eigs=repmat(e1(:,k)+e0(:,k),1,eigs_num).*eig_vectors;
    
  ttA=eig_vectors'*weighted_eigs+scld*eig_values;
  ttb=eig_vectors'*e0(:,k)+scld*sum(eig_vectors'*L,2);
  tA=tA+repmat(ttA,[nclust-1,nclust-1]);
  tb=tb+repmat(ttb,[nclust-1,1]);  

  y=reshape(tA\tb,eigs_num,nclust-1);
  
  x=u(:,1:eigs_num)*y;
  x(:,nclust)=1-sum(x(:,1:nclust-1),2);
  if (itr<maxItr)
    e1=w1.^exp_a*max(abs(x-1),thr_e).^(exp_a-2);
    e0=w0.^exp_a*max(abs(x),thr_e).^(exp_a-2);
  end
end  
