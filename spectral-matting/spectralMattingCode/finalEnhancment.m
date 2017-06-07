function alpha = finalEnhancment(alpha,eigs_num,eig_vectors,eig_values);
% do some iterations of the Newton's method, to get a sparser alpha matte

su = eig_vectors;
sd = eig_values;
scld=1;

if (~exist('exp_a','var'))
exp_a=0.9;
end

thr_e=0.1^10;
w1=0.3; w0=0.3; wa=1;
orig_exp_a = exp_a;

x=alpha(:);
for i=1:30
  e1=w1.^exp_a*max(abs(x-1),thr_e).^(exp_a-2);
  e0=w0.^exp_a*max(abs(x),thr_e).^(exp_a-2);
  ssu=repmat(e1+e0,1,eigs_num).*su;
  y=(su'*ssu+scld*sd)\(su'*e1);
  x=su*y;
  if (i==20)
    exp_a=exp_a-0.2;
  end  
end
[n,m]=size(alpha);
alpha=reshape(x,n,m);

exp_a = orig_exp_a;

