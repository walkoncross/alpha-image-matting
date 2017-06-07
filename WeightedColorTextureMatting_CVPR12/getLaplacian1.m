function [A]=getLaplacian1(I,consts,epsilon,win_size)
  
  if (~exist('epsilon','var'))
    epsilon=0.0000001;
  end  
  if (isempty(epsilon))
    epsilon=0.0000001;
  end
  if (~exist('win_size','var'))
    win_size=1;
  end     
  if (isempty(win_size))
    win_size=1;
  end     

  neb_size=(win_size*2+1)^2; % neighbour size
  [h,w,c]=size(I);
  img_size=w*h;
  consts=imerode(consts,ones(win_size*2+1)); % erode scribbles
  
  indsM=reshape((1:img_size),h,w); % create indices for all pixels in image, column first.

  tlen=sum(sum(1-consts(win_size+1:end-win_size,win_size+1:end-win_size)))*(neb_size^2); % calculate how many unknown pixels should be estimated.

  row_inds=zeros(tlen ,1);
  col_inds=zeros(tlen,1);
  vals=zeros(tlen,1);
  len=0;
  for j=1+win_size:w-win_size
    for i=win_size+1:h-win_size
      if (consts(i,j))
        continue % ignore any pixels under scribble regions.
      end  
      win_inds=indsM(i-win_size:i+win_size,j-win_size:j+win_size); % get neighbour pixel indices of pixel(i,j).
      win_inds=win_inds(:);
      winI=I(i-win_size:i+win_size,j-win_size:j+win_size,:); % get neighbour pixel matrix from original image
      winI=reshape(winI,neb_size,c); % reshape neighbour pixel matrix into a neb_size*c matrix, each column represents a color channel.
      win_mu=mean(winI,1)'; % get mean value of each channel.
      win_var=inv(winI'*winI/neb_size-win_mu*win_mu' +epsilon/neb_size*eye(c)); % see equation 12.
      
      winI=winI-repmat(win_mu',neb_size,1); % see equation 12.
      tvals=(1+winI*win_var*winI')/neb_size; % see equation 12.
      
      row_inds(1+len:neb_size^2+len)=reshape(repmat(win_inds,1,neb_size),...
                                             neb_size^2,1);
      col_inds(1+len:neb_size^2+len)=reshape(repmat(win_inds',neb_size,1),...
                                             neb_size^2,1);
      vals(1+len:neb_size^2+len)=tvals(:);
      len=len+neb_size^2;
    end
  end  
    
  vals=vals(1:len); % cut unused space.
  row_inds=row_inds(1:len); % cut unused space.
  col_inds=col_inds(1:len); % cut unused space.
  A=sparse(row_inds,col_inds,vals,img_size,img_size); % Create the laplacian matrix, holding the sum of second part of equation 12. Note that the row,colomn pairs may have repeated position, so all elements of vals that have duplicate row and column'indices are added together, see the help.
  
  sumA=sum(A,2);
  A=spdiags(sumA(:),0,img_size,img_size)-A; % Sum of ach row of Laplacian matrix is zero.
  
return


