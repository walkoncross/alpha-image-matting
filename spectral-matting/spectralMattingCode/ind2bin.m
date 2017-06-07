function bits = ind2bin(idx, nclust);
 
  

tmp_idx = idx;
bits = zeros(1,nclust);
for i = 1:nclust,
    bits(i) = mod(tmp_idx,2);
    tmp_idx = floor(tmp_idx / 2);
end
