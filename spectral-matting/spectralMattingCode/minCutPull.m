function [assign,e] = minCutPull(Wpairs, W0, W1)
 
Nnodes = length(W0);

tweights = [[1:Nnodes]',W0,W1];
[vi,vj] = meshgrid(1:Nnodes,1:Nnodes);
idxs = find((vi ~= vj).*(Wpairs ~= 0));

edges = [vi(idxs),vj(idxs),Wpairs(idxs),Wpairs(idxs)];

[e,labeling]=maxflow(edges,tweights);

% disp(sprintf('Error of the assignment = %f', e));
assign = labeling'-1;
size(assign)
