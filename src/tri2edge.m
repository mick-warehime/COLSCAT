% given a list of points and elements return the unique edges

function [edgelist,tri2edge] = tri2edge(t)

% all possible edges (including replicates)
ne = size(t,1);
alledges = sort([reshape(t',ne*3,1) reshape(t(:,[2 3 1])',ne*3,1)],2);

% use setdiff to remove duplicate edges
[edgelist,a,el2edge] = unique(alledges,'rows');

% triangulation to edge
tri2edge = reshape(el2edge,3,ne)';
