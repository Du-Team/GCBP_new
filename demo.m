clear;
close all;
src = '';
name = '';
load([src,name]);

% Parameters
no_grid=20;
precentile = 0.1;
eps = 0.05;
max_iters = 12;


% GCBP
clusters = fastgrid(data,no_grid,precentile,eps,max_iters);

if (find(clusters == -1)) % False if not found
    clusters_count = length(unique(clusters)) - 1;
else
    clusters_count = length(unique(clusters));
end

fprintf("The number of clusters: %d \n", clusters_count);
