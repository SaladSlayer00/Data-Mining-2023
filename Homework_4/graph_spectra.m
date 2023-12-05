clear all; close all; clc;

%Read data from the dataset
E = csvread('example1.dat');
%E = csvread('example2.dat');

%Step 1

%Converting Edge list to the adjacency matrix:
%Extract the first column of the data stored in E and assigns it to the variable col1.
col1 = E(:,1);
%Extract the second column of the data stored in E and assigns it to the variable col2.
col2 = E(:,2);

%Find the maximum value between the maximum values of col1 and col2
max_ids = max(max(col1,col2));

%Create a sparse matrix As with dimensions determined by the maximum ID (max_ids) found in the edge list. It uses the values from col1 as row indices, values from col2 as column indices, and assign the value 1 at those positions.
As= sparse(col1, col2, 1, max_ids, max_ids);

%A = full(adjacency(graph(As)));

sigma = 1;

G = graph(col1, col2);
d = distances(G);
%Create Affinity matrix
A = exp(-(d.^2)./(2*sigma^2)) - eye(size(d));

figure(1)
spy(A)
figure(2)
plot(G, 'Layout','force');
figure(3)
h=plot(G,'Layout','force');

k=4;
%k=2;

%Step 2
D=diag(sum(A,2));
L=(D^(-1/2)*A*D^(-1/2));

%Step 3: Find the k largest eigenvectors of L and construct X
[eig_vec, eig_val] = eigs(L);
[eig_val_sorted, idx] = sort(diag(eig_val),'descend');
X = eig_vec(:,idx(1:k));

%Step 4: Renormalize X to have unit length
Y=X./sum(X.*X,2).^(1/2);

%Step 5: Run K-means on Y.
idx = kmeans(Y,k);

%Step 6
cluster_colors=hsv(k);
for i=1:k
    cluster_members=find(idx==i);
    highlight(h, cluster_members , 'NodeColor', cluster_colors(i,:))
end
figure(4)
[F_V,~]=eigs(L,2,'SA');
plot(sort(F_V(:,2)),'-*')