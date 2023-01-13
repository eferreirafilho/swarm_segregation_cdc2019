function [flag_connection]=check_connectivity(L)

%Check if graph is connected (Second Smallest Eigenvalue > 0)
flag_connection=1;

eig_adj_robots_aux=eig(L);
eig_adj_robots(1:length(L))=sort(eig_adj_robots_aux);


%Number of connected components
[eigen_laplacian, alg_multiplicity] = eigval(L);

for i=1:length(eigen_laplacian)
    if eigen_laplacian(i)==0
        zero_eigen_idx=i;
    end
end
        
n_connected=alg_multiplicity(zero_eigen_idx);
%n_connected -> number of connected components on each graph

if size(eig_adj_robots,2)>1
    if eig_adj_robots(2)<=0 || n_connected~=1
        str_warning=['Warning: not connected!'];
        %disp(str_warning)
        %pause;  
        flag_connection=0;    
    end
end