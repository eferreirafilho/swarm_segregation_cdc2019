function [adj_mat]=adj_mat_calculate(rrx,rry,radius)

n_abs=size(rrx,2);
n_robots=size(rrx,1);

%Transforming dist_variables
for k=1:n_abs
rrx_ca1(n_robots*k-(n_robots-1):n_robots*k)=rrx(:,k);
rry_ca1(n_robots*k-(n_robots-1):n_robots*k)=rry(:,k);
end


%Distance between agent i with agent j
for i=1:n_abs*n_robots
for j=1:n_abs*n_robots         
    norm_r(i,j)=sqrt((rrx_ca1(i) - rrx_ca1(j))*(rrx_ca1(i) - rrx_ca1(j)) + (rry_ca1(i) - rry_ca1(j))*(rry_ca1(i) - rry_ca1(j)));
end
end

%Obtaining adjacency matrix of robots (connectivity graph)
adj_mat(n_robots,n_robots)=0;
for i=1:n_robots
     for j=1:n_robots
        if(norm_r(i,j))<radius && i~=j
            adj_mat(i,j)=1;
        end
     end
end
