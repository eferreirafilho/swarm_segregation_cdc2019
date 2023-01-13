%Returns a (Sparse?) random connected Adjacency matrix with connected groups
function [A_full,L_full,A] = rand_AL_sparse(n_abs,n_robots)
L_full=0;
A_full(1:n_robots*n_abs,1:n_robots*n_abs)=0;
flag_connected=0;
flag_connected_full=0;

%Full graph is connected?
while(flag_connected_full~=1)
    i1=randi(n_robots*n_abs);
    i2=randi(n_robots*n_abs);
    if i1~=i2
        A_full(i1,i2)=1;
    end
    
    %Mirroring the matrix
    A_full=triu(A_full)+triu(A_full)';
    %Degree Matrix
    D_full=diag(sum(A_full,2));   
    %Laplacian Matrix
    L_full=D_full-A_full;
    
    %Check full graph connectivity
    flag_connected_full=check_connectivity(L_full);
end

%Groups are connected?
% flag_connected_groups=0;
% while(flag_connected_groups~=(n_abs))
   
    %For each group separatedly
%     for k=1:n_abs
%         %Adjacency Matrix for each group
%         A(:,:,k)=A_full(n_robots*k-(n_robots-1):n_robots*k,n_robots*k-(n_robots-1):n_robots*k);
%         %Degree Matrix for each group
%         D(:,:,k)=diag(sum(A(:,:,k),2));
%         %Laplacian Matrix for each group
%         L(:,:,k)=D(:,:,k)-A(:,:,k);
%         
%         %Is group k connected?
%         if check_connectivity(L(:,:,k))~=1;
%             %If not connected
%             i1=randi([n_robots*k-(n_robots-1) n_robots*k],1,1);
%             i2=randi([n_robots*k-(n_robots-1) n_robots*k],1,1);
%             %Add random connection
%             if i1~=i2
%                 A_full(i1,i2)=1;
%                 A_full(i2,i1)=1;
%             end
%             A_full;
%             %Mirroring the matrix
%             %A_full=triu(A_full)+triu(A_full)';
%              
%         end
%     end
%     flag_connected_groups=0;
%     for k=1:n_abs
%         if check_connectivity(L(:,:,k))~=0
%             flag_connected_groups=flag_connected_groups+1;
%         end
%     end
% end






    

