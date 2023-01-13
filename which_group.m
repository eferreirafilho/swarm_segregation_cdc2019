%which_group(robot,n_abs,n_robots)
%Function that returns the group that a robot belongs
%Works with robot indexes concatenated
%[robot1, robot2, robot3, robot4]
%-----group1-----|----group2------

function [group_index]=which_group(robot,n_abs,n_robots)

    for i=1:n_robots*n_abs
        aux_ca(i)=i;   
    end


  for k=1:n_abs
        aux(n_robots*k-(n_robots-1):n_robots*k)=aux_ca(:,k);
  end
  
  group_index=aux(robot);
end