%Radial consensus with varying topology
%Robots have the knowledge of a fixed reference point
%edsonbffilho@gmail.com
%11-Dec-2018

clear all
close all

% Get the number of available agents from the Robotarium.  We don't need a
% specific value for this algorithm
%N = rb.get_available_agents();
%N=15;
%M = number of groups
n_abs=3;

%number of robots per group
n_robots=4;

N=n_robots*n_abs;

%Simulation or Sending [0] Simulation Mode [1] Send to Robotarium mode
send_robotarium_mode=1;

% Select the number of iterations for the experiment.  This value is
% arbitrary
iterations = 3000;

%% Prealocating
zx_hist(1:n_robots*n_abs,iterations)=0;
zy_hist(1:n_robots*n_abs,iterations)=0;
zvx_hist(1:n_robots*n_abs,iterations)=0;
zvy_hist(1:n_robots*n_abs,iterations)=0;
rrx_ca_hist(1:n_robots*n_abs,iterations)=0;
rry_ca_hist(1:n_robots*n_abs,iterations)=0;
vvx_ca_hist(1:n_robots*n_abs,iterations)=0;
vvy_ca_hist(1:n_robots*n_abs,iterations)=0;
UUx_ca_hist(1:n_robots*n_abs,iterations)=0;
UUy_ca_hist(1:n_robots*n_abs,iterations)=0;
rrx_plot(1:n_robots,1:n_abs)=0;
rry_plot(1:n_robots,1:n_abs)=0;
rrx(1:n_robots,1:n_abs)=0;
rry(1:n_robots,1:n_abs)=0;
zx(1:n_robots,1:n_abs)=0;
zy(1:n_robots,1:n_abs)=0;
R_ca(1:n_robots*n_abs,1)=0;
theta_dot_ca(1:n_robots*n_abs,1)=0;
UUx_ca(1:n_robots*n_abs,1)=0;
UUy_ca(1:n_robots*n_abs,1)=0;
Utheta(1:n_robots,1:n_abs)=0;
A(1:n_robots,1:n_robots,n_abs)=0;
D(1:n_robots,1:n_robots,n_abs)=0;
L(1:n_robots,1:n_robots,n_abs)=0;
tan_vel_ca(1:n_robots*n_abs,1)=0;
tan_vel(1:n_robots,1:n_abs)=0;
all_vec(1:n_robots,1:n_robots,n_abs)=0;
sorted_metric(1:n_robots,1:n_robots,n_abs)=0;
delta2dotX(1:n_robots,1:n_abs)=0;
delta2dotY(1:n_robots,1:n_abs)=0;
h1(1:n_robots,1:n_abs)=0;
h2(1,1:n_abs)=0;
h4(1:n_robots,1:n_abs)=0;
norm_r(1:n_robots*n_abs,1:n_robots*n_abs)=0;
n_error(1:iterations)=0;

%% Get Robotarium object used to communicate with the robots/simulator
r = Robotarium('NumberOfRobots', N, 'ShowFigure', false);
if send_robotarium_mode==1
    r = Robotarium('NumberOfRobots', N, 'ShowFigure', true);
end

%Radius of the robots
radius = 0.5*r.robot_diameter;

%% %Random color
if n_abs>10
    rand_color = rand(n_abs,3);
else
    rand_color(1:10,1:3)=[1 0 0;0 1 0;0 0 1;0 1 1;0 0 0;0.3 0.5 0.7;0.3 0.7 0;1 1 0;0.9 0.9 0.9;0.3 0.3 0.3];
end

%% Experiment constants

% Generate a cyclic graph Laplacian from our handy utilities.  For this
% algorithm, any connected graph will yield consensus
%L = cycleGL(N);

%Fixed Connections
[A_fixed_aux,L_fixed_aux] = rand_AL_sparse(n_abs,n_robots);
G_fix = graph(A_fixed_aux);
min_A=minspantree(G_fix,'Method','sparse');
A_fixed=full(adjacency(min_A));

%Degree Matrix
D_fixed=diag(sum(A_fixed,2));   
%Laplacian Matrix
L_fixed=D_fixed-A_fixed;


%% Radial Consensus Initialize 
UUx(1:n_robots,1:n_abs)=0;
UUy(1:n_robots,1:n_abs)=0;
r_memory(1:n_robots,1:n_abs)=0;
r_memory_local(1:n_abs*n_robots,1:n_abs,1:n_robots*n_abs)=0;
theta_dot(1:n_robots,1:n_abs)=0;
theta_des(1:n_robots,1:n_abs)=0;
theta_dot_des(1:n_robots,1:n_abs)=0;
theta_bar_dot(1:n_robots,1:n_abs)=0;

%Communication radius
c=0.62;

%Desired Radius of each group (Initialized with 1)
d=0.5*c-0.05;

R(1:n_robots,1:n_abs)=d;
hd(n_robots*n_abs,2,n_robots*n_abs)=0;

%Initialize delta
delta(1:2,1:N)=0;
%Initialize theta
theta_ca(1:N,1) = (2*pi).*rand(N,1);
%Transforming variables
theta(1:n_robots,1:n_abs)=0;
for i=1:n_abs
    theta(:,i)=theta_ca(n_robots*i-(n_robots-1):n_robots*i);
end

%[1] - Varying topologie  [2] - Fixed topologie
topology=1;

%% Grab tools we need to convert from single-integrator to unicycle dynamics

% % Gain for the diffeomorphism transformation between single-integrator and
% % unicycle dynamics
% transformation_gain = 0.06;
% [si_to_uni_dyn, uni_to_si_states] = create_si_to_uni_mapping('ProjectionDistance', transformation_gain);
% 
% %safety_radius = 0.20;
% 
% lambda = 0.05;
% safety_radius = 1.5*r.robot_diameter;
% %si_barrier_cert = create_si_barrier_certificate('SafetyRadius', safety_radius);
% 
% % Grab barrier certificates for unicycle dynamics
% uni_barrier_cert = create_uni_barrier_certificate('SafetyRadius', safety_radius, 'ProjectionDistance', lambda);

%% Retrieve tools for single-integrator -> unicycle mapping

% Let's retrieve some of the tools we'll need.  We would like a
% single-integrator position controller, a single-integrator barrier
% function, and a mapping from single-integrator to unicycle dynamics
position_int = create_si_position_controller('XVelocityGain', 1, 'YVelocityGain', 1);
si_barrier_certificate = create_si_barrier_certificate('SafetyRadius', 1.5*r.robot_diameter);
si_to_uni_dyn = create_si_to_uni_mapping2('LinearVelocityGain', 0.75, 'AngularVelocityLimit', pi);

i_step=0.033;

% Initialize velocity vector for agents.  Each agent expects a 2 x 1
% velocity vector containing the linear and angular velocity, respectively.
dxi = zeros(2, N);

% Initialize dxu
dxu = zeros(2, N);

%Iterate for the previously specified number of iterations
for t = 1:iterations
    
    % Retrieve the most recent poses from the Robotarium.  The time delay is
    % approximately 0.033 seconds
    x = r.get_poses();
    
    % Convert to SI states
    %xi = uni_to_si_states(x);
    xi = x(1:2,:);
    
    %Get velocities (Not posible anymore)
    %vvx_ca(:,1)=r.velocities(1,:).*cos(x(3,:));
    %vvy_ca(:,1)=r.velocities(1,:).*sin(x(3,:));
    
    %Get velocities (New way)
    vvx_ca(:,1)=dxu(1,:).*cos(x(3,:));
    vvy_ca(:,1)=dxu(1,:).*sin(x(3,:));
    
    %Transforming variables
    for i=1:n_abs
        rrx_plot(:,i)=xi(1,n_robots*i-(n_robots-1):n_robots*i);
        rry_plot(:,i)=xi(2,n_robots*i-(n_robots-1):n_robots*i);
        rrx(:,i)=x(1,n_robots*i-(n_robots-1):n_robots*i);
        rry(:,i)=x(2,n_robots*i-(n_robots-1):n_robots*i);
        theta_ca(n_robots*i-(n_robots-1):n_robots*i,1)=theta(:,i);
        R_ca(n_robots*i-(n_robots-1):n_robots*i,1)=R(:,i);
        theta_dot_ca(n_robots*i-(n_robots-1):n_robots*i,1)=theta_dot(:,i);
        UUx_ca(n_robots*i-(n_robots-1):n_robots*i,1)=UUx(:,i);
        UUy_ca(n_robots*i-(n_robots-1):n_robots*i,1)=UUy(:,i);
    end
    for i=1:n_abs
        rrx_ca(n_robots*i-(n_robots-1):n_robots*i,1)=rrx(:,i);
        rry_ca(n_robots*i-(n_robots-1):n_robots*i,1)=rry(:,i);
    end 
    
    
    %% Time varying Adjacency matrix
    %For each group separatedly
    for j=1:n_abs
        %Adjacency Matrix for each group
        A(:,:,j)=adj_mat_calculate(rrx(:,j),rry(:,j),c);
        %Degree Matrix for each group
        D(:,:,j)=diag(sum(A(:,:,j),2));
        %Laplacian Matrix for each group
        L(:,:,j)=D(:,:,j)-A(:,:,j);      
    end
    %Adjacency Matrix
    A_var=adj_mat_calculate(rrx_ca,rry_ca,c);   
    
    %% Second Order Consensus
    %Point 
    zx_ca=rrx_ca + R_ca.*cos(theta_ca);
    zy_ca=rry_ca + R_ca.*sin(theta_ca);
    
    %Point Velocity 
    zvx_ca=vvx_ca - R_ca.*theta_dot_ca.*sin(theta_ca);
    zvy_ca=vvy_ca + R_ca.*theta_dot_ca.*cos(theta_ca);
  
    %Velocity Consensus Gain
    k_gamma=8; 
    
    % This allows us to sum over agent i's neighbors
    UUx_ca(1:n_robots*n_abs,1)=0;
    UUy_ca(1:n_robots*n_abs,1)=0;
    
    if topology==1%Varying
        for i=1:n_robots*n_abs
            for j=1:n_robots*n_abs 
                UUx_ca(i)=UUx_ca(i)-A_var(i,j)*((zx_ca(i)-zx_ca(j))+k_gamma*(zvx_ca(i)-zvx_ca(j)));
                UUy_ca(i)=UUy_ca(i)-A_var(i,j)*((zy_ca(i)-zy_ca(j))+k_gamma*(zvy_ca(i)-zvy_ca(j)));
            end
        end
        %Fixed leader in point [0 0]
        UUx_ca=UUx_ca - zx_ca - k_gamma*zvx_ca;
        UUy_ca=UUy_ca - zy_ca - k_gamma*zvy_ca;
    end
    
    if topology==2%Fixed
        for i=1:n_robots*n_abs
            for j=1:n_robots*n_abs 
                UUx_ca(i)=UUx_ca(i)-A_fixed(i,j)*((zx_ca(i)-zx_ca(j))+k_gamma*(zvx_ca(i)-zvx_ca(j)));
                UUy_ca(i)=UUy_ca(i)-A_fixed(i,j)*((zy_ca(i)-zy_ca(j))+k_gamma*(zvy_ca(i)-zvy_ca(j)));
            end
        end
    end
    
    %Save hist
    zx_hist(:,t)=zx_ca;
    zy_hist(:,t)=zy_ca;
    zvx_hist(:,t)=zvx_ca;
    zvy_hist(:,t)=zvy_ca;
    rrx_ca_hist(:,t)=rrx_ca;
    rry_ca_hist(:,t)=rry_ca;
    vvx_ca_hist(:,t)=vvx_ca;
    vvy_ca_hist(:,t)=vvy_ca;
    UUx_ca_hist(:,t)=UUx_ca;
    UUy_ca_hist(:,t)=UUy_ca;
    
    %Transforming variables
    for k=1:n_abs       
        UUx(:,k)=UUx_ca(n_robots*k-(n_robots-1):n_robots*k,1);
        UUy(:,k)=UUy_ca(n_robots*k-(n_robots-1):n_robots*k,1);
        zx(:,k)=zx_ca(n_robots*k-(n_robots-1):n_robots*k,1);
        zy(:,k)=zy_ca(n_robots*k-(n_robots-1):n_robots*k,1);
    end

    %% Radius Heuristics
    %for each robot i
    for i=1:n_abs*n_robots
        %Neighbors of robot i
        vec_A=A_var(i,1:n_abs*n_robots)';%line 3
        %for each neighbor
        for j=1:n_abs*n_robots%line 3
            %If robot i communicates with robot j
            if vec_A(j)==1
            if (which_group(i,n_abs,n_robots)) > (which_group(j,n_abs,n_robots))%line 5
                %Only if the new value is bigger than the one I already have
                if (R_ca(j)+d) > R_ca(i)
                    hd(:,:,i)=0;
                    R_ca(i)=R_ca(j)+d;
                end
            end
            %Robots of the same group
            if (which_group(i,n_abs,n_robots)) == (which_group(j,n_abs,n_robots))%line 9
                %Only if the new value is bigger than the one I already have
                if (R_ca(j)+d) > R_ca(i)
                    hd(:,:,i)=0;
                    R_ca(i)=R_ca(j);
                end
            end
            %Memory save
            if (which_group(j,n_abs,n_robots)) > (which_group(i,n_abs,n_robots)) && R_ca(j)==(R_ca(i)+d)%line 13
                %Only if the new value is bigger than the one I already have
                if (R_ca(j)+d) > R_ca(i)
                    hd(j,1,i)=R_ca(j);%line 15
                    hd(j,2,i)=(which_group(j,n_abs,n_robots));%line 15
                end
            end
            
            %Analyse memory
            flag_memory=0;%If exists
            for k=1:n_abs*n_robots
                if hd(k,2,j) < which_group(i,n_abs,n_robots) && R_ca(i)<=hd(k,1,j)%line 16
                    flag_memory=1;
                end
            end
            if flag_memory==1
                    hd(:,:,i)=0;%line 17
                    R_ca(i)=R_ca(i)+d;% line 18
            end
            end
        end
    end
    
    %Transforming variables
    for k=1:n_abs       
       R(:,k)=R_ca(n_robots*k-(n_robots-1):n_robots*k,1);
    end  
    
    %% Theta controler Utheta
    %True distance from robot to reference point
    trueR=sqrt((rrx-mean(zx_ca)).^2+(rry-mean(zy_ca)).^2); 
    for k=1:n_abs
        trueR_ca(n_robots*k-(n_robots-1):n_robots*k,1)=trueR(:,k);
    end
    
    %Utheta(1:n_robots,1:n_abs)=0;
    
    %Tangencial Velocity (In relation to origin)
    for i=1:n_abs*n_robots
        tan_vel_ca(i,1)=dot((([-rry_ca(i) rrx_ca(i)])/R_ca(i)),[vvx_ca(i) vvy_ca(i)]);
    end
    
    %Transforming variables
    for k=1:n_abs       
       tan_vel(:,k)=tan_vel_ca(n_robots*k-(n_robots-1):n_robots*k,1);
    end  
    
    %Desired Theta Calculations
    % Adjusting Robots Theta
    %Calculating closer robots (by angle)
    for j=1:n_abs
        for i=1:n_robots       
            %Removing robots that I don't see
            vec_aux=wrapTo2Pi(theta(:,j)-theta(i,j)).*A(:,i,j);     
            %Sort only my group
            all_vec(:,i,j)=vec_aux;
            %sorted(:,i,j)]=sort(vec_aux);           
            %Sorted for metric
            vec_aux2=wrapTo2Pi(theta(i,j)-theta(:,j));     
            %Sort only my group
            sorted_metric(:,i,j)=sort(wrapTo2Pi(vec_aux2));  
        end
    end
    %Calculating Neighbors
    all_vec(all_vec == 0) = NaN;
    [min_r,min_r_idx]=min(all_vec);
    [max_r,max_r_idx]=max(all_vec);   
   
    %Checking robots neighbors
    for j=1:n_abs
        for i=1:n_robots
           if all(isnan(all_vec(:,i,j)))~=1
                %If robot has at least one neighbor from the same group
                theta_des(i,j)=((min_r(1,i,j) + max_r(1,i,j)-2*pi)/2);
           else
                %If robots has no neighbors from the same group
                theta_des(i,j)=0;%theta(i,j);
           end       
        end
    end   

    %Theta Controller Gains
    kd_theta=0.005;
    kp_theta=0.001;
    
    %To turn off theta controller
    %kd_theta=0;
    %kp_theta=0;
    
    %Fixed value to make robots rotate   
    %w_des=0.02*R;
    w_des=(0.03)./R;
    
    %Theta Controller
    %Utheta=kd_theta*(-theta_dot)+kp_theta*(theta_des) + (w_des);
    Utheta=kd_theta*(theta_bar_dot)+kp_theta*(theta_des)-0.1*(theta_dot-w_des);
   
%     Utheta(1:n_robots,1:n_abs)=0;
%     for j=1:n_abs
%         for i=1:n_robots
%             if (trueR(i,j)-R(i,j)).^2<0.0001
% %                 %Utheta(i,j)=kd_theta*(theta_bar_dot(i,j))+kp_theta*(theta_des(i,j))-0.1*(theta_dot(i,j)-w_des(i,j));
% %                 disp('Utheta on, robot')
% %                 disp(i)
% %                 disp('group')
% %                 disp(j)
%             end
%         end
%     end
    
    %% Little delta2dot
    for j=1:n_abs
        for i=1:n_robots
            delta2dotX(i,j)=-R(i,j)*(Utheta(i,j)*sin(theta(i,j)+(theta_dot(i,j)^2)*cos(theta(i,j))));
            delta2dotY(i,j)=+R(i,j)*(Utheta(i,j)*cos(theta(i,j)-(theta_dot(i,j)^2)*sin(theta(i,j))));   
        end
    end
    
    %Adding delta2dot to controller  
    UUx=UUx+delta2dotX;
    UUy=UUy+delta2dotY;
    
    %Transforming variables
    for k=1:n_abs       
        UUx_ca(n_robots*k-(n_robots-1):n_robots*k,1)=UUx(:,k);
        UUy_ca(n_robots*k-(n_robots-1):n_robots*k,1)=UUy(:,k);
    end
    
    %% Segregation Error (According to GRoB 2009)
    error=0;
    
    for i=1:n_abs*n_robots
      for j=1:n_abs*n_robots  
        if (which_group(i,n_abs,n_robots)) < (which_group(j,n_abs,n_robots))
                if trueR_ca(i)>=trueR_ca(j)
                    error=error+1;
                end
        else if (which_group(i,n_abs,n_robots)) > (which_group(j,n_abs,n_robots))
                if trueR_ca(i)<=trueR_ca(j)
                    error=error+1;
                end 
              end
        end
      end
    end
    %Worst possible error
    worst_error=(n_abs*n_robots)^2;

    %Normalized error
    n_error(t)=(error-0)/(worst_error-0);
    
    
    %% Transform the double-integrator controller to single-integrator dynamics
    UUx_ca=vvx_ca+i_step*UUx_ca;
    UUy_ca=vvy_ca+i_step*UUy_ca;
    
    %% Applying controllers to the robots
    dxi(1,:)=UUx_ca;
    dxi(2,:)=UUy_ca;  
    
%     %% Avoid errors
%     
%     % To avoid errors, we need to threshold dxi
%     norms = arrayfun(@(x) norm(dxi(:, x)), 1:N);
%     threshold = r.max_linear_velocity/2;
%     to_thresh = norms > threshold;
%     dxi(:, to_thresh) = threshold*dxi(:, to_thresh)./norms(to_thresh);
    
    %%
    % Normalization of controls.  This code ensures that
    %%
    % $$
    %   \|dxu\| \leq dmax
    % $$
    dxmax = 0.75*r.max_linear_velocity;
    for i = 1:N
        if norm(dxi(:,i)) > dxmax
            dxi(:,i) = dxi(:,i)/norm(dxi(:,i))*dxmax;
        end
    end
    
     %% Apply barrier certs. and map to unicycle dynamics
    
    %Ensure the robots don't collide
    dxi = si_barrier_certificate(dxi, x);
    
    % Transform the single-integrator dynamics to unicycle dynamics using a
    % diffeomorphism, which can be found in the utilities
    dxu = si_to_uni_dyn(dxi, x);        

    
    %% Transform the single-integrator to unicycle dynamics using the the
    % transformation we created earlier
    %dxi = si_to_uni_dyn(dxi, x);
    
        %% Utilize barrier certificates  
    %dxi = si_barrier_cert(dxi, xi);
    %dxu = uni_barrier_cert(dxi, xi);
    
    
    %% Draw System
    if send_robotarium_mode==1
    marker_size=70;
    %Draw all robots
    for j=1:n_abs
            for i=1:n_robots
                %h1(i,j)=circle(rrx(i,j),rry(i,j),radius,rand_color(j,1:3));
                if j==1
                   h1(i,j)=plot(rrx(i,j),rry(i,j),'p','color',rand_color(j,1:3),'MarkerSize',marker_size,'MarkerFaceColor',rand_color(j,1:3));
                end
                if j==2
                   h1(i,j)=plot(rrx(i,j),rry(i,j),'d','color',rand_color(j,1:3),'MarkerSize',marker_size,'MarkerFaceColor',rand_color(j,1:3));
                end
                if j==3
                   h1(i,j)=plot(rrx(i,j),rry(i,j),'s','color',rand_color(j,1:3),'MarkerSize',marker_size,'MarkerFaceColor',rand_color(j,1:3));
                end
                if j==4
                   h1(i,j)=plot(rrx(i,j),rry(i,j),'h','color',rand_color(j,1:3),'MarkerSize',marker_size,'MarkerFaceColor',rand_color(j,1:3));
                end
                if j==5
                   h1(i,j)=plot(rrx_plot(i,j),rry_plot(i,j),'s','color',rand_color(j,1:3));
                end
            end
    end           
            
    %if t>iterations-1000
    % Draw reference circles
    
        if topology==1%Varying
            if t<2
                %Reference point
                scatter(0,0,10,[0 0 0])
                %Reference circles
                h2(1)=circleDashed(0,0,(d),rand_color(1,1:3));
                for j=2:n_abs
                    h2(j)=circleDashed(0,0,(d+d*(j-1)),rand_color(j,1:3));
                end
            end
        end
        if topology==2%Varying
            %Reference point
            scatter(mean(zx_ca),mean(zy_ca),10,[0 0 0])
            %Reference circles
            h2(1)=circleDashed(mean(zx_ca),mean(zy_ca),(d),rand_color(1,1:3));
            for j=2:n_abs
                h2(j)=circleDashed(mean(zx_ca),mean(zy_ca),(d+d*(j-1)),rand_color(j,1:3));
            end
        end
    %end   
    
    %Delete plots to plot again
    drawnow
    delete(h1)
    if topology==2
        delete(h2)
    end
    %t
    end
    %Show picture sometimes either way (simulation only)
    if send_robotarium_mode==0
        if mod(t/100,1) == 0 || t==1
        hh=figure(1)
        
        clf(figure(hh))
        axis([-1.65 1.65 -1.05 1.05])
        axis equal
        title(t)
        hold on
            for j=1:n_abs
                for i=1:n_robots
                    circle(rrx(i,j),rry(i,j),radius,rand_color(j,1:3));
                end
            end
            circleDashed(0,0,(d),rand_color(1,1:3));
            for j=2:n_abs
                circleDashed(0,0,(d+d*(j-1)),rand_color(j,1:3));
            end
        end
    end   
    
    %% Angular Integration
    for j=1:n_abs
        for i=1:n_robots
            %theta(i,j)=theta(i,j) + theta_des(i,j)*r.time_step;
            %Angular Position
            theta(i,j)=theta(i,j) + theta_dot(i,j)*r.time_step + 0.5*Utheta(i,j)*r.time_step^2;
            
            %Angular Velocity
            theta_dot(i,j)=theta_dot(i,j)+r.time_step*Utheta(i,j);                       
           
            theta_bar_dot(i,j)=theta_des(i,j)+tan_vel(i,j)*(1/R(i,j))*i_step;
            
            %plot([rrx(i,j) zx(i,j)],[rry(i,j) zy(i,j)],'-.','color',rand_color(j,1:3));    
        end
    end   
     
    %% Send velocities to agents
    
    % Set velocities of agents 1,...,N
    r.set_velocities(1:N, dxu);
    
    % Send the previously set velocities to the agents.  This function must be called!    
    r.step();
    
     %% Verify if agents are outside bounderies
    for i=1:n_abs*n_robots
        if rrx_ca(i)>1.65 || rrx_ca(i)<-1.65 || rry_ca(i)>1.05 || rry_ca(i)<-1.05
            disp('agents outside bounderies')  
            disp(t)
            pause
        end
    end
    
    %Sometimes disp t
    if mod(t/100,1) == 0 || t==1
        t
        n_error(t)
    end
end

%Save data
save('alldata.mat')

% We can call this function to debug our experiment!  Fix all the errors
% before submitting to maximize the chance that your experiment runs
% successfully.
r.debug();