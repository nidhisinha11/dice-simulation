clear all
close all

%% Create the cube
r = 1;
total_mass = 1;
[kmax,lmax,X_rel,jj,kk,M,I] = cube(r,total_mass);

%% Create the ground
%ground will be plane through origin with unit normal [ag,bg,cg]
ag = 0; bg = 0; cg = 1;
nabcg = sqrt(ag^2+bg^2+cg^2);
ag = ag/nabcg; bg = bg/nabcg; cg = cg/nabcg;  %normalize

S_ground = 10^4;    %stiffness of ground (Kg/s^2)
D_ground = 0.1;     %damping of ground (Kg/s^2)
mu_ground = .1;     %friction coefficient of ground

%set cutoff velocity for ground friction to avoid
%division by zero later if tangential velocity is zero
umin = (10^-6)*r*ones(kmax,1);

%% Initial setup
g = 9.8;       
x_cm = [0,0,3];                 %initial position
u_cm = [5,0,5];                 %initial velocity
L = [0.5; 0.5; 0.5];   
%L = rand(3,1)*0.75 + (10^-6);   %random initial angular momentum
omega = I\L;                    %initial angular velocity

%create matrices for node positions and velocities
X = repmat(x_cm,kmax,1) + X_rel;
U = repmat(u_cm,kmax,1);

%% Initial drawing
figure(1)
link = 1:lmax;
h = plot3( [X(jj(link),1), X(kk(link),1)]', [X(jj(link),2), X(kk(link),2)]',...
           [X(jj(link),3), X(kk(link),3)]', 'b' );
%plot settings
grid on
xlabel('x'); ylabel('y')
axis equal; axis([-5,30,-10,10,0,5])
set(gcf, 'Position', get(0, 'Screensize'));
drawnow

%% Run the simulation
%time and recording
nskip = 40;            %frame iteration for animation
tmax = 10;              %duration of simulation (s)
clockmax = 20000;      %number of time steps
dt = tmax/clockmax;     %step size (s)

%for future plotting
X_save = zeros(clockmax,3);
t_save = zeros(clockmax,1);

for clock=1:clockmax
    t = clock*dt;
    
    %% Add gravity to each node
    F = zeros(kmax,3);
    F(:,3) = - (total_mass/8)*g;
    
    %% Add force of ground to each node
    %for each node, evaluate:
    dg = max(0, -(ag*X(:,1)+bg*X(:,2)+cg*X(:,3)) );   %distance into the ground
    un = ag*U(:,1) + bg*U(:,2) + cg*U(:,3);           %velocity normal to ground
    fn = max(0,S_ground*dg - D_ground*un);            %magnitude of normal force
    
    Utan = U - un*[ag,bg,cg];                         %tangential velocity vector
    Utan_norm = vecnorm(Utan,2,2) + umin;             %norm of tangential velocity
    Utan_dir = Utan./[Utan_norm,Utan_norm,Utan_norm]; %normalize
    
    %add normal force and friction force
    F = F + fn*[ag,bg,cg] - mu_ground*[fn,fn,fn].*Utan_dir;
    
    %% Update governing equations
    %update center of mass
    u_cm = u_cm + dt*sum(F,1)./total_mass;
    x_cm = x_cm + dt*u_cm;
    %update angular momentum (changes due to force of ground)
    torque = zeros(3,1);
    for k=1:kmax
        torque = torque + cross(X_rel(k,:),F(k,:))';
    end
    L = L + dt*torque;
    
    %% Calculate for animation  
    %update X_rel due to spinning
    omega = I\L;
    P = (omega./norm(omega))*(omega./norm(omega))';
    omega_cross = [0,-omega(3),omega(2); omega(3),0,-omega(1);
        -omega(2),omega(1),0];
    R = P + cos(norm(omega)*dt)*(eye(3)-P) + ...
        sin(norm(omega)*dt)*(omega_cross/norm(omega));
    for k=1:kmax
        X_rel(k,:) = (R*X_rel(k,:)')';
    end

    %update nodes velocities and positions
    product = zeros(kmax,3);
    for k=1:kmax
    	product(k,:) = cross(omega',X_rel(k,:));
    end
    U = repmat(u_cm,kmax,1) + product;  %needed for force of the ground
    X = repmat(x_cm,kmax,1) + X_rel;    %needed for animation
    
    %% Update animation
    if(mod(clock,nskip)==0)
        c = 0;
        for i=link
            c = c+1;
            h(c).XData = [X(jj(i),1),X(kk(i),1)];
            h(c).YData = [X(jj(i),2),X(kk(i),2)];
            h(c).ZData = [X(jj(i),3),X(kk(i),3)];
        end
        drawnow
    end
    
    %% Store results for future plotting
    X_save(clock,:) = x_cm;
    t_save(clock) = t;
    kinetic_save(clock) = .5*total_mass*norm(u_cm)^2 + .5*norm(L)^2 / I(1,1);
    potential_save(clock) = total_mass*g*x_cm(3);
    energy_save(clock) = kinetic_save(clock) + potential_save(clock);
  
end

%% Graph trajectory
figure(2)
plot(t_save',X_save')
figure(3)
plot(t_save',energy_save')