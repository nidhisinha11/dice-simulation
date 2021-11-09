clear all
close all

%% Create the cube
r = 1;
total_mass = 1;
[kmax,lmax,X_rel,jj,kk,M,I] = cube(r,total_mass);

%% Create the ground
%ground will be plane through origin with unit normal [a,b,c]
ag = 0; bg = 0; cg = 1;
nabcg = sqrt(ag^2+bg^2+cg^2);
ag = ag/nabcg; bg = bg/nabcg; cg = cg/nabcg;  %normalize

    S_ground = 10^4;     %stiffness of ground (Kg/s^2)
    D_ground = 0.2;      %damping of ground (Kg/s^2)
    mu_ground = 0.075;   %friction coefficient of ground

%set cutoff velocity for ground friction to avoid
%division by zero later if tangential velocity is zero
umin = (10^-6)*r*ones(kmax,1);

%% Initial setup
g = 9.8;
x_cm = [0,0,3];     %initial position
u_cm = [5,0,5];     %initial velocity

%create matrices for node positions and velocities
X = repmat(x_cm,kmax,1) + X_rel;
U = repmat(u_cm,kmax,1);

%% Initial drawing
figure(1)
link = 1:lmax;
h = plot3( [X(jj(link),1), X(kk(link),1)]', [X(jj(link),2), X(kk(link),2)]',...
           [X(jj(link),3), X(kk(link),3)]', 'b' );
grid on
xlabel('x')
ylabel('y')
axis equal
axis([-5,30,-5,5,0,5])
drawnow

%% Animation
%time and recording
    nskip = 10;              %frame iteration for animation
    tmax = 10;               %duration of simulation (s)
    clockmax = 10000;        %number of time steps
dt = tmax/clockmax;          %(s)

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
    
    %% Update velicities and positions
    %for center of mass
    u_cm = u_cm + dt*sum(F,1)./total_mass;
    x_cm = x_cm + dt*u_cm;   
    
    %for each node
    X = repmat(x_cm,kmax,1) + X_rel;
    U = repmat(u_cm,kmax,1);

    %% Update animation
    if(mod(clock,nskip)==0)
        c = 0;
        for i=link %I don't understand what this means but it works
            c = c+1;
            h(c).XData = [X(jj(i),1),X(kk(i),1)];
            h(c).YData = [X(jj(i),2),X(kk(i),2)];
            h(c).ZData = [X(jj(i),3),X(kk(i),3)];
        end
        drawnow
    end
    
    %for future plotting
    X_save(clock,:) = x_cm;
    t_save(clock) = t;
    kinetic_save(clock) = .5*total_mass*norm(u_cm)^2;
    potential_save(clock) = total_mass*g*x_cm(3);
    energy_save(clock) = kinetic_save(clock) + potential_save(clock);
end

figure(2)
plot(t_save',X_save')
figure(3)
plot(t_save',energy_save')