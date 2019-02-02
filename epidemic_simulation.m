%Epidemic Simulation -- Project 2
%Introduction to Computer Smulations
%Nadine Soliman

% total number of individuals
N=100;
% infectivity
a=0.9; 
% recovery rate
b=[0.1/10 0.1/10 0.1/10];
%death rate
d = [0.1/2 0.1/2 0.1/2];

replication_rate = 10e10;
length =9600;
mutation_rate = 10e-4;
chance_of_resitance = 1e-12;
viral_load = zeros(N, 3);
critical = [7 7 7];
critical1 = [7 0 0];
critical2= [7 7 0];
alpha = zeros(N, 1);

dt=0.1; 
tmax=100;
clockmax=ceil(tmax/dt);

color = zeros(N, 3);

%black
dead = [0 0 0];

%turquoise
immune = [0 1 1];
% color when individual develops immunity to 1 or 2 strains, but not all
% pink
immune2 = [1 0 1];

%white
alive = [1 1 1];

%red
all_strains = [1 0 0];

%yellow
one_strain  = [1 1 0];

%orange
two_strains = [255 153 51]/255;

%for loop for setting the initial colors to white
for i= 1:N
    color(i,:) = alive;
end

population = zeros(N,3);
immunity = zeros(N, 3);
initial_infected = 2;

% initial conditions
for i= 1:initial_infected
    population(i ,1) = 1;
    color(i ,:) = one_strain;
end

%counts
dead_count =0;
healthy_count =N-initial_infected;
infected_count=initial_infected;
immune_count =0;
partial_immune =0;
cost = 0;
number_of_treatments_per_time = 3;
%average cost per 0.1 of a day for HCV treatment
cost_per_treatment =number_of_treatments_per_time* 50;

%x, y positions of each individual 
x=zeros(1,N);
y=zeros(1,N);

%velocities at which each individual will move
u = zeros(1,N);
v = zeros(1,N);

%the space at which individuals begin moving
space = 100;

% fasted allocated speed for each individual
maxspeed = 5;

%limits where individuals are allowed to be
min= -100;
maximum = 100;

% initializiations for positions and velocities
for i=1:N
    x(i) = - space+ rand*2*space;
    y(i) = - space+ rand*2*space;
    u(i) = - maxspeed+ rand*2*maxspeed;
    v(i) = - maxspeed+ rand*2*maxspeed;
end

%intitializing the vector field for movement of individuals 
rand_value = randn;
[xx,yy] = meshgrid(min:10:maximum,min:10:maximum);
f = sin(2*xx/100 + 2*yy/100)+ rand_value;
x_fac = sin(f);
y_fac = cos(f);
uu = x_fac;
vv = y_fac;

%spread of disease 
for clock = 1: clockmax
    for i = 1:N  
        %first check that they are alive, if not skip them
        if color(i,:) == dead
            continue
        else 
            %updating positions
            ff = sin(2*x(i)/100 + 2*y(i)/100)+ rand_value;
            x_fac = sin(ff);
            y_fac = cos(ff);
            x(i)= x(i) + u(i)*x_fac;
            y(i)= y(i) + v(i)*y_fac;

            %making sure individuals dont exit space // they bounce off edges
            if x(i)> maximum || x(i)<min
                 if x(i)> maximum 
                    x(i) = x(i) -5;
                elseif x(i)<min
                    x(i) = x(i) + 5;
                 end
                v(i)= v(i)*-1;
                u(i)= u(i)*-1;
            end
            if y(i)> maximum || y(i)<min
                if y(i)> maximum 
                    y(i) = y(i) -5;
                elseif y(i)<min
                    y(i) = y(i) + 5;
                end
                    u(i) = u(i)*-1;
                    v(i)= v(i)*-1;
            end
            %mutation
            %skip healthy and fully immune individuals 
            if population(i,:) == [0 0 0]
                continue
            elseif immunity(i,:) == [1 1 1]
                continue
            else
                if rand<0.1*replication_rate *length* mutation_rate* chance_of_resitance
                    population(i,2) = 1;
                    if population(i,:) == [1 1 0]
                        color(i,:) = two_strains;
                    else if population(i,:) == [0 1 0]
                        color(i,:) = one_strain;  
                    end
                end
                if rand<0.1*replication_rate *length* mutation_rate* chance_of_resitance
                    population(i,3) = 1;
                end
            end
            alpha = max(viral_load,[],N);
            % recovery
                %cure for strain a
            if population(i,:) == [1 0 0]
                viral_load(i,1) = viral_load(i,1)+1;
                cost= cost+(cost_per_treatment);
                if rand< b(1)
                    if immunity(i,:) == [0 0 0]
                        partial_immune = partial_immune+1;
                        infected_count = infected_count-1;
                    end
                    immunity(i,:) = [1 0 0];
                    population(i,:) = [0 0 0];
                    %individual has recovered
                    color(i,:) = immune2;
                    viral_load(i,:) = [0 0 0];
                end
%                 %cure for strain b
            elseif population(i,:) == [0 1 0]
                 viral_load(i,2) = viral_load(i,2)+1;
                 cost= cost+(1*cost_per_treatment);
                 if rand< b(2)
                    %individual has recovered
                     if immunity(i,:) == [0 0 0]
                        partial_immune = partial_immune+1;
                        infected_count = infected_count-1;
                    end
                    immunity(i,:) = [0 1 0]
                    population(i,:) = [0 0 0];
                    color(i,:) = immune2;
                    viral_load(i,:) = [0 0 0];
                 end
             %cure for strain c
            elseif population(i,:) == [0 0 1] 
                 viral_load(i,3) = viral_load(i,3)+1;
                 cost= cost+(1*cost_per_treatment);                 
                 if rand< b(3)
                    %individual has recovered
                    if immunity(i,:) == [0 0 0]
                        partial_immune = partial_immune+1;
                        infected_count = infected_count-1;
                    end
                    immunity(i,:) = population(i,:);
                    population(i,:) = [0 0 0];
                    color(i,:) = immune2;
                    viral_load(i,:) = [0 0 0];
                  end
             %cure for strain ab
            elseif logical(population(i,:)) == logical([1 1 0]) 
                 viral_load(i,1:2) = viral_load(i,1:2)+1;
                 cost= cost+(1*cost_per_treatment);
                 if rand< (b(1))
                    %individual has recovered
                    if immunity(i,:) == [0 0 0]
                        partial_immune = partial_immune+1;
                        infected_count = infected_count-1;
                    end
                    immunity(i,:) = population(i,:);
                    population(i,:) = [0 0 0];
                    color(i,:) = immune2;
                    viral_load(i,:) = [0 0 0];
                 end
             %cure for strain ac
            elseif logical(population(i,:)) ==logical([1 0 1]) 
                 viral_load(i,1) = viral_load(i,1)+1;
                 viral_load(i,3) = viral_load(i,3)+1;
                 cost= cost+(1*cost_per_treatment);
                 if rand< (b(1))
                    %individual has recovered
                    if immunity(i,:) == [0 0 0]
                        partial_immune = partial_immune+1;
                        infected_count = infected_count-1;
                    end
                    immunity(i,:) = [1 0 1];
                    population(i,:) = [0 0 0];
                    color(i,:) = immune2;
                    viral_load(i,:) = [0 0 0];
                 end
% %             %cure for strain bc
            elseif logical(population(i,:)) == logical([0 1 1]) 
                viral_load(i,2:3) = viral_load(i,2:3)+1;
                cost= cost+(1*cost_per_treatment);
                 if rand< (b(2))
                    %individual has recovered
                    if immunity(i,:) == [0 0 0]
                        partial_immune = partial_immune+1;
                        infected_count = infected_count-1;
                    end
                    immunity(i,:) = population(i,:);
                    population(i,:) = [0 0 0];
                    color(i,:) = immune2;
                    viral_load(i,:) = [0 0 0];
                 end
           %cure for strain abc
            elseif logical(population(i,:)) == logical([1 1 1]) 
                 viral_load(i,1:3) = viral_load(i,1:3)+1;
                 cost= cost+(1*cost_per_treatment);
                 if rand< (b(1))
                    if immunity(i,:) == [0 0 0]
                        partial_immune = partial_immune+1;
                        infected_count = infected_count-1;
                    end
                    %individual has recovered
                    immunity(i,:) = population(i,:);
                    population(i,:) = [0 0 0];
                    color(i,:) = immune;
                    viral_load(i,:) = [0 0 0];
                 end
            end    

            %death
            %death for strain a
            if color(i,:)== dead
                continue
            end
            if logical(population(i,:)) == logical([1 0 0])
                if rand< d(1)
                    population(i,:) = [0 0 0];
                    color(i,:) = dead;
                    infected_count = infected_count-1;
                    dead_count = dead_count+1;
                end
            %death for strain b
            elseif logical(population(i,:)) == logical([0 1 0])
                 if rand< d(2)
                    color(i,:) = dead;
                    infected_count = infected_count-1;
                    dead_count = dead_count+1;
                 end
            %death for strain c
            elseif logical(population(i,:)) == logical([0 0 1]) 
                 if rand< d(3)
                    population(i,:) = [0 0 0];
                    color(i,:) = dead;
                    infected_count = infected_count-1;
                    dead_count = dead_count+1;
                 end
            %death for strain ab
            elseif logical(population(i,:)) == logical([1 1 0]) 
                 if rand< (d(2))
                    population(i,:) = [0 0 0];
                    color(i,:) = dead;
                    infected_count = infected_count-1;
                    dead_count = dead_count+1;
                 end
            %death for strain ac
            elseif logical(population(i,:)) == logical([1 0 1])
                 if rand< (d(3))
                    population(i,:) = [0 0 0];
                    color(i,:) = dead;
                    infected_count = infected_count-1;
                    dead_count = dead_count+1;
                 end
            %death for strain bc
            elseif logical(population(i,:)) == logical([0 1 1]) 
                 if rand< (d(3))
                    population(i,:) = [0 0 0];
                    color(i,:) = dead;
                    infected_count = infected_count-1;
                    dead_count = dead_count+1;
                 end
            %death for strain abc
            elseif logical(population(i,:) == [1 1 1]) 
                 if rand< (d(3))
                    population(i,:) = [0 0 0];
                    color(i,:) = dead;
                    infected_count = infected_count-1;
                    dead_count = dead_count+1;
                 end
            end    
        for j= 1:N
            if color(j,:) == dead
                continue
            end
            if immunity(j,:) == [1 1 1]
                continue
            else
                dx = x(i) - x(j);
                dy = y(i) - y(j);
                rsquared = dx^2+ dy^2;

                % spreading the infection
                if rsquared<100
                    if immunity(i,:) == [1 1 1]
                        continue
                    elseif abs(population(j,:) + immunity(j,:) - population(i,:)) == [0 0 0]
                        continue
                    %check if they actually transmit if probability is 
                    %greater than infectivity
                    %the condition is set so that they have an 
                    %"a"(variable above) chance to transmit any strain
                    else if (rand<a) 
                            if max(viral_load(i,:))> 7
                                if abs(logical(population(j,:) +immunity(j,:)- population(i,:))) == logical([1 1 1])
                                    %uninfected person came into contact with 
                                    %a completely resitant strain
                                    if population(i,:)== [0 0 0]
                                        infected_count= infected_count+1;
                                    elseif population(j,:)== [0 0 0]
                                        infected_count= infected_count+1;
                                    end
                                    population(i,:) = [1 1 1];
                                    population(j,:) = [1 1 1];
                                    color(i,:) = all_strains;
                                    color(j,:) = all_strains;
                                
                                elseif abs((population(j,:)+immunity(j,:) - population(i,:))) == ([1 0 0])
                                    % person not infected with strain a came into contact 
                                    % with someone who has strain a
                                    if population(i,:)== [0 0 0]
                                        infected_count= infected_count+1;
                                    elseif population(j,:)== [0 0 0]
                                        infected_count= infected_count+1;
                                    end
                                    population(i,1) = 1;
                                    population(j,1) = 1;
                                    color(i,:) = one_strain;
                                    color(j,:) = one_strain;
                                elseif abs((population(j,:) +immunity(j,:) - population(i,:))) == ([0 1 0])
                                    % person not infected with strain b came into contact 
                                    % with someone who has strain b
                                    if population(i,:)== [0 0 0]
                                        infected_count= infected_count+1;
                                    elseif population(j,:)== [0 0 0]
                                        infected_count= infected_count+1;
                                    end
                                    population(i,2) = 1;
                                    population(j,2) = 1;
                                    color(i,:) = one_strain;
                                    color(j,:) = one_strain;                         
                                elseif abs((population(j,:)+immunity(j,:) - population(i,:))) ==([0 0 1])
                                    % person not infected with strain c came into contact 
                                    % with someone who has strain c
                                    if population(i,:)== [0 0 0]
                                        infected_count= infected_count+1;
                                    elseif population(j,:)== [0 0 0]
                                        infected_count= infected_count+1;
                                    end
                                    population(i,3) = 1;
                                    population(j,3) = 1;
                                    color(i,:) = one_strain;
                                    color(j,:) = one_strain;
                                elseif abs((population(j,:) +immunity(j,:)- population(i,:))) == ([1 1 0])
                                    % person not infected with strain ab came into contact 
                                    % with someone who has strain ab
                                    if population(i,:)== [0 0 0]
                                        infected_count= infected_count+1;
                                    elseif population(j,:)== [0 0 0]
                                        infected_count= infected_count+1;
                                    end
                                    population(i,1) = 1;
                                    population(j,1) = 1;
                                    population(i,2) = 1;
                                    population(j,2) = 1;
                                    color(i,:) = two_strains;
                                    color(j,:) = two_strains;
                                elseif abs((population(j,:) +immunity(j,:)- population(i,:))) == ([1 0 1])
                                    % person not infected with strain ac came into contact 
                                    % with someone who has strain ac
                                    if population(i,:)== [0 0 0]
                                        infected_count= infected_count+1;
                                    elseif population(j,:)== [0 0 0]
                                        infected_count= infected_count+1;
                                    end
                                    population(i,1) = 1;
                                    population(j,1) = 1;
                                    population(i,3) = 1;
                                    population(j,3) = 1;
                                    color(i,:) = two_strains;
                                    color(j,:) = two_strains;
                                elseif abs((population(j,:) +immunity(j,:)- population(i,:))) == ([0 1 1])
                                    % person not infected with strain abc came into contact 
                                    % with someone who has strain bc
                                    if population(i,:)== [0 0 0]
                                        infected_count= infected_count+1;
                                    elseif population(j,:)== [0 0 0]
                                        infected_count= infected_count+1;
                                    end
                                    population(i,2) = 1;
                                    population(j,2) = 1;
                                    population(i,3) = 1;
                                    population(j,3) = 1;
                                    color(i,:) = two_strains;
                                    color(j,:) = two_strains;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    end
end
    %plotting
    %plotting the vector field 
    quiver(xx,yy,uu,vv, 'FaceColor', 'b')
    hold on
    
    %plotting the individuals
    scatter(x,y, 30, color, 'filled') 
        %plotting the radius of infectivity
        for i= 1:N
            if color(i,:) == dead
                continue
            elseif immunity(i,:) == [1 1 1]
                continue
            end
            t = 0 : 0.5 : 2 * pi;
            r = 10;
            x_circ = x(i)+(r*cos(t));
            y_circ = y(i)+(r*sin(t));
            if alpha(i)>=7
                alpha(i) = 6;
            end
            patch( x_circ,y_circ,color(i,:),'FaceAlpha',alpha(i)/7,'EdgeColor', 'none');
        end
    hold off
    healthy_count = N-(infected_count+ dead_count + immune_count+partial_immune);
%     whitebg(1,'k')
    set(gca,'color','black')
    set(gcf,'color','black')
    axis equal
    %displaying counts
    title({['Susceptible Count: ' ,num2str(healthy_count)];['Partially Immune Count: ' , num2str(partial_immune)];['Fully Immune Count: ' , num2str(immune_count)];['Infected Count :', num2str(infected_count)]; ['Death Count: ', num2str(dead_count)]; ['Cost of Treatment: ', num2str(cost)]});
    axis([min-40, maximum+40, min-40, maximum+40])
    drawnow
    axis manual
    
end

