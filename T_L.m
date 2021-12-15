% TEACHING-LEARNING ALGORITHM

clc;
clear all;
close all;

% Choose the function (Ex. 0, 1, 2... etc.)
f_objetivo = 1;

switch f_objetivo
    case 0
        %Griewank
       f = @(x,y) ((x.^2 + y.^2)/4000)-cos(x).*cos(y/sqrt(2)) + 1;
       U = [10 10];
       L = [-10 -10];
       %Minimum = 0; x= 0, y= 0
       
    case 1
        %Rastrigin
       f = @(x,y)  10*2 + x.^2 - 10 .*cos(pi*x)+ y.^2 - 10 .* cos(pi*y);
       U = [5 5];
       L = [-5 -5];
       %Minimum = 0; x = 0, y = 0 
       
      
    case 2
         %DropWave
    f = @(x,y) - ((1 + cos(12*sqrt(x.^2+y.^2))) ./ (0.5 *(x.^2+y.^2) + 2));
    U = [2 2];
    L = [-2 -2];
       %Minimum = -1; x = 0, y =0
      
    case 3
        %Esphere
        f = @(x,y) (x-2).^2 + (y-2).^2;
        U = [5 5];
        L = [-5 -5];
        %Minimum = 0; x = 2, y = 2 
       
    otherwise
        disp("Choose a valid value of the function")
        return
end

% % % % % % % % % % % % % % % % % % % % % % % % % % 
%contour of the function
[X,Y] = meshgrid(L(1):0.2:U(1), L(2):0.2:U(2));
Z = f(X,Y);
grid on;
contour(X,Y,Z,40);
hold on
% % % % % % % % % % % % % % % % % % % % % % % % % % 


D = 2;
N =35;
G = 150;

%List to find the best solution
lista = zeros(1,G);

pob = zeros(D,N);
teacher = zeros(D,1);

fitness = zeros(1,N);

% Initialization
r = rand(D,N);
%Asigning random values to the population
for i=1:D
    pob(i,:) = L(i) + (U(i) - L(i)) .* r(i,:);
end

%Calculating fitness
fitness = f(pob(1,:),pob(2,:));


for i=1:G

 %Aplying both phases of the algorithm for each individual
 for k=1:N
      % Defining the teacher
    [~,t] = min(fitness);
    teacher = pob(:,t);
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
     % TEACHING PHASE
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    %CALCULATING LEARNING FACTOR
    Tf = randi([1,2]);
    
    x_mean = zeros(1,2);
    c = zeros(2,1);
    
    for j=1:D
      % CALCULATING AVERAGE OF EACH VARIABLE OF THE POPULATION
        x_mean = mean(pob(j,:));

        r =  rand();
        % NEW VARIABLES FOR A NEW CANDIDATE SOLUTION

        c(j) = pob(j,k) + r * (teacher(j) - Tf * x_mean);

    end
    
    %If the new solution is doing betetr
    if f(c(1),c(2)) < f(pob(1,k),pob(2,k))
        pob(:,k) = c;
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
     % lEARNING PHASE
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    %Choosing a different individual
    n = k;
    while n == k
        n =  randi([1,N]);
    end
    
    % Ealuating each solution
    a = fitness(k);
    b = fitness(n);
    
    % if actual individual has a better solution
    if a < b
        %Adjusting direction of the candidate function to the better
        %individual
        for j=1:D
            r = rand();
            c(j) = pob(j,k) + r * (pob(j,k) - pob(j,n));
        end
        % else
    else
       %Adjusting direction of the candidate function to the better
        for j=1:D
            r = rand();
            c(j) = pob(j,k) + r * (pob(j,n) - pob(j,k));
        end
        
    end
    % If the candidate solution is better than the actual individual
    if f(c(1),c(2)) < f(pob(1,k),pob(2,k))
        pob(:,k) = c;
    end
     fitness = f(pob(1,:),pob(2,:));
    
 end
 
 %Calculating fitness to choose better solution
 fitness = f(pob(1,:),pob(2,:));
 [value,v] = min(fitness);
 %Add better solution to the list
 lista(i) = value;
 
 %Graphics
 axis([L(1) U(1) L(2) U(2)]);
 h = plot(pob(1,:),pob(2,:), 'rx');
 pause(0.05);
 delete(h);

end

h = plot(pob(1,:),pob(2,:), 'rx');

figure
plot(lista);

%The better solution would be the last one that we added to the list
best = pob(:,v)
f(best(1),best(2))






