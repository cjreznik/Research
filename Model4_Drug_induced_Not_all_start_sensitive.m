clear all; 
close all;
rng(3);
%%Clock function and saves best fit each time
cl = clock; clN = 0;
for ii = 2:5
    clN = floor(100*clN + cl(ii));
end
path = ['Best_Fit_' , num2str(clN)];
if exist(path, 'dir') ~= 7
    mkdir(path)
end

%% Variables
filename = 'CTXData.xlsx';
A = xlsread(filename);
A(A == 0) = NaN; %first row is not numbers
%stop
A(1,1) = 0; %fills in first cell as Day 0
%stop
z = [];   %error sum
 newParam_r_s = {};
 newParam_r_r = {};
 newParam_g = {};
 newParam_lambda_s = {};
 newParam_gamma = {};
 newParam_lambda_r = {};
y0 = {};
SSE = {}; %cell arrays for newParam, initial cond & SSE variables
Volume = {}; 
Days = {};
%soln_all = {};
solution_all = {};
time_all = {}; 

%% Reading in data for Mouse m
num_mice = size(A,2);%2 or 3?

Best_r_s = zeros(num_mice-1,1);
Best_r_r = zeros(num_mice-1,1);
Best_g = zeros(num_mice-1,1);
Best_lambda_s = zeros(num_mice-1,1);
Best_gamma = zeros(num_mice-1,1);
Best_lambda_r = zeros(num_mice-1,1);
Best_SSE = zeros(num_mice-1,1);
Best_y0 = zeros(num_mice-1,1);
Best_y0r = zeros(num_mice-1,1);

for m = 2:num_mice
   
   Volume_orig = A(:,m);
   Days_orig = A(:,1);
   Volume{m-1}(:) = rmmissing(Volume_orig); %removes any entry that contains missing data
   length(Volume{m-1}(:));
   Days{m-1}(:) = Days_orig(1:length(Volume{m-1}(:))); 
 
end

%%%%%%%%%%%%%%QMC test code%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%initialize min and max values   
min_r_s = 0; max_r_s = .1;
min_r_r = 0; max_r_r = .1;
min_g = 0; max_g = .1;
min_lambda_s = 0; max_lambda_s = .1;
min_gamma = 0; max_gamma = .1;
min_lambda_r = 0; max_lambda_r = .1;

for m = 2:num_mice
    % Scale Volume to be between 50% and 150%
    
    min_y0 = 0; %0.5*Volume{m-1}(1);
    max_y0 = Volume{m-1}(1); %1.5*Volume{m-1}(1);
    
    u = 1;
    i = 1;
    %y0{m-1}(i) =  Volume{m-1}(1);
    n_pts = 10000 
    n_skip = 1000;
    n_leap = 0;
    %create sobol points and scale them
    uniform_sobol = sobolset(8,'Skip',n_skip,'Leap',n_leap)
    uniform_sobol = net(uniform_sobol,n_pts);


    uniform_sobol_scaled(:,1) = (max_r_s-min_r_s)*uniform_sobol(:,1) + min_r_s; %r_s
    uniform_sobol_scaled(:,2) = (max_r_r-min_r_r)*uniform_sobol(:,2) + min_r_r; %r_r
    uniform_sobol_scaled(:,3) = (max_g-min_g)*uniform_sobol(:,3) + min_g; %g
    uniform_sobol_scaled(:,4) = (max_lambda_s-min_lambda_s)*uniform_sobol(:,4) + min_lambda_s; %lambda_s
    uniform_sobol_scaled(:,5) = (max_lambda_r-min_lambda_r)*uniform_sobol(:,5) + min_lambda_r; %lambda_r
    uniform_sobol_scaled(:,6) = (max_gamma-min_gamma)*uniform_sobol(:,6) + min_gamma;  %gamma
    uniform_sobol_scaled(:,7) = (max_y0-min_y0)*uniform_sobol(:,7) + min_y0; % y0
    uniform_sobol_scaled(:,8) = (max_y0-min_y0)*uniform_sobol(:,8) + min_y0;
    SSE_QMC = [];  
    

    count = 0;
    for u = 1:n_pts % do what we did below: solve DE using these parameters and get SSE
       
        if ((uniform_sobol_scaled(u,4) > uniform_sobol_scaled(u,5))&&(uniform_sobol_scaled(u,1)>=uniform_sobol_scaled(u,2)))
           count = count + 1;
        
           if mod(u,100) == 0
            fprintf('Up to %d\n',u); 
           end    
        
        dosage =  ceil((Days{1,m-1}(end)/7)); %calculates dosage for each mouse based on days

         d = 2; %initialize d
        %create first soln vector
        soln{d-1} = ode23s(@(t,x)  ExpDrugModel_Drug_causes_resistance(t,x,uniform_sobol_scaled(u,1),uniform_sobol_scaled(u,2),uniform_sobol_scaled(u,3), uniform_sobol_scaled(u,4),uniform_sobol_scaled(u,5),uniform_sobol_scaled(u,6)),0:7,[uniform_sobol_scaled(u,7),uniform_sobol_scaled(u,8),5]);

        t1 = 7;
        t2 = 14;

       %loop through d times to create d sol vecs
        for d = 2:dosage

          soln{d} = ode23s(@(t,x)  ExpDrugModel_Drug_causes_resistance(t,x,uniform_sobol_scaled(u,1),uniform_sobol_scaled(u,2),uniform_sobol_scaled(u,3), uniform_sobol_scaled(u,4),uniform_sobol_scaled(u,5),uniform_sobol_scaled(u,6)),t1:t2,[soln{d-1}.y(1,end),soln{d-1}.y(2,end),soln{d-1}.y(3,end)+5]);

          t1 = t1 + 7;

          t2 = t2 + 7;

          d = d + 1;

        end
    %create initial time and solution vecs to concatenate
        d = 1;
        time = soln{d}.x;
        solution = [soln{d}.y(1,:)+soln{d}.y(2,:)];
        d = d+1; 

        %concatenating%
        for d = 2:dosage
          temp = [soln{d}.y(1,2:end)+soln{d}.y(2,2:end)]; %concatenate the 2 y soln vecs together
          time = [time soln{d}.x(2:end);]; %concatenate time
          solution = [solution temp]; % solution vec being concatenated
        end
        %concatenated vectors of solution and time
        solution_all{m-1} = solution; 
        time_all{m-1} = time; 

        d = 1;
        t1 = 7;
        t2 = 7;

       
        %initialSSE = 0; 
        newz = 0;
        soln_index = [];
        
        %instead of using "deval" index to nearest point on soultion to get
        %SSE%
        for r = 1:length(Days{m-1})

        %       timevec = time_all{m-1};
        %       solvec = solution_all{m-1};

              [c , index] = min(abs(time-Days{m-1}(r)));
              soln_index(r)= solution(index);
              error = (solution(index) - Volume{m-1}(r))^2;
              newz = newz + error;
              %tempindx{m-1}(r) = index;    

        end
          
          SSE_QMC(u) = newz;
        else
            
            SSE_QMC(u) = 10^10;
%           fprintf('r_r = %f, r_s = %f, g = %f, lambda_r = %f, lambda_s = %f, gamma = %f, SSE = %f\n',...
%             uniform_sobol_scaled(u,1), uniform_sobol_scaled(u,2),uniform_sobol_scaled(u,3),...
%             uniform_sobol_scaled(u,4),uniform_sobol_scaled(u,5),uniform_sobol_scaled(u,6), SSE_QMC(u)); 
       end
       end

    [c, index] = min(SSE_QMC)

  %create a loop that takes the uniform sobol pts at the min SSE index and save it          
                    
   i = 1;        
  
   tspan = 0:Days{m-1}(end);
   newParam_r_s{m-1}(i) =  uniform_sobol_scaled(index,1); 
   newParam_r_r{m-1}(i) = uniform_sobol_scaled(index,2);
   newParam_g{m-1}(i) = uniform_sobol_scaled(index,3);
   newParam_lambda_s{m-1}(i) = uniform_sobol_scaled(index,4);
   newParam_lambda_r{m-1}(i) = uniform_sobol_scaled(index,5); 
   newParam_gamma{m-1}(i) = uniform_sobol_scaled(index,6);
   y0{m-1}(i)=  uniform_sobol_scaled(index,7); 
   y0r{m-1}(i) = uniform_sobol_scaled(index,8);
%    % plot model at these points 
%     IC = [y0{m-1}(1),0,5];
%    tspan = 0:Days{m-1}(end);  
%    sol = ode23(@(t,x) ExpDrugModel(t,x,newParam_r_s{m-1}(i),newParam_r_r{m-1}(i),newParam_g{m-1}(i),newParam_lambda_s{m-1}(i),newParam_gamma{m-1}(i),newParam_lambda_r{m-1}(i)),tspan,IC);
%    
%       
%for m = 2:3 %num_mice
    sol = {};
    solution = [];
    time = [];
    i = 1;    
   
%%%%%%%%%%%%%%Drug Dosing with best QMC points%%%%%%%%%%%%%%%

   dosage =  ceil((Days{1,m-1}(end)/7)); %calculates dosage for each mouse based on days
  
   d = 2; %initialize d
   
   %create initial solution for sol(d)
   
   
   soln{d-1} = ode23s(@(t,x)  ExpDrugModel_Drug_causes_resistance(t,x,newParam_r_s{m-1}(i),newParam_r_r{m-1}(i),newParam_g{m-1}(i), newParam_lambda_s{m-1}(i),newParam_lambda_r{m-1}(i),newParam_gamma{m-1}(i)),0:7,[y0{m-1}(i),y0r{m-1}(i),5]);
   %fprintf('Ok 1\n'); 
 
   t1 = 7;
   t2 = 14;
   
   %loop through d times to create d sol vecs
        for d = 2:dosage
              soln{d} = ode23s(@(t,x) ExpDrugModel_Drug_causes_resistance(t,x,newParam_r_s{m-1}(i),newParam_r_r{m-1}(i),newParam_g{m-1}(i), newParam_lambda_s{m-1}(i),newParam_lambda_r{m-1}(i),newParam_gamma{m-1}(i)),t1:t2,[soln{d-1}.y(1,end),soln{d-1}.y(2,end),soln{d-1}.y(3,end)+5]);

              t1 = t1 + 7;

              t2 = t2 + 7;

              d = d + 1;

        end

%     %test plot%
    d = 1;
    figure;
    plot(Days{m-1}(:),Volume{m-1}(:),'o'); hold on;
    for d = 1:dosage

        plot(soln{d}.x, soln{d}.y(1,:)+soln{d}.y(2,:)+soln{d}.y(3,:));
        hold on; 

    end     

    d = 1;

    time = soln{d}.x;
    solution = [soln{d}.y(1,:)+soln{d}.y(2,:)];
    d = d+1; 

    %concatenating%
      for d = 2:dosage
          temp = [soln{d}.y(1,2:end)+soln{d}.y(2,2:end)];
          time = [time soln{d}.x(2:end);];
          solution = [solution temp];
          d = d + 1;

      end

     solution_all{m-1} = solution; 
     time_all{m-1} = time; 

      d = 1;
      t1 = 7;
      t2 = 14;

      %deval%
      initialSSE = 0; 
      newz = 0;
      for r = 1:length(Days{m-1})
      
%       timevec = time_all{m-1};
%       solvec = solution_all{m-1};
      
      [c , index] = min(abs(time-Days{m-1}(r)));
      error = initialSSE + (solution(index) - Volume{m-1}(r))^2;
      newz = newz + error;
      %tempindx{m-1}(r) = index;    
      
  end
  
  
  i = 1;


 %newz = 0;
   SSE{m-1}(i) = newz;   
    convornot = 0;
    check = 0;  

%for m = 2:3
%stop 


%     %while (convornot < 5) %while loop with convergence criteria
    for j = 1:10000%000 %(used to test if loop works)
        %fprintf('At step %d\n',j);
        %% step 1 : generate random r in (-1,1) 
        a = -1;
        b = 1;
        r = ((b-a)*(rand)+a)*0.01; 
        s = ((b-a)*(rand)+a)*0.01;
        q = ((b-a)*(rand)+a)*0.01;
        s2 = ((b-a)*(rand)+a)*0.01;
        r2 = ((b-a)*(rand)+a)*0.01;
        q2 = ((b-a)*(rand)+a)*0.01;
        v = ((b-a)*(rand)+a)*10;
        v2 = ((b-a)*(rand)+a)*10;
        
%         i
%         newParam_r_s{m-1}(i)
       
        tempnewparam_r_s = newParam_r_s{m-1}(i)+r;
        tempnewparamr_r = newParam_r_r{m-1}(i) + s;
        tempnewparam_g = newParam_g{m-1}(i)+q;
        tempnewparam_lambda_s = newParam_lambda_s{m-1}(i)+s2;
        tempnewparam_gamma = newParam_gamma{m-1}(i) + r2;
        tempnewparam_lambda_r  = newParam_lambda_r{m-1}(i) + q2;
        tempnewy0 = y0{m-1}(i) + v;
        tempnewy0r = y0r{m-1}(i) + v2;
   
    
 
        while tempnewy0 < 0 || tempnewparam_lambda_s < 0 || tempnewparam_g < 0 ||  tempnewparam_gamma < 0  || tempnewparam_lambda_r <0 || tempnewparam_lambda_r > tempnewparam_lambda_s

                    a = -1;
                    b = 1;
                    r = ((b-a)*(rand)+a)*0.01;
                    s = ((b-a)*(rand)+a)*0.01;
                    q = ((b-a)*(rand)+a)*0.01;
                    s2 = ((b-a)*(rand)+a)*0.01;
                    r2 = ((b-a)*(rand)+a)*0.01;
                    q2 = ((b-a)*(rand)+a)*0.01;
                    v = ((b-a)*(rand)+a)*10;
                    v2 = ((b-a)*(rand)+a)*10;
                    tempnewparam_r_s = newParam_r_s{m-1}(i)+r;
                    tempnewparamr_r = newParam_r_r{m-1}(i) + s;
                    tempnewparam_g = newParam_g{m-1}(i)+q;
                    tempnewparam_lambda_s = newParam_lambda_s{m-1}(i)+s2;
                    tempnewparam_lambda_r  = newParam_lambda_r{m-1}(i) + q2;
                    tempnewparam_gamma = newParam_gamma{m-1}(i) + r2;
                    tempnewy0 = y0{m-1}(i) + v;
                    tempnewy0r = y0r{m-1}(i) + v2;
   
                end

    
           %fprintf('tempnewparam_r_s = %f\n tempnewparamr_r= %f\n tempnewparamr_g = %f\n tempnewparam_lambda_s= %f\n tempnewparam_gamma= %f\n tempnewparam_lambda_r= %f\n tempnewy0= %f\n',tempnewparam_r_s,tempnewparamr_r, tempnewparam_g,tempnewparam_lambda_s, tempnewparam_gamma,tempnewparam_lambda_r, tempnewy0  ) 


            % step 2 :solve DE using ode23 and dosing

            dosage =  ceil((Days{1,m-1}(end)/7));
             d = 2;
             sol = {};
             sol{d-1} = ode23s(@(t,x)  ExpDrugModel_Drug_causes_resistance(t,x,tempnewparam_r_s,tempnewparamr_r,tempnewparam_g, tempnewparam_lambda_s,tempnewparam_lambda_r,tempnewparam_gamma),0:7,[tempnewy0,tempnewy0r,5]);
             t1 = 7;
             t2 = 14;
             
            %loop through d times to create d sol vecs
                for d = 2:dosage

                     sol{d} = ode23s(@(t,x)  ExpDrugModel_Drug_causes_resistance(t,x,tempnewparam_r_s,tempnewparamr_r,tempnewparam_g, tempnewparam_lambda_s,tempnewparam_lambda_r,tempnewparam_gamma),t1:t2,[sol{d-1}.y(1,end),sol{d-1}.y(2,end),sol{d-1}.y(3,end)+5]);

                    %sol{d} = ode23s(@(t,x) ExpDrugModel(t,x,newParam_r_s{m-1}(i),newParam_r_r{m-1}(i),newParam_g{m-1}(i), newParam_lambda_s{m-1}(i),newParam_lambda_r{m-1}(i),newParam_gamma{m-1}(i)),t1:t2,[sol{d-1}.y(1,end),sol{d-1}.y(2,end),sol{d-1}.y(3,end)+5]);


                    t1 = t1 + 7;

                    t2 = t2 + 7;

                    d = d + 1;

                end
        
            d = 1;
            time = []; 
            time = sol{d}.x;
            solution = []; 
            solution = [sol{d}.y(1,:)+sol{d}.y(2,:)];
            d = d+1; 

    %concatenating?%
      for d = 2:dosage
          temp = [sol{d}.y(1,2:end)+sol{d}.y(2,2:end)];
          time = [time sol{d}.x(2:end);];
          solution = [solution temp];
          d = d + 1;

      end

        
        %sol = ode23s(@(t,x) ExpDrugModel(t,x,tempnewparam_r_s,tempnewparamr_r,tempnewparam_g, tempnewparam_lambda_s,tempnewparam_lambda_r,tempnewparam_gamma),0:7,[tempnewy0,0,5]);
         
        %y = deval(sol,Days{m-1},1)';
            newz = 0; 
    
              d = 1;
              t1 = 7;
              t2 = 14;

          %deval%
          initialSSE = 0; 
          newz = 0;
          for r = 1:length(Days{m-1})

        %       timevec = time_all{m-1};
        %       solvec = solution_all{m-1};

              [c , index] = min(abs(time-Days{m-1}(r)));
              error = initialSSE + (solution(index) - Volume{m-1}(r))^2;
              newz = newz + error;
              %tempindx{m-1}(r) = index;    

          end

        % step 3: Calculate goodness of fit
%             for w = 1:size(Days{m-1}(:), 1) 
%                 error = (Volume{m-1}(w) - y(w))^2;
%                 newz = newz + error;
%             end
%             
        b = 1000; %beta value
        d = newz - SSE{m-1}(i); %delta value
                % step 4: prob of acceptance
            if d < 0
               probability = 1;
            else 
               probability = exp(-b*d);
            end
        %fprintf('Current SSE = %f gives prob = %f\n',newz,probability); 
        % step 5: Randomly perturbate to get new r and new fit
        newr = rand;
            if (newr <= probability)
                fprintf('New parameters have SSE = %f, compared to old SSE = %f\n',newz, SSE{m-1}(i)); 

            i = i + 1;
            SSE{m-1}(i) = newz;
            newParam_r_s{m-1}(i) = tempnewparam_r_s;
            newParam_r_r{m-1}(i) = tempnewparamr_r;
            newParam_g{m-1}(i) = tempnewparam_g;
            newParam_lambda_s{m-1}(i) = tempnewparam_lambda_s;
            newParam_gamma{m-1}(i) = tempnewparam_gamma;
            newParam_lambda_r{m-1}(i) = tempnewparam_lambda_r;
            y0{m-1}(i) = tempnewy0; 
            y0r{m-1}(i) = tempnewy0r;

        
       % fprintf('\tAccepting these parameters, tempnewparam_r_s = %f\n ,tempnewparamr_r= %f\n, tempnewparam_g = %f\n, tempnewparam_lambda_s= %f\n, tempnewparam_gamma= %f\n, tempnewparam_lambda_r= %f\n ,tempnewy0= %f\n, d = %f, and probability = %f\n', tempnewparam_r_s, tempnewparamr_r, tempnewparam_g, tempnewparam_lambda_s , tempnewparam_gamma, tempnewparam_lambda_r, tempnewy0, d, probability)
        
            r_s_error = abs(newParam_r_s{m-1}(i) - newParam_r_s{m-1}(i-1))/abs(newParam_r_s{m-1}(i-1));
            r_r_error = abs(newParam_r_r{m-1}(i) - newParam_r_r{m-1}(i-1))/abs(newParam_r_r{m-1}(i-1));
            g_error = abs(newParam_g{m-1}(i) - newParam_g{m-1}(i-1))/abs(newParam_g{m-1}(i-1));
            lambda_s_error = abs(newParam_lambda_s{m-1}(i) - newParam_lambda_s{m-1}(i-1))/abs(newParam_lambda_s{m-1}(i-1));
            gamma_error = abs(newParam_gamma{m-1}(i) - newParam_gamma{m-1}(i-1))/abs(newParam_gamma{m-1}(i-1));
            lambda_r_error = abs(newParam_lambda_r{m-1}(i) - newParam_lambda_r{m-1}(i-1))/abs(newParam_lambda_r{m-1}(i-1));
            y0_error = abs(y0{m-1}(i) - y0{m-1}(i-1))/abs(y0{m-1}(i-1));
            y0r_error = abs(y0r{m-1}(i) - y0r{m-1}(i-1))/abs(y0r{m-1}(i-1));
                if r_s_error < .01 && r_r_error < .01 && g_error < .01 && lambda_s_error < .01 && gamma_error < .01 && lambda_r_error <.01 && y0_error < .01 && y0r_error < .01
                    convornot = convornot + 1;
                else
                    convornot = 0;

                end
         
            else
       % fprintf('\tNot accepting tempnewparam_r_s = %f\n , tempnewparamr_r= %f\n , tempnewparam_g = %f\n , tempnewparam_lambda_s= %f\n , tempnewparam_gamma= %f\n , tempnewparam_lambda_r= %f\n , tempnewy0= %f\n , d = %f , and probability = %f\n', tempnewparam_r_s, tempnewparamr_r, tempnewparam_g, tempnewparam_lambda_s , tempnewparam_gamma, tempnewparam_lambda_r, tempnewy0, d, probability);
        
            end
    
   %check = check + 1;
    
       end


  %end

  
%for m = 14:15
  
    
    r_s = newParam_r_s{m-1}(end);
    r_r = newParam_r_r{m-1}(end);
    g = newParam_g{m-1}(end);
    lambda_s = newParam_lambda_s{m-1}(end);
    gamma = newParam_gamma{m-1}(end);
    lambda_r = newParam_lambda_r{m-1}(end);
    IC = [y0{m-1}(1),y0r{m-1}(1),5];
    tspan = 0:Days{m-1}(end);
    %Days2 = Days{m-1}(end);
    %Volume2 = Volume{m-1}(end);
     
    d = 2;
    
    dosage = ceil((Days{1,m-1}(end)/7));
    
    sol{d-1} = ode23s(@(t,x)  ExpDrugModel_Drug_causes_resistance(t,x,newParam_r_s{m-1}(end),newParam_r_r{m-1}(end),newParam_g{m-1}(end), newParam_lambda_s{m-1}(end),newParam_lambda_r{m-1}(end),newParam_gamma{m-1}(end)),0:7,[y0{m-1}(1),y0r{m-1}(1),5]);
    
    
    t1 = 7;
    t2 = 14;
    for d = 2:dosage
    
     sol{d} = ode23s(@(t,x)  ExpDrugModel_Drug_causes_resistance(t,x,newParam_r_s{m-1}(end),newParam_r_r{m-1}(end),newParam_g{m-1}(end), newParam_lambda_s{m-1}(end),newParam_lambda_r{m-1}(end),newParam_gamma{m-1}(end)),t1:t2,[soln{d-1}.y(1,end),soln{d-1}.y(2,end),soln{d-1}.y(3,end)+5]);
      
      t1 = t1 + 7;
      
      t2 = t2 + 7;
      
      d = d + 1;
    
    end
    
    
     d = 1;
            time = []; 
            time = sol{d}.x;
            solution = []; 
            solution = [sol{d}.y(1,:)+sol{d}.y(2,:)];
            d = d+1; 
    
    
    figure;
     for d = 2:dosage
          temp = [sol{d}.y(1,2:end)+sol{d}.y(2,2:end)];
          time = [time sol{d}.x(2:end);];
          solution = [solution temp];
          d = d + 1;

      end
    
    plot(time,solution); hold on ; plot(Days{m-1}(:),Volume{m-1}(:),'o'); hold off
  
    Best_r_s(m-1,1) = newParam_r_s{m-1}(end);
    Best_r_r(m-1,1) = newParam_r_r{m-1}(end);
    Best_g(m-1,1) = newParam_g{m-1}(end);
    Best_lambda_s(m-1,1) = newParam_lambda_s{m-1}(end);
    Best_gamma(m-1,1) = newParam_gamma{m-1}(end);
    Best_lambda_r(m-1,1) = newParam_lambda_r{m-1}(end);
    Best_SSE(m-1,1) = SSE{m-1}(end);
    Best_y0(m-1,1) = y0{m-1}(end);
    Best_y0r(m-1,1) = y0r{m-1}(end);

    fprintf('For mouse %d, Best_r_s = %f, Best_r_r = %f, Best_g = %f, Best_lambda_s = %f , Best_gamma = %f, Best_lambda_r = %f, Best_SSE = %f, Best_y0 = %f, Best_y0r = %f',m-1, Best_r_s(m-1), Best_r_r(m-1), Best_g(m-1), Best_lambda_s(m-1), Best_gamma(m-1), Best_lambda_r(m-1), Best_SSE(m-1),Best_y0(m-1),Best_y0r(m-1))
      
        
        
    fname_fig = [path '/mouse' num2str(m-1)];        
     saveas(gcf,[fname_fig,'.fig'])
     saveas(gcf,[fname_fig,'.png']);

     
     
     
end
    
     filename = [path '/BestParams.mat']; 
  
    save(filename, 'Best_r_s', 'Best_r_r' , 'Best_g',  'Best_lambda_s' , 'Best_gamma' , 'Best_lambda_r','Best_SSE','Best_y0','Best_y0r') 

     


