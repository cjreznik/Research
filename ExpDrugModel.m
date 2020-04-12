function xp = ExpDrugModel(t,x,r_s,r_r,g,lambda_s,lambda_r,gamma)

%x(1) is Sensitive cells, x(2) is Resistant cells, x(3) is Drug

xp = zeros(3,1);



xp(1) = r_s*x(1)-g*x(1)-(lambda_s)*x(3)*x(1); %sensitive
xp(2) = r_r*x(2)+g*x(1)-(lambda_r)*x(3)*x(2); %resistant
xp(3) = -(gamma)*x(3); %drug
%Volume = xp(1) + xp(2);






%How to run:
 %[t,x] = ode45('DrugModel',[t0,tf],[x(1)0,x(2)0,x(3)0]);
 
% How to plot:
% figure; plot(t,x(:,1),'b'); hold on; plot(t,x(:,2),'r'); hold off
%Plot w Volume
% figure; plot(t,x(:,1),'b'); hold on; plot(t,x(:,2),'r'); hold on; plot(t,(x(:,1)+x(:,2)),'g');hold off
 