% trim the reachable workspace around geometric workspace
clear all;close all;clc
load('April4th_res_10.mat')
count = 0;
Nvalid = [];
cvalid = [];
for i=1:size(Nf,3)
if(equilibriumtest(Nf(:,:,i))==1)
    count=count+1;
    Nvalid(:,:,end+1) = Nf(:,:,i);
    N = Nf(:,:,i);
    cvalid = [cvalid double((1/3).*(N(:,4)+N(:,5)+N(:,6)))];
end
end
fprintf("End of program! No of valid points : ")+string(count);
save('valid_10_res_equilibrium_points_actualprism.mat');
scatter3(cvalid(1,:),cvalid(2,:),cvalid(3,:),'filled')

function flag = equilibriumtest(N)
Cb = [-eye(3) eye(3)];
Cs = [0 0 0 0 -1 1; % 5-6
      0 0 0 1 0 -1; % 6-4
      0 0 0 -1 1 0; % 4-5
      0 -1 0 0 0 1; % 2-6
      0 0 -1 1 0 0; % 3-4
      -1 0 0 0 1 0]; % 1-5
  nb = size(Cb,1);
  ns = size(Cs,1);
  S = N*Cs';
  B = N*Cb';
W = [0 0 0 0 0 0;
     0 0 0 0 0 0;
     0 0 -0 -5 -5 -5];
% split into pinned and unpinned
Nx = N(:,4:end);
N0 = N(:,1:3);
Cbx = Cb(:,4:end); 
Cb0 = Cb(:,1:3);
Csx = Cs(:,4:end);
Cs0 = Cs(:,1:3);
Wx = W(:,4:6);
W0 = W(:,1:3);
EYE = eye(size(Csx,2));
Ax=[];Bx=[];

for i=1:size(EYE,2)
    Ax = [Ax ;S*diag(Csx*EYE(:,i)) -B*diag(Cbx*EYE(:,i))];
    Bx = [Bx ;Wx*EYE(:,i)];
end
f = [diag(S'*S)' diag(B'*B)'];
options = optimoptions('linprog','Algorithm','interior-point','Display','off');
[x,fval,flag] = linprog(f',[],[],Ax,Bx,zeros(ns+nb,1),[],options);
if flag==1
    fprintf('EQUILIBRIUM \n')
end
end