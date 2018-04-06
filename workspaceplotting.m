clear all;clc;close all;
d2r = pi/180;

fprintf('Program initialized\n');
tic
%% initialize values and compute the workspace

dth=30;
res = 50;
th1v = linspace(-10,60,res);th2v=linspace(-10,60,res);th3v=linspace(-10,60,res);
ph1v = linspace(10,90,res);ph2v=linspace(10,90,res);ph3v=linspace(10,90,res);
th1v=d2r.*th1v;th2v=d2r.*th2v;th3v=d2r.*th3v;
ph1v=d2r.*ph1v;ph2v=d2r.*ph2v;ph3v=d2r.*ph3v;

bar_l = 10;
edge_l = 7.05;

j=1;

% [x1,x2,x3,x4,x5,x6] = getpoints(10*d2r,10*d2r,10*d2r,60*d2r,60*d2r,60*d2r,2,2);
% N = [x1,x2,x3,x4,x5,x6]

% figure(1)
% [x1p,x2p,x3p,x4p,x5p,x6p] = getpoints(th1v,th2v,th3v,ph1v,ph2v,ph3v,bar_l,edge_l);
% line([x1p(1) x4p(1)],[x1p(2) x4p(2)],[x1p(3) x4p(3)],'Color','black');hold on;
% line([x2p(1) x5p(1)],[x2p(2) x5p(2)],[x2p(3) x5p(3)],'Color','black');hold on;
% line([x3p(1) x6p(1)],[x3p(2) x6p(2)],[x3p(3) x6p(3)],'Color','black');hold on;
%                 
% axis([-2 2 -2 2 0 3])

c=[];Nf = [];
% the great search starts
for i1=1:1:numel(th1v)
    for i2 = 1:1:numel(th2v)
        for i3 = 1:1:numel(th3v)
            for i4 = 1:1:numel(ph1v)
                for i5 = 1:1:numel(ph2v)
                    for i6 = 1:1:numel(ph3v)
                       [x1p,x2p,x3p,x4p,x5p,x6p] = getpoints(th1v(i1),th2v(i2),th3v(i3),ph1v(i4),ph2v(i5),ph3v(i6),bar_l,edge_l);
                       val = check(x1p,x2p,x3p,x4p,x5p,x6p);
                       N = [x1p,x2p,x3p,x4p,x5p,x6p];
                       eqb = equilibriumtest(N);
                        if (val==1 && eqb==1) 
                        % save the points in the workspace variable
                        c = [c double((1/3).*(x4p+x5p+x6p))];
                        Nf(:,:,end+1) = N;
                        end
                    end
                end
            end
        end
    end
end
save('April6th_res_10.mat');
    
% line([x1p(1) x4p(1)],[x1p(2) x4p(2)],[x1p(3) x4p(3)],'Color','black');hold on;
% line([x2p(1) x5p(1)],[x2p(2) x5p(2)],[x2p(3) x5p(3)],'Color','black');hold on;
% line([x3p(1) x6p(1)],[x3p(2) x6p(2)],[x3p(3) x6p(3)],'Color','black');hold on;
%    

%plot the workspace
% scatter3(c(1,:),c(2,:),c(3,:),'filled')

fprintf('\n\nProgram Completed successfully.\n ')
toc

%% function declarations
function [x1,x2,x3,x4,x5,x6] = getpoints(th1,th2,th3,ph1,ph2,ph3,b,s)
d2r=pi/180;
x1 = [-s/2;-s/(2*sqrt(3));0];
x2 = [0;s/sqrt(3);0];
x3 = [s/2;-s/(2*sqrt(3));0];
x4 = [b*cos(ph1)*cos(th1);        b*cos(ph1)*sin(th1);         b*sin(ph1)]+x1;
x5 = [b*cos(ph2)*cos(240*d2r+th2);b*cos(ph2)*sin(240*d2r+th2); b*sin(ph2)]+x2;
x6 = [b*cos(ph3)*cos(120*d2r+th3);b*cos(ph3)*sin(120*d2r+th3); b*sin(ph3)]+x3;
end


function v = unit(v)
    v = v./norm(v);
end

function sp = check(x1,x2,x3,x4,x5,x6)
v1= x4-x1; v2=x5-x2; v3=x6-x3;
n1 = cross(v1,v2);
n2 = cross(v2,v3);
n3 = cross(v3,v1);

a1 = pinv([-(x4-x1) x5-x2 -cross(x4-x1,x5-x2)])*(x1-x2);
p1=x1+a1(1)*(x4-x1);
p2=x2+a1(2)*(x5-x2);
chkp1 = unit(p1-p2)./unit(n1);


a2 = pinv([-(x5-x2) x6-x3 -cross(x5-x2,x6-x3)])*(x2-x3);
p1=x2+a2(1)*(x5-x2);
p2=x3+a2(2)*(x6-x3);
chkp2 = unit(p1-p2)./unit(n2);

a3 = pinv([-(x6-x3) x4-x1 -cross(x6-x3,x4-x1)])*(x3-x1);
p1 = x3+a3(1)*(x6-x3);
p2 = x1+a3(2)*(x4-x1);
chkp3 = unit(p1-p2)./unit(n3);
Mat = double([chkp1,chkp2,chkp3]);
if vpa(sum(sum(Mat))) == 9.0
    sp = 1;
else
    sp= -1;
end
end

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