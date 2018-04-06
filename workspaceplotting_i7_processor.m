% calculate workspace of Prism
% (0,0,0) origin point is at the center
clear all;clc;close all;
d2r = pi/180;
% syms b s th1 th2 th3 ph1 ph2 ph3 real 
% syms t1 t2 t3 positive

% x1 = [0;s/sqrt(3);0];
% x2 = [-s/2;-s/(2*sqrt(3));0];
% x3 = [s/2;-s/(2*sqrt(3));0];
% x4 = [b*cos(pi-ph1)*cos(240*d2r+th1);b*cos(pi-ph1)*sin(240*d2r+th1);b*sin(pi-ph1)]+x1;
% x5 = [b*cos(ph2)*cos(th2);b*cos(ph2)*sin(th2);b*sin(ph2)]+x2;
% x6 = [b*cos(ph3)*cos(120*d2r+th3);b*cos(ph3)*sin(120*d2r+th3);b*sin(ph3)]+x3;

fprintf('Program initialized');

tic
%% initialize values and compute the workspace

dth=30;
th1v = [-10 0 10 15 20 30];th2v=[-10 0 10 15 20 30];th3v=[-10 0 10 15 20 30 ];
ph1v = [15 30 45 60 75];ph2v=[15 30 45 60 75];ph3v=[15 30 45 60 75];
th1v=d2r.*th1v;th2v=d2r.*th2v;th3v=d2r.*th3v;
ph1v=d2r.*ph1v;ph2v=d2r.*ph2v;ph3v=d2r.*ph3v;

bar_l = 2;
edge_l = 2;
j=1;


% % false condition
% th1v = 20;th2v=20;th3v=20;
% ph1v = 10;ph2v=30;ph3v=30;

% % true condition
% th1v = 25;th2v=25;th3v=25;
% ph1v = 20;ph2v=20;ph3v=20;
% % 
% th1v=d2r.*th1v;th2v=d2r.*th2v;th3v=d2r.*th3v;
% ph1v=d2r.*ph1v;ph2v=d2r.*ph2v;ph3v=d2r.*ph3v;

% %%
% [x1p,x2p,x3p,x4p,x5p,x6p] = getpoints(th1v,th2v,th3v,ph1v,ph2v,ph3v,bar_l,edge_l);
% line([x1p(1) x4p(1)],[x1p(2) x4p(2)],[x1p(3) x4p(3)],'Color','black');hold on;
% line([x2p(1) x5p(1)],[x2p(2) x5p(2)],[x2p(3) x5p(3)],'Color','black');hold on;
% line([x3p(1) x6p(1)],[x3p(2) x6p(2)],[x3p(3) x6p(3)],'Color','black');hold on;
                
figure(1)
axis([-3 3 -3 3 0 3])
c=[];
% the great search starts
for i1=1:1:numel(th1v)
    for i2 = 1:1:numel(th2v)
        for i3 = 1:1:numel(th3v)
            for i4 = 1:1:numel(ph1v)
                for i5 = 1:1:numel(ph2v)
                    for i6 = 1:1:numel(ph3v)
                       [x1p,x2p,x3p,x4p,x5p,x6p] = getpoints(th1v(i1),th2v(i2),th3v(i3),ph1v(i4),ph2v(i5),ph3v(i6),bar_l,edge_l);
                       val = check(x1p,x2p,x3p,x4p,x5p,x6p);
                        if val==1
                        c = [c double((1/3).*(x4p+x5p+x6p))];
                        end
                    end
                end
            end
        end
    end
end
save('March10th.mat','c');
fprintf('Saved the values for d_theta = 5');
%     
%     line([x1p(1) x4p(1)],[x1p(2) x4p(2)],[x1p(3) x4p(3)],'Color','black');hold on;
%     line([x2p(1) x5p(1)],[x2p(2) x5p(2)],[x2p(3) x5p(3)],'Color','black');hold on;
%     line([x3p(1) x6p(1)],[x3p(2) x6p(2)],[x3p(3) x6p(3)],'Color','black');hold on;
%     end

%plot the workspace
scatter3(c(1,:),c(2,:),c(3,:),'filled')

fprintf('\n\nProgram Completed successfully.\n ')
toc

%% function declarations
function [x1,x2,x3,x4,x5,x6] = getpoints(th1,th2,th3,ph1,ph2,ph3,b,s)
d2r=pi/180;
x1 = [0;s/sqrt(3);0];
x2 = [-s/2;-s/(2*sqrt(3));0];
x3 = [s/2;-s/(2*sqrt(3));0];
x4 = [b*cos(ph1)*cos(240*d2r+th1);b*cos(ph1)*sin(240*d2r+th1); b*sin(ph1)]+x1;
x5 = [b*cos(ph2)*cos(th2);        b*cos(ph2)*sin(th2);         b*sin(ph2)]+x2;
x6 = [b*cos(ph3)*cos(120*d2r+th3);b*cos(ph3)*sin(120*d2r+th3); b*sin(ph3)]+x3;
end

% function [n1,n2,n3] = get_normals(th1,th2,th3,ph1,ph2,ph3,b,s)
% n1 = [                      - b^2*cos(ph2)*sin(ph1)*sin(th2) - b^2*cos(ph1)*sin(ph2)*sin(th1 + (4*pi)/3);
%                               b^2*cos(ph2)*cos(th2)*sin(ph1) + b^2*cos(ph1)*sin(ph2)*cos(th1 + (4*pi)/3);
%  b^2*cos(ph1)*cos(ph2)*cos(th2)*sin(th1 + (4*pi)/3) - b^2*cos(ph1)*cos(ph2)*cos(th1 + (4*pi)/3)*sin(th2)];
% 
% n2 = [                        b^2*cos(ph2)*sin(ph3)*sin(th2) - b^2*cos(ph3)*sin(ph2)*sin(th3 + (2*pi)/3);
%                               b^2*cos(ph3)*sin(ph2)*cos(th3 + (2*pi)/3) - b^2*cos(ph2)*cos(th2)*sin(ph3);
%  b^2*cos(ph2)*cos(ph3)*cos(th2)*sin(th3 + (2*pi)/3) - b^2*cos(ph2)*cos(ph3)*cos(th3 + (2*pi)/3)*sin(th2)];
% 
% n3 = [                         b^2*cos(ph1)*sin(ph3)*sin(th1 + (4*pi)/3) + b^2*cos(ph3)*sin(ph1)*sin(th3 + (2*pi)/3);
%                              - b^2*cos(ph1)*sin(ph3)*cos(th1 + (4*pi)/3) - b^2*cos(ph3)*sin(ph1)*cos(th3 + (2*pi)/3);
%  b^2*cos(ph1)*cos(ph3)*cos(th1 + (4*pi)/3)*sin(th3 + (2*pi)/3) - b^2*cos(ph1)*cos(ph3)*cos(th3 + (2*pi)/3)*sin(th1 + (4*pi)/3)];
% n1 = n1/norm(n1);
% n2 = n2/norm(n2);
% n3 = n3/norm(n3);
% end

function v = unit(v)
    v = v./norm(v);
end

% function val = check(x1,x2,x3,x4,x5,x6)
% syms k t1 t2 t3 real
% v1= x4-x1; v2=x5-x2; v3=x6-x3;
% n1 = cross(v1,v2);
% n2 = cross(v2,v3);
% n3 = cross(v3,v1);
% 
% L1 = x1 + t1*(x4-x1);
% L2 = x2 + t2*(x5-x2);
% L3 = x3 + t3*(x6-x3);
% 
% eqn1 = L2-L1- k*n1 == 0;
% [A,B] = equationsToMatrix(eqn1, [t1,t2,k]);
% sol=double(A^-1*B);
% p1 = double(subs(L1,t1,sol(1)));
% p2 = double(subs(L2,t2,sol(2)));
% chkp1 = unit(p1-p2)./unit(n1);
% 
% eqn2 = L3-L2- k*n2 == 0;
% [A,B] = equationsToMatrix(eqn2, [t2,t3,k]);
% sol=double(A^-1*B);
% p1 = double(subs(L2,t2,sol(1)));
% p2 = double(subs(L3,t3,sol(2)));
% chkp2 = unit(p1-p2)./unit(n2);
% 
% eqn3 = L1-L3- k*n3 == 0;
% [A,B] = equationsToMatrix(eqn3, [t3,t1,k]);
% sol=double(A^-1*B);
% p1 = double(subs(L3,t3,sol(1)));
% p2 = double(subs(L1,t1,sol(2)));
% chkp3 = unit(p1-p2)./unit(n3);
% Mat = double([chkp1,chkp2,chkp3]);
% if vpa(sum(sum(Mat))) == 9.0
%     val = 1;
% else
%     val = 0;
% end
% 
% end

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