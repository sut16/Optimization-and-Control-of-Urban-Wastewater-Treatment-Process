
function x_out=cmopMOEADARP(u)
load("C:\Users\匡\Desktop\污水处理\BSM11\trainedECEQ_bp_1.mat");
load("C:\Users\匡\Desktop\污水处理\BSM11\ECEQPS1_1",'PS1');
load("C:\Users\匡\Desktop\污水处理\BSM11\ECEQPS2_1",'PS2');
global u_nor ps1min_arr ps1rang_arr ps2min_arr ps2rang_arr M1 
ps1min_arr=PS1.xmin;
ps1rang_arr=PS1.xrange;
ps2min_arr=PS2.xmin;
ps2rang_arr=PS2.xrange;
M1=net;
u=u';
%u=[20496	,30.7798500000000	,61.0327800000000,	2.33262550577549,	15.0105992860157];
%u=[21874.9999996160	,26.4558899910947	,51.2655199708211];
 %u=[18612.9999643516	,28.2786999634876	,50.1706299403334];
% 27780.0000566720	46.6146800842669	73.7548301780374	0.263902449321137	14.1418799438871
% 29551	49.2480200000000	79.3185000000000	0.239395864470295	14.0063679528508
% 31861.0000198400	49.9994099519110	82.7012900510707	0.232769888287721	13.8628509600825
% 32171.0000002880	49.2480198877927	83.4992699090304	0.257388268130415	13.7354539878173
% 32180	45.7415400000000	80.6564700000000	0.339110372860929	13.6405871100740
% 31868.9997553920	40.8575098833395	75.8099998041767	0.521137725595573	13.6185843503718
% 28046.9999567360	39.0346899300931	72.7502598997472	0.829834559449627	13.6842566093755
% 26695	36.8501000000000	69.6173600000000	1.26207569446677	13.8399807066980
% 26689.9999284480	34.2480698664199	66.1858697991066	1.82736197758988	14.1139046578533
% 25571.9999975680	32.1608799558071	63.0469099355479	2.53174921279223	14.5178699648484
% 25496	30.7798500000000	61.0327800000000	3.33262550577549	15.0105992860157
uin=u;ps1min_arru= ps1min_arr(3:end);ps1rang_arru=ps1rang_arr(3:end);
length_u=length(uin);
u_nor=data_handling(length_u,uin,ps1min_arru,ps1rang_arru);
[Dec,obj1]=platemo('objFcn',{@f1,@f2},'conFcn',@g1,'lower',[0.4,0.3],'upper',[3,2],'algorithm',@MOEADNCS,'N',50,'maxFE',2500);
[obj1,Dec]=sortpopandDEC(obj1,Dec);
% plot(obj1(:,1),obj1(:,2),'b.');
% hold on
[kneepoint,num]=findkneepoint(obj1);
% plot(kneepoint(1),kneepoint(2),'r.');
% [val,num]=min(obj1);
x_out=Dec(num,:);

end



function y=f1(x)
global u_nor ps1min_arr ps1rang_arr M1 ps2min_arr ps2rang_arr tran tranSNH predtotN
x(1)=(x(1)-ps1min_arr(1))/ps1rang_arr(1);
x(2)=(x(2)-ps1min_arr(2))/ps1rang_arr(2);
x=[x,u_nor];
%O2 = M1.GetOutput(x);
O2 = M1(x');
O2=O2';
%predSNH=O2(:,3)*ps2rang_arr(3)+ps2min_arr(3);
predtotN=O2(:,4)*ps2rang_arr(4)+ps2min_arr(4);
 %pealtySNHfac=400000*max(predSNH-4,0);pealtytotNfac=1000*max(predtotN-18,0);
tranSNH=O2(:,3)*ps2rang_arr(3)+ps2min_arr(3);
 y=O2(:,1)*ps2rang_arr(1)+ps2min_arr(1);
 tran=O2(:,2)*ps2rang_arr(2)+ps2min_arr(2);


end

function y=f2(x)
%global u_nor ps1min_arr ps1rang_arr M1 ps2min_arr ps2rang_arr 
global tran
% x(1)=(x(1)-ps1min_arr(1))/ps1rang_arr(1);
% x(2)=(x(2)-ps1min_arr(2))/ps1rang_arr(2);
% x=[x,u_nor];
% O2 = M1(x');
% O2=O2';
% y=O2(:,2)*ps2rang_arr(2)+ps2min_arr(2);
y=tran;
end

function output_data=data_handling(handling_length,u,ps1min_arr,ps1rang_arr)
min_arr=ps1min_arr;
rang_arr=ps1rang_arr;
for index=1:handling_length
    u(index)=(u(index)-min_arr(index))/rang_arr(index);
end
output_data=u;
end
% 
function y=g1(x)
global tranSNH
y=tranSNH-40;
end
% function y=g2(x)
% global predtotN
% y=predtotN-inf;
% end

function [kneepoint,num]=findkneepoint(popobj)
%%%%找到离极端点连成直线最远距离的点即可
[N,~]=size(popobj);
f1minpoint=popobj(1,:);f2minpoint=popobj(end,:);
popdis_store=zeros(1,N);
for index=1:N
     popdis_store(index)= (det([f2minpoint-f1minpoint;popobj(index,:)-f1minpoint]))/norm(f2minpoint-f1minpoint);
end
[~,num]=min(popdis_store);
kneepoint=popobj(num,:);
end

function [Functionval,dec1]=sortpopandDEC(popobj,dec)
[temp,index]=sort(popobj(:,1));
popobj(:,2)=popobj(index,2);
Functionval=[temp, popobj(:,2)];
dec1=dec(index,:);

end
