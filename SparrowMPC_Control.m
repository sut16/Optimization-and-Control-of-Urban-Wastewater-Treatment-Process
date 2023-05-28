%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 滚动优化算法 by kzy 2022.03.16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = SparrowMPC_Control(u)
load("C:\Users\匡\Desktop\污水处理\BSM11\trained_bp_288.mat");
load("C:\Users\匡\Desktop\污水处理\BSM11\PS1_288",'PS1');
load("C:\Users\匡\Desktop\污水处理\BSM11\PS2_288",'PS2');

% load("C:\Users\匡\Desktop\污水处理\BSM11\rain_trained_bpnew_288.mat");
% load("C:\Users\匡\Desktop\污水处理\BSM11\rain_PS1new_288",'PS1');
% load("C:\Users\匡\Desktop\污水处理\BSM11\rain_PS2new_288",'PS2');


% load("C:\Users\匡\Desktop\污水处理\BSM11\storm_trained_bp_288.mat");
% load("C:\Users\匡\Desktop\污水处理\BSM11\storm_PS1_288",'PS1');
% load("C:\Users\匡\Desktop\污水处理\BSM11\storm_PS2_288",'PS2');

u=u';
ps1min_arr=PS1.xmin;
ps1rang_arr=PS1.xrange;
ps2min_arr=PS2.xmin;
ps2rang_arr=PS2.xrange;
u_length=length(u); 
so_error=u(u_length-3); 
sno_error=u(u_length-2);
So5Setting=u(u_length-1);
SNo2Setting=u(u_length);
u_InputWaterDataHanding=u(1:u_length-4); 
u_InputWaterDatalength=length(u_InputWaterDataHanding);
ps1min_arrin=ps1min_arr(1:u_InputWaterDatalength);
ps1rang_arrin=ps1rang_arr(1:u_InputWaterDatalength);
NorInputWaterData=data_handling(u_InputWaterDatalength,u_InputWaterDataHanding,ps1min_arrin,ps1rang_arrin);
FrontKla=u(u_length-5);
FrontQa=u(u_length-4);
ps1min_controlarr=ps1min_arr(u_length-4+1:end);
ps1min_controlrang=ps1rang_arr(u_length-4+1:end);
rate_so=1;
rate_sno=1;
SearchAgents_no=100;
Max_iteration=20;
np=1;
cycle_time=1;
for j=1:cycle_time
    [lb,ub,dim,fobj]=Get_Functions_details(np,FrontKla,FrontQa,sno_error); 
    X2=sobol_init(cycle_time*SearchAgents_no,dim,ub,lb);
    [SSA_curve2,value,~,pred_val]=SSA5(X2,SearchAgents_no,Max_iteration,lb,ub,dim,fobj,NorInputWaterData,net,ps1min_controlarr,ps1min_controlrang,ps2min_arr,ps2rang_arr,np,so_error,sno_error,So5Setting,SNo2Setting,rate_so,rate_sno,FrontKla,FrontQa);
    y=[SSA_curve2(1);SSA_curve2(2);value;pred_val(1);pred_val(2)];
end
end

%%% 麻雀优化算法滚动优化器

function [Best_pos2,Best_score2,curve4,pred_val]=SSA5(X0,pop,Max_iter,lb,ub,dim,fobj,NorInputWaterData,net,ps1min_controlarr,ps1min_controlrang,ps2min_arr,ps2rang_arr,np,so_error,sno_error,So5Setting,SNo2Setting, rate_so,rate_sno,FrontKla,FrontQa)

ST = 0.8;%预警值
SD = 0.1;%意识到有危险麻雀的比重
SDNumber = pop*SD;%意识到有危险麻雀数量
if(max(size(ub)) == 1)
   ub = ub*ones(1,dim);
   lb = lb*ones(1,dim);  
end
X = X0;
%计算初始适应度值
X3=sin_chaos_init(pop,dim,ub,lb);
X2=[X0;X3];   
[~,fitness_temp] =  fobj(np,net,NorInputWaterData,ps1min_controlarr,ps1min_controlrang,ps2min_arr,ps2rang_arr,X2,so_error,sno_error,So5Setting,SNo2Setting, rate_so,rate_sno,FrontKla,FrontQa,2*pop);
[fitness_temp, index]= sort(fitness_temp);%排序
GBestF = fitness_temp(1);%全局最优适应度值
for i = 1:pop
    X(i,:) = X2(index(i),:);
end
    [~,fitness] =  fobj(np,net,NorInputWaterData,ps1min_controlarr,ps1min_controlrang,ps2min_arr,ps2rang_arr,X,so_error,sno_error,So5Setting,SNo2Setting,rate_so,rate_sno,FrontKla,FrontQa,pop);
    GBestX = X(1,:);%全局最优位置
X_new = X;
%}
curve4=inf(1,Max_iter);
for i = 1: Max_iter
    
    BestF = fitness(1);
    WorstF = fitness(end);
    PD =exp(-1.5*i/Max_iter)-0.1;
    PDNumber = ceil(pop*PD); %自适应发现者 
    R2 = rand(1);
    c=exp(1-(i+Max_iter)/(Max_iter-i));
    for j = 1:PDNumber
        d=exp(-j/PDNumber);
         if(R2<ST)
            if(j==1)
                X_new(j,:) = X(j,:);
            else
                X_new(j,:) =  X(j,:)+randn().*pdist2(X(j,:),X(1,:))*ones(1,dim);
            end
        else
            if(j==1)
                X_new(j,:) = X(j,:);
            else
                step_1=(1-d)*pdist2(X(j,:),lb)*c;
               step_2=(1-d)*pdist2(X(j,:),X(1,:));
                if (step_1<step_2)
                    X_new(j,:) = X(j,:) + step_2*randn().*((X(j,:)-X(1,:))/(norm((X(j,:)-X(1,:)))));
                else
                     X_new(j,:) = X(j,:) + step_1*rand().*get_point(dim);
                end
            end
        end
        %}
    end
    to_one_fit=zeros(1,PDNumber);
    fit_value=zeros(1,PDNumber);
    fitness_new=zeros(1,pop);
    PDNumber1=ceil(0.5*c*PDNumber);
    for j = 1:PDNumber1
        to_one_fit(j)=(fitness(PDNumber1)-fitness(j))/(fitness(PDNumber1)-BestF);
        if(j==1)
            fit_value(j)=to_one_fit(j);
        else
            fit_value(j)=fit_value(j-1)+to_one_fit(j);
        end
    end
    for j = 1:PDNumber1
        fit_value(j)=fit_value(j)/fit_value(PDNumber1);
    end
    for j = PDNumber+1:pop
        if(j>(pop - PDNumber)/2 + PDNumber)
            X_new(j,:)= randn().*exp((X(end,:) - X(j,:))/j^2);      
        else 
            rand_arr=zeros(1,dim);
            for index=1:dim
                rand_arr(index)=2*rand()-1;
            end
            if(PDNumber1==1)
                ff=abs(X(j,:) - X(1,:));
                X_new(j,:)=X(1,:) + rand_arr.*ff;
            else
                for iter_t = 1:PDNumber1
                    if(c<fit_value(iter_t))
                        ff=abs(X(j,:) - X(iter_t,:));
                        X_new(j,:)=X(iter_t,:) +rand_arr.*ff;
                    break;
                    end
                end
            end
        end
    end
    Temp = randperm(pop);
    SDchooseIndex = Temp(1:SDNumber); 
    for j = 1:SDNumber
        if(fitness(SDchooseIndex(j))>BestF)
            X_new(SDchooseIndex(j),:) = X(1,:) + randn().*abs(X(SDchooseIndex(j),:) - X(1,:));
        end
    end
   store=zeros(1,dim);

%    
%    for j = 1
%         for index2=1:dim
%             if(X_new(j,index2) >(ub(dim)+lb(dim))/2)
%                 store(index2)=abs(X_new(j,index2))+abs(ub(dim)) ;
%             else
%                 store(index2)=abs(X_new(j,index2))+abs(lb(dim)) ;
%             end
%         end
%         a=norm(store)/(norm(X_new(j,:))+10^-90 );
%         temp=X_new(j,:);
%         val=(3*(power(dim,0.5)+1)/2);
%         X_new(j,:)=X_new(j,:)+a/val*randn()*X_new(j,:);%a/(3*(power(dim,0.5)+1)/2)*randn()*    
%         for q = 1: dim
%             if(X_new(j,q)>ub(q))
%                 X_new(j,q) =ub(q);%rand()*(ub(a)-lb(a))+lb(a);
%             end
%             if(X_new(j,q)<lb(q))
%                 X_new(j,q) =lb(q);%rand()*(ub(a)-lb(a))+lb(a);
%             end
%         end
%         X_in=[X_new(j,:);temp];
%         [~,value_X_in]=fobj(np,net,NorInputWaterData,ps1min_controlarr,ps1min_controlrang,ps2min_arr,ps2rang_arr,X_in,so_error,sno_error,So5Setting,SNo2Setting,rate_so,rate_sno,FrontKla,FrontQa,2);
%         %[~,value_2]=fobj(np,net,NorInputWaterData,ps1min_controlarr,ps1min_controlrang,ps2min_arr,ps2rang_arr,temp,so_error,sno_error,So5Setting,SNo2Setting,1);
%         value_1=value_X_in(1);value_2=value_X_in(2);
%         if(value_2<value_1)
%             X_new(j,:)=temp;
%         end
%    end   



%        
%    for index=1:dim
%       temp=X_new(1,:);
%       if(i<=Max_iter/4)
%           length2=( ub(index)+lb(index) )/2 ;
%           length1=2*abs(X_new(1,index)-(ub(index)+lb(index) )/2 ) ;
%       else
%           length1= abs(  X_new(2,index) - X_new(1,index)   );
%           length2=  (X_new(2,index)   + X_new(1,index) )/2;
%       end
%       if(X_new(1,index)<length2)
%           if(i>Max_iter/4)
%               X_new(1,index)=length2+rand()*(length1);
%           else
%               k=200;%exp(-0.5*(i-Max_iter/8));
%               X_new(1,index)= (ub(index)+lb(index) )/2+(ub(index)+lb(index) )/(2*k)-X_new(1,index)/k;
%           end
%       else
%           if(i>Max_iter/4)
%                   X_new(1,index)=length2-rand()*(length1);
%           else
%                   k=200;%exp(-0.5*(i-Max_iter/8));
%                   X_new(1,index)= (ub(index)+lb(index) )/2+(ub(index)+lb(index) )/(2*k)-X_new(1,index)/k;
%           end
%       end
%             if(X_new(1,index)>ub(index))
%                 X_new(1,index)=rand()*(ub(index)-lb(index))+lb(index);
%             end
%             if(X_new(1,index)<lb(index))
%                 X_new(1,index)=rand()*(ub(index)-lb(index))+lb(index);
%             end
% 
% 
% 
%             X_in_2=[temp;X_new(1,:)];
%             [~,value_in_2]=fobj(np,net,NorInputWaterData,ps1min_controlarr,ps1min_controlrang,ps2min_arr,ps2rang_arr,X_in_2,so_error,sno_error,So5Setting,SNo2Setting,rate_so,rate_sno,FrontKla,FrontQa,2);
%             t1=value_in_2(1);t2=value_in_2(2);
%             if(t1<t2)
%                 X_new(1,:)=temp;
%             end
%    end  
%        
   %边界控制
    for j = 1:pop
       for a = 1: dim
           if(X_new(j,a)>ub(a))
               X_new(j,a) =rand()*(ub(a)-lb(a))+lb(a);
           end
           if(X_new(j,a)<lb(a))
               X_new(j,a) =rand()*(ub(a)-lb(a))+lb(a);
           end
       end
    end 
   %更新位置
        [~,fitness_new] = fobj(np,net,NorInputWaterData,ps1min_controlarr,ps1min_controlrang,ps2min_arr,ps2rang_arr,X_new,so_error,sno_error,So5Setting,SNo2Setting,rate_so,rate_sno,FrontKla,FrontQa,pop);
   for j = 1:pop
       if(fitness_new(j) < GBestF)
            GBestF = fitness_new(j);
            GBestX = X_new(j,:);   
        end
   end
   X = X_new;
   fitness = fitness_new;
    %排序更新
   [fitness, index]= sort(fitness);%排序
   BestF = fitness(1);
   WorstF = fitness(end);
   for j = 1:pop
      X(j,:) = X_new(index(j),:);
   end
   curve4(i) = GBestF;
   
end
Best_pos2 =GBestX;
Best_score2 = curve4(end);
[pred_val,~] = fobj(np,net,NorInputWaterData,ps1min_controlarr,ps1min_controlrang,ps2min_arr,ps2rang_arr,GBestX,so_error,sno_error,So5Setting,SNo2Setting,rate_so,rate_sno,FrontKla,FrontQa,1);
end

%产生sobol序列
function  Positions=sobol_init(length_arr,dim,ub,lb)
Positions=zeros(length_arr,dim);
 p = sobolset(dim);
Boundary_no= size(ub,2); 
if Boundary_no==1
    for index=1:length_arr
        Positions(index,:)=(p(index+1,:).*(ub-lb))+lb;
    end
end
 if Boundary_no>1
     for index=1:length_arr
         for i=1: dim
            ub_i=ub(i);
            lb_i=lb(i);
            Positions(index,i)=p(index,i).*(ub_i-lb_i)+lb_i;
         end
     end
 end
end

%%%%sin混沌映射
function  Positions=sin_chaos_init(pop,dim,ub,lb)
Positions=inf(pop,dim);
init_value=2*rand()-1;
while(init_value==0)
    init_value=2*rand()-1;
end
Boundary_no= size(ub,2); 
if Boundary_no==1
    for j=1:dim
        for index=1:pop
            a=sin(2/init_value);
            b=ub-lb;
        Positions(index,j)=(a+1)/2*b+lb;
        init_value=sin(2/init_value);
        end
    end
end
 if Boundary_no>1
         for j=1:dim
              for index=1:pop
                ub_i=ub(j);
                lb_i=lb(j);
                a=sin(2/init_value);
                b=ub_i-lb_i;
                Positions(index,j)=(a+1)/2*b+lb_i;
                init_value=sin(2/init_value);
              end
        end
 end
end

%%%n维超球面均匀分布单位向量
function x=get_point(dim)
x=zeros(1,dim);
for index=1:dim
    x(index)=randn();
end
m=norm(x);
x=x/m;
end

function [lb,ub,dim,fobj] = Get_Functions_details(np,FrontKla,FrontQa,sno_error)
        fobj = @F1;
        % lb=zeros(1,np*2);
        if sno_error==0
            deltaQa=50000;
        else
            deltaQa=10000;
        end
          lb(1)=0;
          lb(2)=FrontQa-deltaQa;
          if(lb(2)<30)
                lb(2)=30;
         end
         ub=zeros(1,np*2);
        for index=1:np
            ub(2*(index-1)+1)=240;
            ub(2*index)=FrontQa+deltaQa;
            if(ub(2*index)>92230)
                ub(2*index)=92230;
            end
        end
        dim=2*np; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 函数F1   输入：各个变量   输出：预测值和此时对应的f
%%% x是控制变量 u里面的最后2个是设定 倒数第三第四是So5和Sno2的误差 
%%% 倒数第五第六是上一步的溶解氧系数和内回流量
%%% np为预测步长，net为训练好的神经网络,u为优化算法输入
%%% ps1min_arr、ps1rang_arr为了变量归一化  ps2min_arr、ps2rang_arr为了反归一化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%,ps1min_controlarr,ps1min_controlrang

function [ret_va,f] = F1(np,net,NorInputWaterData,ps1min_controlarr,ps1min_controlrang,ps2min_arr,ps2rang_arr,x,so_error,sno_error,So5Setting,SNo2Setting,rate_so,rate_sno,FrontKla,FrontQa,NumData)
x_nor=x;
% for index=1:NumData
%     x_nor(index,:)=ControlDataHanding(np,2,x(index,:),ps1min_controlarr,ps1min_controlrang);
% end
x_nor(:,1)=(x_nor(:,1)-repmat(ps1min_controlarr(1),NumData,1))/(ps1min_controlrang(1));
x_nor(:,2)=(x_nor(:,2)-repmat(ps1min_controlarr(2),NumData,1))/(ps1min_controlrang(2));
input_val=[repmat(NorInputWaterData,NumData,1) x_nor]; %神经网络输入数据
%control_store(1,:)=x_arr(1,:); %%%储存第一个控制变量
%%%得到预测值   %%% x是控制变量 u里面的最后2个是设定 倒数第三第四是修正
store=get_pred_value(input_val,net,ps2min_arr,ps2rang_arr,so_error,sno_error,rate_so,rate_sno,NumData);
f=zeros(NumData,1);
%%%求取f
rateKla=0.0;
rateQa=0.0;
    for index=1:NumData
        f(index)=f(index)+(store(index,1)/So5Setting-1)^2+(store(index,2)/SNo2Setting-1)^2+rateKla*((FrontKla-x(1))/240)^2+rateQa*((FrontQa-x(2))/92230)^2;
    end
   %f=(store(:,1)-So5Setting).^2+4*(store(:,2)-(SNo2Setting)).^2;
    ret_va=store;%%返回后一步预测值
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 函数get_pred_value   输入：各个变量   输出：预测值
%%% net为训练好的神经网络,u为优化算法输入
%%% ps1min_arr、ps1rang_arr为了变量归一化  ps2min_arr、ps2rang_arr为了反归一化
%%% so_error为上一步的溶解氧预测误差 sno_error为上一步的硝态氮的预测误差
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pred_value=get_pred_value(test_xarr,net,ps2min_arr,ps2rang_arr,so_error,sno_error,rate_so,rate_sno,NumData)
O2 = net(test_xarr'); %%%得到预测值
O2=O2';
OUT2=zeros(NumData,2);
OUT2(:,1)=O2(:,1)*ps2rang_arr(1)+repmat(ps2min_arr(1)+rate_so*so_error,NumData,1);
OUT2(:,2)=O2(:,2)*ps2rang_arr(2)+repmat(ps2min_arr(2)+rate_sno*sno_error,NumData,1);
pred_value=OUT2;
end

%%%入水数据归一化处理
function output_data=data_handling(handling_length,u,ps1min_arr,ps1rang_arr)
min_arr=ps1min_arr;
rang_arr=ps1rang_arr;
for index=1:handling_length
    u(index)=(u(index)-min_arr(index))/rang_arr(index);
end
output_data=u;
end

                 



