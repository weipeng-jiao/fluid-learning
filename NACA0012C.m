%NACA0024 C型网格生成


clear all
clc
%物理域绘制
% 画出左侧和上侧外界区域
x=-30:0.8:50;
for i=1:length(x)
    if i<=31
       y(i)=(900-x(i)^2)^0.5;
    elseif i<=length(x)
        y(i)=30;
    end
end
x0=x;y0=y;

x22=0:0.0125:1;  %画翼型(翼型单侧一共80个点)
y22=0.1781*x22.^0.5-0.0756*x22-0.2122*x22.^2+0.1705*x22.^3-0.0609*x22.^4;

%翼型上的数据储存(81个)
x_wall=[x22];
y_wall=[y22];
y_wall(:,end)=0;
l_wall=length(y_wall);
%画翼面以后的网格
x_low=3.33:2.33:50;
y_low=zeros(size(x_low));
l_low=length(y_low);

%左侧上侧边界数据储存（暂时只存到上边界结束）
x_left=[x0(1:l_wall)];
y_left=[y0(1:l_wall)];
l_left=length(y_left);
%上侧的后半部分
x_up=[x0(l_wall+1:l_wall+l_low-1)];
y_up=[y0(l_wall+1:l_wall+l_low-1)];
l_up=length(y_up);

n=9;%插值层数

%储存位置的大数组（暂时只存到上边界结束）
X=zeros(n+1,l_wall+l_low);
Y=zeros(n+1,l_wall+l_low);

%给X和Y的第一层（也就是第一行）赋值（分成n层，实际有n+1层数据）
X(1,1:l_wall)=x_wall;Y(1,1:l_wall)=y_wall;
X(1,l_wall+1:l_wall+l_low)=x_low;Y(1,l_wall+1:l_wall+l_low)=y_low;
X(n+1,1:l_wall)=x_left;Y(n+1,1:l_wall)=y_left;
X(n+1,l_wall+1:l_wall+l_up)=x_up;Y(n+1,l_wall+1:l_wall+l_up)=y_up;

    for j=1:l_wall+l_up
        a1=Y(n+1,j);
        a2=X(n+1,j);
        a=(Y(n+1,j)-Y(1,j))/(X(n+1,j)-X(1,j));
        b=Y(1,j);
        h=(X(n+1,j)-X(1,j))/n;
        for i=1:n   
            xi=X(1,j)+i*h;yi=a*(xi-X(1,j))+b;
            X(i+1,j)=xi;Y(i+1,j)=yi;
        end
    end
X(:,l_wall+l_up)=50*ones(n+1,1);
X(:,end)=[];Y(:,end)=[];
%物理域所有位置的坐标矩阵
% figure
X0=[fliplr(X) X];Y0=[fliplr(Y) -Y];
Z0=ones(size(X0));
figure
mesh(X0,Y0,Z0);
view(2)

[M,N]=size(X0);

%初始化
X2=X0;X2(:,N/2+1)=[];Y2=Y0;Y2(:,N/2+1)=[];X3=X2;Y3=Y2;

[M,N]=size(X2);

%规定误差
e0=0.001;

%进行迭代
for k=1:3000
    for i=2:M-1;
        for j=2:N-1
            a=((X2(i,j+1)-X2(i,j-1))^2+(Y2(i,j+1)-Y2(i,j-1))^2)/4;
            b=-(((X2(i+1,j)-X2(i-1,j))*(X2(i,j+1)-X2(i,j-1)))/4+((Y2(i+1,j)-Y2(i-1,j))*(Y2(i,j+1)-Y2(i,j-1)))/4);
            g=((X2(i+1,j)-X2(i-1,j))^2+(Y2(i+1,j)-Y2(i-1,j))^2)/4;
            X2(i,j)=(a*(X2(i+1,j)+X2(i-1,j))+0.5*b*(X2(i-1,j+1)+X2(i+1,j-1)-X2(i+1,j+1)-X2(i-1,j-1))+g*(X2(i,j+1)+X2(i,j-1)))...
               /(2*a+2*g) ;
           Y2(i,j)=(a*(Y2(i+1,j)+Y2(i-1,j))+0.5*b*(Y2(i-1,j+1)+Y2(i+1,j-1)-Y2(i+1,j+1)-Y2(i-1,j-1))+g*(Y2(i,j+1)+Y2(i,j-1)))...
               /(2*a+2*g) ;          
        end
    end
    %边界迭代
    %上段
    for i=1
        for j=2:n
            m=(N-1)-(j-1)*2;
            m=m+j;
            a=((X2(i,j+1)-X2(i,j-1))^2+(Y2(i,j+1)-Y2(i,j-1))^2)/4;
            b=-(((X2(i+1,j)-X2(i+1,m))*(X2(i,j+1)-X2(i,j-1)))/4+((Y2(i+1,j)-Y2(i+1,m))*(Y2(i,j+1)-Y2(i,j-1)))/4);
            g=((X2(i+1,j)-X2(i+1,m))^2+(Y2(i+1,j)-Y2(i+1,m))^2)/4;
            X2(i,j)=(a*(X2(i+1,j)+X2(i+1,m))+0.5*b*(X2(i+1,m+1)+X2(i+1,j-1)-X2(i+1,j+1)-X2(i+1,m-1))+g*(X2(i,j+1)+X2(i,j-1)))...
               /(2*a+2*g) ;
        end
    end
    %下段
    X2(1,N-l_up+1:N)=fliplr(X2(1,1:l_up));
    EX=abs(X3-X2);EY=abs(Y3-Y2);
    a=max(max(EX));b=max(max(EY));
    e=max(a,b);
    if e<e0
        break
    else 
        X3=X2;Y3=Y2;
    end
end
%     for m=2:20
%             X2(1,m)=(X2(1,m-1)+X2(1,m+1)+X2(2,m)+X2(1,203-m))/4;
%             Y2(1,m)=(Y2(1,m-1)+Y2(1,m+1)+Y2(2,m)+Y2(1,203-m))/4;
%         end
%         for n=183:202
%             X2(1,m)=(X2(1,m-1)+X2(1,m+1)+X2(2,m)+X2(1,203-m))/4;
%             Y2(1,m)=(Y2(1,m-1)+Y2(1,m+1)+Y2(2,m)+Y2(1,203-m))/4;
%         end          
%     figure
%     Z2=ones(size(X2));
%     mesh(X2,Y2,Z2);
%     view(2)
%     
%     EX=abs(X3-X2);EY=abs(Y3-Y2);
%     a=max(max(EX));b=max(max(EY));
%     e=max(a,b);
%     if e<e0
%         break
%     else 
%         X3=X2;Y3=Y2;
%     end
% end
% X2
figure
Z2=ones(size(X2));
mesh(X2,Y2,Z2);
view(2)
