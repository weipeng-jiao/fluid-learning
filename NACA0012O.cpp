/*-----------------------------------------------/
program:NACA0024-O-mesh-cp 
author:tony jiao
CopyRight:MIT
/-——————————————————————————————————————————————*/


#include<iostream>
#include<fstream>
#include<cmath>
using namespace std;
#define pi 3.141592
#define e 1e-7   //迭代误差
#define Vecioty 10 //速度
#define m 101 //离散点数
#define n 80  //层数
#define w 1.988   //松弛因子
#define q 1.94 //计算空间Zeta和Eta的值

double A(double x10,double x12,double y10,double y12)
{
	double sa;
    sa=(x12-x10)*(x12-x10)+(y12-y10)*(y12-y10);
	return sa;
}
double B(double x01,double x21,double x10,double x12,double y01,double y21,double y10,double y12)
{
	double sb;
    sb=(x21-x01)*(x12-x10)+(y21-y01)*(y12-y10);
	return sb;
}
double C(double x01,double x21,double y01,double y21)
{
	double sc;
    sc=(x21-x01)*(x21-x01)+(y21-y01)*(y21-y01);
	return sc;
}
void main()
{	
	int i,j,s=0,s1=0;
  double x[m][n],y[m][n],dx[m][n],df[m][n],dy[m][n],a[m][n],b[m][n],maxa,maxb,maxf,f[m][n],ef[m][n],u[m][n],v[m][n],beta[m][n],grma[m][n],Jcb[m][n],cp[m];
  
  //网格生成
	cout<<"网格生成"<<endl;
	
	
	//翼型坐标
  for(i=0;i<m;i++)
	{
		x[i][0]=0.499*(1+cos(2.0*i*pi/(m-1)));
		if(i<0.5*(m-1))
		{
				y[i][0]=-0.6*(-0.1015*pow(x[i][0],4)+0.284*pow(x[i][0],3)-0.3576*pow(x[i][0],2)-0.1221*x[i][0]+0.2969*sqrt(x[i][0]));
		}
		else 
		{  y[i][0]=0.6*(-0.1015*pow(x[i][0],4)+0.284*pow(x[i][0],3)-0.3576*pow(x[i][0],2)-0.1221*x[i][0]+0.2969*sqrt(x[i][0]));}
	
	}
	
	//远场坐标
    for(i=0;i<m;i++)
	{
		  x[i][n-1]=10*cos(2.0*i*pi/(m-1))+0.5;
		  f[i][n-1]=Vecioty*x[i][n-1];
		 if(i<0.5*(m-1)) 
		 {
      y[i][n-1]=-fabs(5*sin(2.0*i*pi/(m-1)));
		 }
		 else 
		 {
      y[i][n-1]=fabs(5*sin(2.0*i*pi/(m-1)));
     }
		    		 
	}
	
	//场内坐标初始化
   	for(i=0;i<m;i++)
	{	                  
		for(j=1;j<n-1;j++)
		{  
			x[i][j]=x[i][0]+j*(x[i][n-1]-x[i][0])/(n-1);
	        y[i][j]=y[i][0]+j*(y[i][n-1]-y[i][0])/(n-1);
			
		}
	}
	cout<<"正在生成网格数据中******"<<endl;
	cout<<"......"<<endl;
	 
	 
	//高斯-赛德尔迭代求解场内坐标
	do            
    {
	  
		for(i=1;i<m-1;i++)
		{
			for(j=1;j<n-1;j++)
			{  
			a[i][j]=x[i][j];
			b[i][j]=y[i][j];
		 	
		 	x[i][j]=0.5*(A(x[i][j-1],x[i][j+1],y[i][j-1],y[i][j+1])*(x[i-1][j]+x[i+1][j])-0.5*(B(x[i-1][j],x[i+1][j],x[i][j-1],x[i][j+1],y[i-1][j],y[i+1][j],y[i][j-1],y[i][j+1]))*(x[i+1][j+1]+x[i-1][j-1]-x[i-1][j+1]-x[i+1][j-1])+C(x[i-1][j],x[i+1][j],y[i-1][j],y[i+1][j])*(x[i][j-1]+x[i][j+1]))/(A(x[i][j-1],x[i][j+1],y[i][j-1],y[i][j+1])+C(x[i-1][j],x[i+1][j],y[i-1][j],y[i+1][j]));
		 	y[i][j]=0.5*(A(x[i][j-1],x[i][j+1],y[i][j-1],y[i][j+1])*(y[i-1][j]+y[i+1][j])-0.5*(B(x[i-1][j],x[i+1][j],x[i][j-1],x[i][j+1],y[i-1][j],y[i+1][j],y[i][j-1],y[i][j+1]))*(y[i+1][j+1]+y[i-1][j-1]-y[i-1][j+1]-y[i+1][j-1])+C(x[i-1][j],x[i+1][j],y[i-1][j],y[i+1][j])*(y[i][j-1]+y[i][j+1]))/(A(x[i][j-1],x[i][j+1],y[i][j-1],y[i][j+1])+C(x[i-1][j],x[i+1][j],y[i-1][j],y[i+1][j]));
		  x[0][j]=0.5*(A(x[0][j-1],x[0][j+1],y[0][j-1],y[0][j+1])*(x[m-2][j]+x[1][j])-0.5*(B(x[m-2][j],x[1][j],x[0][j-1],x[0][j+1],y[m-2][j],y[1][j],y[0][j-1],y[0][j+1]))*(x[1][j+1]+x[m-2][j-1]-x[m-2][j+1]-x[1][j-1])+C(x[m-2][j],x[1][j],y[m-2][j],y[1][j])*(x[0][j-1]+x[0][j+1]))/(A(x[0][j-1],x[0][j+1],y[0][j-1],y[0][j+1])+C(x[m-2][j],x[1][j],y[m-2][j],y[1][j]));
			y[0][j]=0.5*(A(x[0][j-1],x[0][j+1],y[0][j-1],y[0][j+1])*(y[m-2][j]+y[1][j])-0.5*(B(x[m-2][j],x[1][j],x[0][j-1],x[0][j+1],y[m-2][j],y[1][j],y[0][j-1],y[0][j+1]))*(y[1][j+1]+y[m-2][j-1]-y[m-2][j+1]-y[1][j-1])+C(x[m-2][j],x[1][j],y[m-2][j],y[1][j])*(y[0][j-1]+y[0][j+1]))/(A(x[0][j-1],x[0][j+1],y[0][j-1],y[0][j+1])+C(x[m-2][j],x[1][j],y[m-2][j],y[1][j]));
			
			x[m-1][j]=x[0][j];
			y[m-1][j]=y[0][j];
			x[i][j]=a[i][j]*(1-q)+x[i][j]*q;
			y[i][j]=b[i][j]*(1-q)+y[i][j]*q;
			dx[i][j]=fabs(a[i][j]-x[i][j]);
	    dy[i][j]=fabs(b[i][j]-y[i][j]);
			}		    
		}
			maxa=dx[1][1],maxb=dy[1][1];
		for(i=1;i<m-1;i++)
		{
			for(j=1;j<n-1;j++)
			{
				if(dx[i][j]>=maxa)
				{ 
					maxa=dx[i][j];
				}
			    if(dy[i][j]>=maxb)
				{
					maxb=dy[i][j];
				}
			}
		}
		s+=1;		
	}while(maxa>e||maxb>e); 
	cout<<"网格收敛迭代次数: "<<s<<endl;
	
	//求解速度势方程
	
	//初始化内场
	for(i=0;i<m;i++)
	{	                  
		for(j=0;j<n-1;j++)
		{  
			f[i][j]=x[i][j]+y[i][j];
	 	}
	}
	
	//计算雅可比矩阵和alpf、beta、gamma系数
	for(i=1;i<m-1;i++)  
	{
		beta[i][0]=(x[i+1][0]-x[i-1][0])*(4*x[i][1]-3*x[i][0]-x[i][2])+(y[i+1][0]-y[i-1][0])*(4*y[i][1]-3*y[i][0]-y[i][2]);
		beta[0][0]=(x[1][0]-x[m-2][0])*(4*x[0][1]-3*x[0][0]-x[0][2])+2*(y[1][0]-y[m-2][0])*(4*y[0][1]-3*y[0][0]-y[0][2]);
		beta[m-1][0]=beta[0][0];
	  Jcb[i][0]=(x[i+1][0]-x[i-1][0])*(4*y[i][1]-3*y[i][0]-y[i][2])-(y[i+1][0]-y[i-1][0])*(4*x[i][1]-3*x[i][0]-x[i][2]);
		Jcb[0][0]=(x[1][0]-x[m-2][0])*(4*y[0][1]-3*y[0][0]-y[0][2])-(y[1][0]-y[m-2][0])*(4*x[0][1]-3*x[0][0]-x[0][2]);
		Jcb[m-1][n-1]=Jcb[0][n-1];
		Jcb[m-1][0]=Jcb[0][0];
		grma[i][0]=(x[i+1][0]-x[i-1][0])*(x[i+1][0]-x[i-1][0])+(y[i+1][0]-y[i-1][0])*(y[i+1][0]-y[i-1][0]);
    grma[0][0]=(x[1][0]-x[m-2][0])*(x[1][0]-x[m-2][0])+(y[1][0]-y[m-2][0])*(y[1][0]-y[m-2][0]);
    grma[m-1][0]=grma[0][0];
		
	}
  cout<<"求解流场"<<endl;
   
  do                        
	{
	   for(j=1;j<n-1;j++) 
		{
			for(i=1;i<m-1;i++)
			{
				ef[i][j]=f[i][j];
				f[0][0]=(4*f[0][1]-f[0][2]-beta[0][0]/grma[0][0]*(f[1][0]-f[m-2][1]))/3;
        f[i][0]=(4*f[i][1]-f[i][2]-beta[i][0]/grma[i][0]*(f[i+1][0]-f[i-1][1]))/3;
			  f[m-1][0]=f[0][0]; 
        f[i][j]=0.5*(A(x[i][j-1],x[i][j+1],y[i][j-1],y[i][j+1])*(f[i-1][j]+f[i+1][j])-0.5*(B(x[i-1][j],x[i+1][j],x[i][j-1],x[i][j+1],y[i-1][j],y[i+1][j],y[i][j-1],y[i][j+1]))*(f[i+1][j+1]+f[i-1][j-1]-f[i-1][j+1]-f[i+1][j-1])+C(x[i-1][j],x[i+1][j],y[i-1][j],y[i+1][j])*(f[i][j-1]+f[i][j+1]))/(A(x[i][j-1],x[i][j+1],y[i][j-1],y[i][j+1])+C(x[i-1][j],x[i+1][j],y[i-1][j],y[i+1][j]));
        f[0][j]=0.5*(A(x[0][j-1],x[0][j+1],y[0][j-1],y[0][j+1])*(f[m-2][j]+f[1][j])-0.5*(B(x[m-2][j],x[1][j],x[0][j-1],x[0][j+1],y[m-2][j],y[1][j],y[0][j-1],y[0][j+1]))*(f[1][j+1]+f[m-2][j-1]-f[m-2][j+1]-f[1][j-1])+C(x[m-2][j],x[1][j],y[m-2][j],y[1][j])*(f[0][j-1]+f[0][j+1]))/(A(x[0][j-1],x[0][j+1],y[0][j-1],y[0][j+1])+C(x[m-2][j],x[1][j],y[m-2][j],y[1][j]));
			  f[m-1][j]=f[0][j];
				f[i][j]=ef[i][j]*(1-w)+f[i][j]*w;
     		df[i][j]=fabs(f[i][j]-ef[i][j]);
			}
	   }
    maxf=0;
	
	for(i=1;i<m-1;i++)
		{
			for(j=1;j<n-1;j++)
			{ 
								
				if(df[i][j]>=maxf)
				{ 
					maxf=df[i][j];
				}

			 }
		}
		s1+=1;		
	}while(maxf>=e);
	
	cout<<"计算流场收敛迭代次数: "<<s1<<endl;
   
   
   //求解流场速度分布
	for(i=1;i<m-1;i++)
		{
			for(j=1;j<n-1;j++)
			{ 
			  Jcb[i][j]=(x[i+1][j]-x[i-1][j])*(y[i][j+1]-y[i][j-1])-(y[i+1][j]-y[i-1][j])*(x[i][j+1]-x[i][j-1]);
        Jcb[0][j]=(x[1][j]-x[m-2][j])*(y[0][j+1]-y[0][j-1])-(y[1][j]-y[m-2][j])*(x[0][j+1]-x[0][j-1]);
			  u[i][j]=((f[i+1][j]-f[i-1][j])*(y[i][j+1]-y[i][j-1])-(f[i][j+1]-f[i][j-1])*(y[i+1][j]-y[i-1][j]))/Jcb[i][j];
			  v[i][j]=((f[i][j+1]-f[i][j-1])*(x[i+1][j]-x[i-1][j])-(f[i+1][j]-f[i-1][j])*(x[i][j+1]-x[i][j-1]))/Jcb[i][j];
			  u[i][0]=((f[i+1][0]-f[i-1][0])*(4*y[i][1]-3*y[i][0]-y[i][2])-(4*f[i][1]-3*f[i][0]-f[i][2])*(y[i+1][0]-y[i-1][0]))/Jcb[i][0];
			  v[i][0]=((4*f[i][1]-3*f[i][0]-f[i][2])*(x[i+1][0]-x[i-1][0])-(f[i+1][0]-f[i-1][0])*(4*x[i][1]-3*x[i][0]-x[i][2]))/Jcb[i][0];
        u[0][j]=((f[1][j]-f[m-2][j])*(y[0][j+1]-y[0][j-1])-(f[0][j+1]-f[0][j-1])*(y[1][j]-y[m-2][j]))/Jcb[0][j];
			  v[0][j]=((f[0][j+1]-f[0][j-1])*(x[1][j]-x[m-2][j])-(f[1][j]-f[m-2][j])*(x[0][j+1]-x[0][j-1]))/Jcb[0][j];
        u[m-1][j]=u[0][j];
			  v[m-1][j]=v[0][j];
			  u[0][0]=((f[1][0]-f[m-2][0])*(4*y[0][1]-3*y[0][0]-y[0][2])-(4*f[0][1]-3*f[0][0]-f[0][2])*(y[1][0]-y[m-2][0]))/Jcb[0][0];
			  v[0][0]=((4*f[0][1]-3*f[0][0]-f[0][2])*(x[1][0]-x[m-2][0])-(f[1][0]-f[m-2][0])*(4*x[0][1]-3*x[0][0]-x[0][2]))/Jcb[0][0];
        u[i][n-1]=10;
			  u[0][n-1]=10;
        v[0][n-1]=0;
        v[i][n-1]=0;
        u[m-1][n-1]=10;
			  v[m-1][n-1]=0;
        u[m-1][0]=u[0][0];
			  v[m-1][0]=v[0][0];
			}
		}
    ofstream outfile;
		ofstream outfile1;
    outfile.open("网格.dat");
    outfile1.open("-Cp与x关系图.dat");
    outfile<<"VARIABLES="<<"X"<<","<<"Y"<<","<<"u"<<","<<"v"<<endl;
    outfile<<"ZONE I="<<n<<","<<"J="<<m<<endl;
    
    for(i=0;i<m;i++)
		{
			for(j=0;j<n;j++)
			{
			  outfile<<x[i][j]<<"  "<<y[i][j]<<" "<<u[i][j]<<" "<<v[i][j]<<endl;
            }
		}
        for(i=(m-1)/2;i<m-1;i++)
		{
	           cp[i]=(u[i][0]*u[i][0]+v[i][0]*v[i][0])/Vecioty/Vecioty-1;
			   outfile1<<x[i][0]<<"  "<<cp[i]<<endl;
		}
       outfile.close();
       outfile1.close();
       cout<<"计算结束"<<endl;
	}
