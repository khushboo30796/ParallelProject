#include <stdlib.h>  
#include <stdio.h>
#include <math.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <omp.h>
using namespace std;
using namespace std::chrono;
struct vector
{	
	float x,y,z;
	vector operator+(vector r)
	{
		return vector(x+r.x,y+r.y,z+r.z);
	}
	vector operator*(float r)
	{
		return vector(x*r,y*r,z*r);
	}
	float operator%(vector r)
	{
		return x*r.x+y*r.y+z*r.z;
	}
	vector()
	{}
	vector operator^(vector r)
	{
		return vector(y*r.z-z*r.y,z*r.x-x*r.z,x*r.y-y*r.x);
	}
	vector(float a,float b,float c)
	{
		x=a;y=b;z=c;
	}
	vector operator!()
	{
		return*this*(1/sqrt(*this%*this));
	}
	void init(float a,float b,float c)
	{
		x=a;y=b;z=c;
	}
};
int * G;
void allocG(int n)
{
	G=(int *)calloc(10,sizeof(int));
	bool arr[10][20];
	for(int i=0;i<10;i++)
		for(int j=0;j<20;j++)
		{
			arr[i][j]=false;
			G[i]=0;
		}
	int count=0;
	while(count<n)
	{
		int r=rand()%10;
		int c=rand()%20;
		if(arr[r][c]==false)
		{
			arr[r][c]=true;
			count++;
		}
	}
	for(int i=0;i<10;i++)
		for(int j=0;j<20;j++)
			if(arr[i][19-j]==true)
				G[i]+=(int)pow(2,j);
}
float Random()
{
	return(float)rand()/RAND_MAX;
}
   //The intersection test for line [o,v].
   // Return 2 if a hit was found (and also return distance t and bouncing ray n).
   // Return 0 if no hit was found but ray goes upward
   // Return 1 if no hit was found but ray goes downward
int Trace(vector o,vector d,float&t,vector&n)
{
	t=1e9;
	int m=0;
	float p=-o.z/d.z;
	if(.01<p)
	{
		t=p,n=vector(0,0,1),m=1;
	}
	for(int k=20;k--;)
		for(int j=10;j--;)
			if(G[j]&1<<k)
			{
				vector p=o+vector(-k,0,-j-4);
				float b=p%d;
				float c=p%p-1;
				float q=b*b-c;
				if(q>0)
				{
					float s=-b-sqrt(q);
					if(s<t&&s>.01)
					{
						t=s,n=!(p+d*t),m=2;
					}
				}
			}
		return m;
}
vector Sample(vector o,vector d)
{
	float t;
	vector n;
	int m=Trace(o,d,t,n);
	if(!m)
		return vector(.7,.6,1)*pow(1-d.z,4);
	vector h=o+d*t;
	vector l=!(vector(9+Random(),9+Random(),16)+h*-1);
	vector r=d+n*(n%d*-2);
	float b=l%n;
	if(b<0||Trace(h,l,t,n))
		b=0;
	float p=pow(l%r*(b>0),99);
	if(m&1)
	{
		h=h*.2;
		return((int)(ceil(h.x)+ceil(h.y))&1?vector(3,1,1):vector(3,3,3))*(b*.2+.1);
	}
	return vector(p,p,p)+Sample(h,r)*.5;
}
vector * rayTracerParallel()
{
	vector * result=(vector *)calloc(512*512,sizeof(vector));	
	vector g=!vector(-12,-32,0);
	vector a=!(vector(0,0,1)^g)*.002;
	vector b=!(g^a)*.002;
	vector c=(a+b)*-256+g;
	#pragma omp parallel for
	for(int y=512;y>0;y--)
		for(int x=512;x>0;x--)
		{
			vector p(13,13,13);
			for(int r=64;r--;)
			{
				vector t=a*(Random()-.5)*99+b*(Random()-.5)*99;
				p=Sample(vector(17,16,8)+t,!(t*-1+(a*(Random()+x)+b*(y+Random())+c)*16))*3.5+p;
			}
			result[(x-1)*512+y-1]=p;
		}
	return result;
}
vector * rayTracer()
{
	vector * result=(vector *)calloc(512*512,sizeof(vector));
	vector g=!vector(-6,-16,0);
	vector a=!(vector(0,0,1)^g)*.002;
	vector b=!(g^a)*.002;
	vector c=(a+b)*-256+g;
	for(int y=512;y>0;y--)
		for(int x=512;x>0;x--)
		{
			vector p(13,13,13);
			for(int r=64;r--;)
			{
				vector t=a*(Random()-.5)*99+b*(Random()-.5)*99;
				p=Sample(vector(17,16,8)+t,!(t*-1+(a*(Random()+x)+b*(y+Random())+c)*16))*3.5+p;
			}
			result[(x-1)*512+y-1]=p;
		}
	return result;
}
int main()
{
	int n=100;
	allocG(n);
	printf("P6 512 512 255 ");
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
		vector * result1=rayTracer();
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>( t2 - t1 ).count();
	high_resolution_clock::time_point t1p = high_resolution_clock::now();
		vector * result2=rayTracerParallel();
	high_resolution_clock::time_point t2p = high_resolution_clock::now();
	auto duration2 = duration_cast<microseconds>( t2p - t1p ).count();
	for(int y=512;y>0;y--)
			for(int x=512;x>0;x--)
				printf("%c%c%c",(int)(result2[(x-1)*512+y-1].x),(int)(result2[(x-1)*512+y-1].y),(int)(result2[(x-1)*512+y-1].z));

	fstream f;
	f.open("speedupcopy.txt",ios::out|ios::app);
	f<<n<<" "<<duration<<" "<<duration2<<" "<<(double)duration/duration2<<'\n';
	f.close();
}