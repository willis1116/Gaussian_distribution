#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

using namespace std;

class range
{
public:
    range(double x,double y,double z,double dx,double dy,double dz,double intheta,double inu)
    {
        theta=(intheta*M_PI)/180;
        u=inu;
        l=(x/dx)+1;
        w=(y/dy)+1;
        h=(z/dz)+1;
        t=l*w*h;
        X=new double[t];
        Y=new double[t];
        Z=new double[t];
        int i,j,k,n;
        for(k=0; k<h; k++) // build the grids, from the bottom to top
        {
            for(j=0; j<w; j++)
            {
                for(i=0; i<l; i++)
                {
                    n=i+l*j+(w*w)*k;
                    X[n]=i*dx;
                    Y[n]=j*dy;
                    Z[n]=k*dz;
                }
            }
        }
    }
    range (const range&i) //宣告複製建構值
    {
        theta=i.theta;
        u=i.u;
        l=i.l;
        w=i.w;
        h=i.h;
        t=i.t;
        ww=i.ww;
        ll=i.ll;
        X=new double[t];
        Y=new double[t];
        Z=new double[t];
        int o;
        for(o=0; o<t; o++)
        {
            X[o]=i.X[o];
            Y[o]=i.Y[o];
            Z[o]=i.Z[o];
        }
    }

    int getl()
    {
        ll=l;
        return ll;
    }
    int getw()
    {
        ww=w;
        return ww;
    }
    int gett()
    {
        return t;
    }
    double *getX()
    {
        return X; //回傳X位址
    }
    double *getY() //宣告取值函式getY
    {
        return Y;  //回傳Y位址
    }
    double *getZ() //宣告取值函式getZ
    {
        return Z; //回傳Z位址
    }
    double gettheta() //宣告取值函式gettheta
    {
        return theta;
    }
    double getu()
    {
        return u;
    }
private:
    int t,ww,ll;
    double l, w, h, theta, u;
    double *X,*Y,*Z;
};

class CHIM //chimney
{
public:
    CHIM()
    {
        cx=10.0;
        cy=10.0;
        ch=1.0;
        Q=2.0;
    }
    CHIM(double incx,double incy,double inch,double inQ) //data from the input file
    {
        cx=incx;
        cy=incy;
        ch=inch;
        Q=inQ;
    }
    CHIM(const CHIM&i)  //宣告複製建構值
    {
        cx=i.cx; //將原本已經存在的函式中的cx存入新的物件的cx中
        cy=i.cy; //將原本已經存在的函式中的cy存入新的物件的cy中
        ch=i.ch; //將原本已經存在的函式中的ch存入新的物件的ch中
        Q=i.Q; //將原本已經存在的函式中的Q存入新的物件的Q中
    }
    void operator=(const CHIM&i) //宣告等號重載運算值
    {
        cx=i.cx; //將原本已經存在的函式中的cx存入新的物件的cx中
        cy=i.cy; //將原本已經存在的函式中的cy存入新的物件的cy中
        ch=i.ch;
        Q=i.Q;
    }
    void set_value(double incx,double incy,double inch,double inQ)
    {
        cx=incx;
        cy=incy;
        ch=inch;
        Q=inQ;
    }

    double getcx()
    {
        return cx;
    }
    double getcy()
    {
        return cy;
    }
    double getch()
    {
        return ch;
    }
    double getQ()
    {
        return Q;
    }
protected:
    double cx,cy,ch,Q;
};

class chimney : public CHIM
{
public:
    chimney()
    {
    }
    chimney(int t)
    {
        grid_nums=t;
        C=new double [grid_nums]();
    }
    void init(double incx,double incy,double inch,double inQ,int t)
    {
        cx=incx;
        cy=incy;
        ch=inch;
        Q=inQ;
        grid_nums=t;
        C=new double [grid_nums];
    }
    chimney(const chimney&i)  //宣告複製建構值
    {
        cx=i.cx;
        cy=i.cy;
        ch=i.ch;
        Q=i.Q;
        grid_nums=i.grid_nums;
        C=new double [grid_nums];
        for(int o=0; o<grid_nums; o++)
            C[o]=i.C[o];
    }
    void operator=(const chimney&i) //宣告等號重載運算值
    {
        cx=i.cx;
        cy=i.cy;
        ch=i.ch;
        Q=i.Q;
        grid_nums=i.grid_nums;
        int o;
        C=new double [grid_nums];
        for(o=0; o<grid_nums; o++)
            C[o]=i.C[o];
    }
    chimney operator+(chimney Chim) //宣告加號重載運算值
    {
        chimney calC=*this; //宣告新的物件calC,並用複製建構值將他初始化等於呼叫加號重載運算值的物件
        for(int i=0; i<grid_nums; i++) //用for迴圈使舊的濃度值存入新物件的濃度格中
        {
            calC.C[i]=calC.C[i]+Chim.C[i];
        }
        return calC;  //回傳calC物件
    }
    double *getC()
    {
        return C;
    }
    void cal(double *X,double *Y,double *Z,double theta,double &cx,double &cy,double &cz,double Q,double &u)
    {
        int i;
        double sigma_y, sigma_z, a, b, c, d;
        double *x_shift, *y_shift;
        x_shift=new double [grid_nums];
        y_shift=new double [grid_nums];
        for(i=0; i<grid_nums; i++) //平移旋轉座標軸
        {
            x_shift[i]=cos(theta)*(X[i]-cx)+sin(theta)*(Y[i]-cy);
            y_shift[i]=-sin(theta)*(X[i]-cx)+cos(theta)*(Y[i]-cy);
        }
        for(i=0; i<grid_nums; i++) //Steady–State Gussian Plume Model
        {
            if(x_shift[i]>0) //confirm the position is positive
            {
                sigma_y=(0.22*x_shift[i]/sqrt(1+0.0001*x_shift[i]));
                sigma_z=0.001*x_shift[i];
                a=(Q/(2*M_PI*u*sigma_y*sigma_z));
                b=(-y_shift[i]*y_shift[i]/(2*sigma_z*sigma_z));
                c=(-(Z[i]-cz)*(Z[i]-cz)/(2*sigma_z*sigma_z));
                d=(-(Z[i]+cz)*(Z[i]+cz)/(2*sigma_z*sigma_z));
                C[i]=C[i]+a*exp(b)*(exp(c)+exp(d));
            }
        }
    }
private: //宣告私有資料成員
    double *C; //宣告濃度指標陣列C>>存放所有濃度值
    int grid_nums; //宣告整數total>>存放所有值
};

void Max_Concentration(double *C, int &tt, double &maxC, int i) //使用 recursive 算最大濃度值
{
    if(i<tt)
    {
        if(maxC<C[i])
        {
            maxC=C[i];
        }
        Max_Concentration(C,tt,maxC,i+1); //recursive
    }
}

int main()
{
    double x, y, z, dx, dy, dz, wind_direction, wind_speed, theta, u;
    double incx ,incy ,inch, inq, cx, cy, cz, Q, maxC;
    int i, t, w, l,nc, tt;
    double *X, *Y, *Z, *CC;
    ifstream input_file("input.txt"); //open the file
    input_file >> x; //x distance
    input_file >> y;
    input_file >> z;
    input_file >> dx; //grid size
    input_file >> dy;
    input_file >> dz;
    input_file >> wind_speed;
    input_file >> wind_direction;
    range area(x,y,z,dx,dy,dz,wind_direction,wind_speed); //build the range class object <area>
    X=area.getX();
    Y=area.getY();
    Z=area.getZ();
    t=area.gett();
    w=area.getw();
    l=area.getl();
    CC=new double [t];
    theta=area.gettheta();
    u=area.getu();
    chimney *C;  //宣告chimney物件,指標*C
    chimney totalC(t);
    chimney calCC(t);
    input_file >> nc;
    C=new chimney [nc];  //給指標物件CC所有數
    for(i=0; i<nc; i++)
    {
        input_file >> incx;
        input_file >> incy;
        input_file >> inch;
        input_file >> inq;
        C[i].init(incx,incy,inch,inq,t); //呼叫chimney.C中的設值函式init
        cx=C[i].getcx();
        cy=C[i].getcy();
        cz=C[i].getch();
        Q=C[i].getQ();
        C[i].cal(X,Y,Z,theta,cx,cy,cz,Q,u);
        if(i==0)
        {
            calCC=C[i];
            totalC=calCC;

        }
        else
        {
            calCC=totalC;
            totalC=calCC+C[i];
        }
        CC=totalC.getC();
    }
    input_file.close();
    int a;
    ofstream output_file("output.csv");
    output_file <<"distance" << "," << "concentration" << endl;
    for(i=0; i<w; i++)
    {
        a=i*52;
        output_file << sqrt((X[a]*X[a])+(Y[a]*Y[a])) << "," << CC[a] << endl;
    }
    output_file.close();
    i=0; //設i為零以控制遞迴
    maxC=0.0;
    tt=l*w;
    Max_Concentration(CC,tt,maxC,i);// maximum of ground concentration
    cout << "Maximum of ground concentration = " << maxC << endl;
    system("pause");
    return 0;
}



