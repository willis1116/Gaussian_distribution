#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>

using namespace std;

void Dispersion_Model(double *X,double *Y,double *Z,int &t,double *C,double theta,double &cx,double &cy,double &cz,double Q,double &u)
{
    int i;
    double sigma_y, sigma_z, a, b, c, d;
    double *x_shift, *y_shift;
    x_shift=new double [t];
    y_shift=new double [t];
    for(i=0;i<t;i++) //平移旋轉座標軸
    {
        x_shift[i]=cos(theta)*(X[i]-cx)+sin(theta)*(Y[i]-cy);
        y_shift[i]=-sin(theta)*(X[i]-cx)+cos(theta)*(Y[i]-cy);
    }
    for(i=0; i<t; i++) //Steady–State Gussian Plume Model
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
    int i, j, k, t, n, int_w, tt;
    double x, y, z, dx, dy, dz, l, w, h, maxC;
    double *X,*Y,*Z,*A,*B,*C;
    double  nc, theta,  u, cx, cy, cz, Q;
    ifstream input("input.txt"); //open the file and give the variable name
    input >> x;
    input >> y;
    input >> z;
    input >> dx;
    input >> dy;
    input >> dz;
    l=(x/dx)+1;
    w=(y/dy)+1;
    h=(z/dz)+1;
    t=l*w*h;
    tt=l*w;
    int_w=int(w);
    X=new double [t];
    Y=new double [t];
    Z=new double [t];
    C=new double [t]();
    A=new double [int_w];
    B=new double [int_w];
    for(k=0;k<h;k++) // build the grids, from the bottom to top
    {
        for(j=0;j<w;j++)
        {
            for(i=0;i<l;i++)
            {
                n=i+l*j+l*w*k;
                X[n]=dx*i;
                Y[n]=dy*j;
                Z[n]=dz*k;
            }
        }
    }
    input >> u;
    input >> theta;
    input >> nc;
    theta=(theta*M_PI)/180;
    for(i=0;i<nc;i++)
    {
        input >> cx;
        input >> cy;
        input >> cz;
        input >> Q;
        Dispersion_Model(X,Y,Z,t,C,theta,cx,cy,cz,Q,u);
    }
    input.close();

    j=0;
    for(i=0;i<int_w;i++) //方形區域45度對角線
    {
        A[i]=sqrt((X[j]*X[j])+(Y[j]*Y[j]));
        B[i]=C[j];
        j=j+52;
    }
    ofstream output_file("output.csv");
    output_file <<"distance" << "," << "concentration" << endl;
    for(i=0;i<int_w;i++)
    {
         output_file << A[i] << "," << B[i] << endl;
    }
    output_file.close();
    i=0;
    maxC=0;
    Max_Concentration(C,tt,maxC,i);// maximum of ground concentration
    cout << "Maximum of ground concentration = " << maxC << endl;
    system("pause");
    return 0;
}
