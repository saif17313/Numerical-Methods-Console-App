#include<bits/stdc++.h>
using namespace std;
typedef long double ld;
vector<ld> coefficient;
int degree=0;
ld a=0,b=0,c=0,d=0;
int ft;
int chk;

ld func(ld x)
{
    if(ft==1){
        ld result = 0.0;
         for (int i = 0; i <= degree; i++) {
        result += coefficient[i] * pow(x, degree-i);
      }  
      return result;  
    }
    else
    {
        return a * sin(x) + b * cos(x) + c * tan(x) + d;
    }
}

ld derivative(ld x)
{
   if(ft==1){
    ld result=0.0;
    for (int i = 0; i <= degree; i++) {
        result += (degree-i)*coefficient[i] * pow(x, degree-i-1);
      } 
      return result; 

   }
   else{
   return a * cos(x) - b * sin(x) + c * (1 / pow(cos(x), 2));
   }
}

void bisection()
{
    ld p=-1000,q=-999;
    while(func(p)*func(q)>=0.0)
    {
       if(func(p)==0.0)
       {
           cout<<"\033[1;31mOne of the value is : \033[0m"<<p<<endl;
           return;
       }
       else if(func(q)==0.0)
       {
         cout<<"\033[1;31mOne of the value is : \033[0m"<<q<<endl;
           return;
       }
       p++;
       q++;
    }
    ld x,x1=0,d1=1000;
    int count=0;
    while(d1>0.00001)
    {
       count++;
       x=(p+q)/2;

    if(func(x)==0.0)
        break;

    if(func(x)*func(p)<0.0)
    {
       q=x;
    }
    else
    {
      p=x;
    }
    d1=fabs(x-x1);
    x1=x;
    }
    cout<<setprecision(10)<<fixed;
    cout<<"\033[1;31mOne of the value is : \033[0m"<<x<<endl;
    cout<<"\033[1mThe total needed iteration is : \033[0m"<<count<<endl;
    coefficient.clear();
    a=b=c=d=0;
    return;
}

void false_position()

{
    ld p=-1000,q=-999;
    while(func(p)*func(q)>=0.0)
    {
       if(func(p)==0.0)
       {
           cout<<"\033[1;31mOne of the value is : \033[0m"<<p<<endl;
           return;
       }
       else if(func(q)==0.0)
       {
         cout<<"\033[1;31mOne of the value is : \033[0m"<<q<<endl;
           return;
       }
       p++;
       q++;
    }
    ld x,x1=0,d=1000;
    int count=0;
    while(d>0.0001)
    {
        count++;
       x=(p*func(q)-q*func(p))/(func(q)-func(p));

    if(func(x)==0.0)
        break;

    if(func(x)*func(p)<0.0)
    {
       q=x;
    }
    else
    {
      p=x;
    }
    d=fabs(x-x1);
    x1=x;
    }
    cout<<setprecision(10)<<fixed;
    cout<<"\033[1;31mOne of the value is : \033[0m"<<x<<endl;
    cout<<"\033[1mThe total needed iteration is : \033[0m"<<count<<endl;
    coefficient.clear();
    a=b=c=d=0;
    return;
}
void gauss_seidel(int var)
{
    ld matx [var][var+ 1];
    cout<<"\033[1mEnter the equations below :\033[0m\n";
    for(int i=0;i<var;i++)
    {
        for(int j=0;j<=var;j++)
        {
            cin>>matx[i][j];

        }
    }
    int itr=50;
  
    cout<<endl;
    vector<ld>curV(var,0.0);
    vector<ld>prev(var,0.0);
    vector<ld>error(var,0.0);
    int count=0;
    const ld threshold = 0.00001;
    while(itr--)
    {
        bool converged = true;
      for(int i=0;i<var;i++)
      { 
        ld sum=matx[i][var];
        for(int j=0;j<var;j++)
        {
           
            if(j!=i)    
            {
                sum-=matx[i][j]*curV[j];
            }

        }
        curV[i]=sum/matx[i][i];
        error[i]=fabs(curV[i]-prev[i]);
         if (error[i] > threshold) converged = false;
        
        prev[i]=curV[i];
      }
      count++;
      if (converged) break;
      
    }
    cout<<"The roots for linear equations\n";
   for(ld root : curV)
   {
    cout<< root <<" ";
   }

    cout<<"\nTotal iterations: "<<count<<endl;
}
void gaussElimination(int n) {
    vector<vector<ld>> matrix(n,vector<ld>(n+1));
    cout<<"Enter the augmented matrix:"<<endl;
    for (int i=0;i<n;i++) {
        for (int j=0;j<=n;j++) {
            cin>>matrix[i][j];
        }
    }
    for (int i=0;i<n;i++)
    {
        int max=i;
        for (int k=i+1;k<n;k++) 
        {
            if (fabs(matrix[k][i])>fabs(matrix[max][i]))
            {
                max=k;
            }
        }
        swap(matrix[i],matrix[max]);
        
        if (fabs(matrix[i][i])<1e-10) 
        {
            cout<<"No unique solution exists"<<endl;
            return;
        }
        for (int k=i+1;k<=n;k++) 
        {
            matrix[i][k]/=matrix[i][i];
        }
        matrix[i][i] = 1.0;

        for (int j=0;j<n;j++) 
        {
            if (j!=i) 
            {
                ld factor=matrix[j][i];
                for (int k = i; k <= n; k++) {
                    matrix[j][k] -= factor * matrix[i][k];
                }
            }
        }
    }

    cout << "Solution:"<<endl;
    for (int i = 0; i < n; i++) {
        cout <<"x"<<i+1<<" = "<<matrix[i][n]<< endl;
    }
}

void gaussJordanElimination(int n) {
    vector<vector<ld>>matrix(n,vector<ld>(n+1));
    cout<<"Enter the augmented matrix:"<<endl;
    for(int i=0;i<n;i++) 
    {
        for (int j=0;j<=n;j++) 
        {
            cin>>matrix[i][j];
        }
    }
    for(int i=0;i<n;i++)
    {
        int max=i;
        for (int k=i+1;k<n;k++)
        {
            if(fabs(matrix[k][i])>fabs(matrix[max][i]))
            {
                max=k;
            }
        }
        swap(matrix[i], matrix[max]);
        if (fabs(matrix[i][i]) < 1e-10) {
            cout << "No unique solution exists."<<endl;
            return;
        }
        ld pivot=matrix[i][i];
        for (int k=0;k<=n;k++) {
            matrix[i][k]/= pivot;
        }

        for (int j=0;j<n;j++)
        {
            if (j!=i) {
                ld factor=matrix[j][i];
                for (int k=0; k<=n;k++) {
                    matrix[j][k]-=factor*matrix[i][k];
                }
            }
        }
    }

    cout << "Solution:\n"<<endl;
    for (int i=0;i<n;i++) {
        cout <<"x"<<i+1<<" = "<<matrix[i][n]<< endl;
    }
}

void inversion(int n) {
    vector<vector<double>>matrix(n,vector<double>(n+n));
        cout<<"Enter elements \033[1;35ma1, a2, ..., an\033[0m\n";

    for(int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            cin>>matrix[i][j];
        }
    }
    for(int i=0;i<n;i++)
    {
         for (int j=n;j<n+n;j++)
        {
            if((j-i)==n)
            matrix[i][j]=1;
            else
            matrix[i][j]=0;
        }
    }

    for(int i=0;i<n;i++)
    {
        int max=i;
        for (int k=i+1;k<n;k++)
        {
            if(fabs(matrix[k][i])>fabs(matrix[max][i]))
            {
                max=k;
            }
        }
        swap(matrix[i], matrix[max]);
        if (fabs(matrix[i][i]) < 1e-10) {
            cout<<"\033[1;31mNo unique solution exists\033[0m"<<endl;
            return;
        }
        double pivot=matrix[i][i];
        for (int k=0;k<n+n;k++) {
            matrix[i][k]/= pivot;
        }

        for (int j=0;j<n;j++)
        {
            if (j!=i) {
                double factor=matrix[j][i];
                for (int k=0; k<n+n;k++) {
                    matrix[j][k]-=factor*matrix[i][k];
                }
            }
        }
    }

    cout << "\033[1;35mThe Inverse Matrix is :\033[0m"<<endl;
    for(int i=0;i<n;i++)
    {
        for(int j=n;j<n+n;j++)
        {
            cout<<matrix[i][j]<<"  ";
        }
        cout<<endl;
    }
}


ld secant()
{
    ld x1=0,x2=5,x3=0,f1,f2;
int count=0;
while (true)
{
    count++; 
     f1=func(x1);
     f2=func(x2);
     if(f2-f1==0)
     {
         cout<<"division by 0,undefined"<<endl;
         return x3;
     }
     x3=x2-((f2*(x2-x1))/(f2-f1));
    
    if(fabs(func(x3))<1e-6)
    break;
    x1=x2;
    x2=x3;
}
    
     cout<<"Total needed iteration for secant method"<< count<< endl;
     return x3;
}

int main()
{
    int w = 80; 
    string message = "\033[1;34mKUET Numerical Methods Console Application\033[0m"; 
    string linear="\033[1;32mMethods for solving linear equations :\033[0m\n\033[1;32m-->\033[0m For \033[1;33mJacobi iterative method\033[0m enter \033[1;33m1\033[0m\n\033[1;32m-->\033[0m For \033[1;33mGauss-Seidel iterative method\033[0m enter \033[1;33m2\033[0m\n\033[1;32m-->\033[0m For \033[1;33mGauss elimination\033[0m enter \033[1;33m3\033[0m\n\033[1;32m-->\033[0m For \033[1;33mGausss-Jordan elimination\033[0m enter \033[1;33m4\033[0m\n\033[1;32m-->\033[0m For \033[1;33mLU factorization\033[0m enter \033[1;33m5\033[0m\n\033[1;32m-->\033[0m For \033[1;33mMain Menu\033[0m enter \033[1;33m0\033[0m\n\033[1;32m-->\033[0m To  \033[1;31mEXIT\033[0m enter \033[1;31m-1\033[0m\n";
    string error="\003[1;31mInvalid Input ! Try again..\033[0m\n";
    string option="\033[1;32m-->\033[0m For \033[1;33mSolution of Linear Equations\033[0m enter \033[1;33m1\033[0m\n\033[1;32m-->\033[0m For \033[1;33mSolution of Non-Linear Equations\033[0m enter \033[1;33m2\033[0m\n\033[1;32m-->\033[0m For \033[1;33mSolution of Differential Equations\033[0m enter \033[1;33m3\033[0m\n\033[1;32m-->\033[0m For \033[1;33mMatrix Inversion\033[0m enter \033[1;33m4\033[0m\n\033[1;32m-->\033[0m To  \033[1;31mEXIT\033[0m enter \033[1;31m-1\033[0m\n";
    string choose="\033[1;36m# Enter an option to select :\033[0m\n";
    string nonlinear="\033[1;32mMethods for solving non-linear equations :\033[0m\n\033[1;32m-->\033[0m For \033[1;33mBi-section method\033[0m enter \033[1;33m1\033[0m\n\033[1;32m-->\033[0m For \033[1;33mFalse position method\033[0m enter \033[1;33m2\033[0m\n\033[1;32m-->\033[0m For \033[1;33mSecant method\033[0m enter \033[1;33m3\033[0m\n\033[1;32m-->\033[0m For \033[1;33mNewton-Raphson method\033[0m enter \033[1;33m4\033[0m\n\033[1;32m-->\033[0m For \033[1;33mMain Menu\033[0m enter \033[1;33m0\033[0m\n\033[1;32m-->\033[0m To  \033[1;31mEXIT\033[0m enter \033[1;31m-1\033[0m\n";
    string inverseoption="\033[1;32m-->\033[0m For \033[1;33mInversion of Matrix\033[0m enter \033[1;33m1\033[0m\n\033[1;32m-->\033[0m For \033[1;33mMain Menu\033[0m enter \033[1;33m0\033[0m\n\033[1;32m-->\033[0m To  \033[1;31mEXIT\033[0m enter \033[1;31m-1\033[0m\n";
    string rangeoption="\033[1;32m-->\033[0m For \033[1;33mRange-Kutta Method\033[0m enter \033[1;33m1\033[0m\n\033[1;32m-->\033[0m For \033[1;33mMain Menu\033[0m enter \033[1;33m0\033[0m\n\033[1;32m-->\033[0m To  \033[1;31mEXIT\033[0m enter \033[1;31m-1\033[0m\n";
    string again="\033[1mDo you want to continue ?\033[0m\nPress \033[1;35m1\033[0m for \033[1;35mYES\033[0m\nPress \033[1;31mAny Integer\033[0m for \033[1;31mNO\033[0m :\n";
    cout << setw((w + message.length() - 9) / 2) << message << endl;
    cout<<"\n";
    int s;
    main:
    cout<<option;
    cout<<choose;
    cin>>s;
    switch(s){
    case -1:
        return 0;
    case 1:
        linearr:
         int l1,n;
         cout<<linear;
         cout<<choose;
         cin>>l1;
         if(l1!=0 && l1!=-1)
         {
         cout<<"Enter \033[1;35mvariable no.\033[0m:\n";
         cin>>n;
         }
         switch(l1)
         {
            case 0:
              goto main;
              break;
            case -1:
              return 0;
            case 1:
              jacobian(n);
              cout<<again;
              cin>>chk;
              if(chk==1)
                goto main;
              else
                return 0;
              break;
            case 2:
              gauss_seidel(n);
              cout<<again;
              cin>>chk;
              if(chk==1)
                goto main;
              else
                return 0;
              break;
            case 3:
              gaussElimination(n);
              cout<<again;
              cin>>chk;
              if(chk==1)
                goto main;
              else
                return 0;
              break;
            case 4:
              gaussJordanElimination(n);
              cout<<again;
              cin>>chk;
              if(chk==1)
                goto main;
              else
                return 0;
              break;
            case 5:
             //lu
             cout<<again;
              cin>>chk;
              if(chk==1)
                goto main;
              else
                return 0;
              break;
            default:
              cout<<error;
              goto linearr;  
         }
         break;
    case 2:
        non_linear:
         int l2;
         cout<<nonlinear;
         cout<<choose;
         cin>>l2;
        if(l2!=-1 && l2!=0)
        {
           func:
           cout<<"\033[1;32mFunction Type :\033[0m\n";
           cout<<"\033[1;31m-->\033[0m For \033[1mAlgebric Function\033[0m enter \033[1m1\033[0m\n";
           cout<<"\033[1;31m-->\033[0m For \033[1mTrigonometric Function\033[0m enter \033[1m2\033[0m\n";
           cout<<choose;
           cin>>ft;
           if(ft!=1 && ft!=2)
           {
             cout<<error;
             goto func;
           }
           if(ft==1)
           {
             cout<<"\033[1mEnter \033[35mdegree\033[0m of the algebric function :\033[0m\n";
             cin>>degree;
             coefficient.resize(degree+1);
             cout<<"\033[1mEnter the co-efficients for the function :\033[0m\n";
             for(int i=0;i<=degree;i++)
             {
                cin>>coefficient[i];
             }
           }
           else
           {
              cout<<"\033[1mEnter \033[1;35ma, b, c, d\033[0m in the \033[1;35m asin(x) + bcos(x) + ctan(x) +d\033[0m :\n";
              cin>>a>>b>>c>>d;
           }
        }
        switch(l2){
            case 0:
              goto main;
              break;
            case -1:
              return 0;
            case 1:
              bisection();
              cout<<again;
              cin>>chk;
              if(chk==1)
                goto main;
              else
                return 0;
              break;
            case 2:
              false_position();
              cout<<again;
              cin>>chk;
              if(chk==1)
                goto main;
              else
                return 0;
              break;
            case 3:
              secant();
              cout<<again;
              cin>>chk;
              if(chk==1)
                goto main;
              else
                return 0;
              break;
            case 4:
              //newton_raphsonAlgebric();
              cout<<again;
              cin>>chk;
              if(chk==1)
                goto main;
              else
                return 0;
              break;
            default:
              cout<<error;
              goto non_linear;  
        }
        break;
        case 3:
        range:
          int l4,n3;
          cout<<rangeoption;
          cout<<choose;
          cin>>l4;
          if(l4!=0 && l4!=-1)
          {

          }
          switch(l4){
            case 1:
              //range_kutta();
              cout<<again;
              cin>>chk;
              if(chk==1)
                goto main;
              else
                return 0;
              break;
            case 0:
               goto main;
               break;
            case -1:
               return 0;
            default:
              cout<<error;
              goto inv;
              break;
              
         }
         break;
        case 4:
        inv:
         int l3,n2;
         cout<<inverseoption;
         cout<<choose;
         cin>>l3;
         if(l3!=0 && l3!=-1)
         {
         cout<<"Enter \033[1;35mn\033[0m for \033[1;35m(n*n)\033[0m matrix:\n";
         cin>>n2;
         }
         else{
            cout<<error;
            goto inv;
         }
         switch(l3){
            case 1:
              inversion(n2);
              cout<<again;
              cin>>chk;
              if(chk==1)
                goto main;
              else
                return 0;
              break;
            case 0:
               goto main;
               break;
            case -1:
               return 0;
            default:
              cout<<error;
              goto inv;
              break;
              
         }
         break;
    default:
          cout<<error;
          goto main;
          break;
    }
    return 0;
}
