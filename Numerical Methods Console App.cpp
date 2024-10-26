#include<bits/stdc++.h>
using namespace std;
typedef long double ld;
float* coeff;
vector<ld> coefficient;
int degree=0;
ld a=0,b=0,c=0,d=0,ft=0;

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
float fx_forNewtonRaphson (float n) {
    float res = 0;

    for (int i = 0; i <= degree; i++) 
        res += coeff[i] * (pow(n, degree - i));
    
    return res;   
}
float fdx_forNewtonRaphson(float n)
{
    float ress=0;
    int dd=degree;
    for(int i = 0; i <degree; i++)
    {
        ress+= (degree - i) * coeff[i] * pow(n, degree - i - 1);
    }
    return ress;
}
void newton_raphsonAlgebric()
{
    coeff = new float[degree + 1];
    int dd=degree;
    for (int i = 0; i <= degree; i++) {
        float Coeff;
        cout << "Coeff: ";
        cin >> Coeff;
        coeff[i] = Coeff;
    }
    float st,end;
    st=sqrt(pow (coeff[1] / coeff[0], 2) - 2 * (coeff[2] / coeff[0]));
    end=-1*st;
    cout<<st<<" "<<end<<endl;
    vector<float>roots;
    cout<<st<<" "<<end<<endl;
    while(dd>0)
    {

    int itr=0;
    cout<<"Enter Initial Guess"<<endl;
    float x1,x,ffx,fxd;
    cin>>x1;
    
    bool  is_duplicate = false;
    do{
         x=x1;
         ffx=fx_forNewtonRaphson(x);
         if(ffx==0){
            break;
         }
     
         fxd=fdx_forNewtonRaphson(x);
         if(fxd==0)
         {
             cout << "Zero derivative; try a different initial guess." << endl;
             break;
         }
         x1=x-(ffx/fxd);
         itr++;
         
         cout<<fabs(x-x1)<<endl;


    }while(fabs(x-x1)>=0.00001);
    for(float root : roots)
    {
        if(root==x1)
        {  is_duplicate=true;
        cout<<"Iteration not counted"<<endl;
        }

    }
     if (!is_duplicate) {
            roots.push_back(x1);
            cout<<"iteration "<<itr<<endl;
            cout << "Root found: " << x1 << endl;
            dd--;
        }
   
    }
    delete[] coeff;
   

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
    float matx [var][var+ 1];
    cout<<"\033[1mEnter the equations below :\033[0m\n";
    for(int i=0;i<var;i++)
    {
        for(int j=0;j<=var;j++)
        {
            cin>>matx[i][j];

        }
    }
    int itr;
    cout<<"Enter  the number of iterations: ";
    cin>>itr;
    cout<<endl;
    vector<float>curV(var,0.0);
    vector<float>prev(var,0.0);
    vector<float>error(var,0.0);
    int count=0;
    while(itr--)
    {
      for(int i=0;i<var;i++)
      { 
        float sum=matx[i][var];
        for(int j=0;j<var;j++)
        {
           
            if(j!=i)    
            {
                sum-=matx[i][j]*curV[j];
            }

        }
        curV[i]=sum/matx[i][i];
        curV[i]=fabs(curV[i]);
        cout<<curV[i]<<" ";
        error[i]=fabs(curV[i]-prev[i]);
        prev[i]=curV[i];
      }
      cout<<endl;
      if(error[0]<0.00001)break;
      count++;
    }
    cout<<"Total iterations: "<<count;
}
void gaussElimination(int n) {
    vector<vector<double>> matrix(n,vector<double>(n+1));
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
                double factor=matrix[j][i];
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
    vector<vector<double>>matrix(n,vector<double>(n+1));
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
        double pivot=matrix[i][i];
        for (int k=0;k<=n;k++) {
            matrix[i][k]/= pivot;
        }

        for (int j=0;j<n;j++)
        {
            if (j!=i) {
                double factor=matrix[j][i];
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
ld func(ld x)
{
    return a*pow(x,5)+b*pow(x,4)+c*pow(x,3)+d*pow(x,2)+e*x+f;
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
    cout<<"\n";
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
         cout<<"Enter variable no.:\n";
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
              break;
            case 2:
              gauss_seidel(n);
              break;
            case 3:
              gaussElimination(n);
              break;
            case 4:
              gaussJordanElimination(n);
              break;
            case 5:
             //lu
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
              cout<<"\033[1mEnter \033[1;35ma, b, c\033[0m in the \033[1;35m asin(x) + bcos(x) + ctan(x) +d\033[0m :\n";
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
              break;
            case 2:
              false_position();
              break;
            case 3:
              //secant();
              break;
            case 4:
              newton_raphsonAlgebric();
              break;
            case 5:
              break;
            default:
              cout<<error;
              goto non_linear;  
        }
        break;
    default:
          cout<<error;
          goto main;
          break;

         
    }
    return 0;
}
