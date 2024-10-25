#include<bits/stdc++.h>
using namespace std;
float *coeff;
int degree=0;
int a=0,b=0,c=0;

double t_func(x)
{
    return a * sin(x) + b * cos(x) + c * tan(x);
}

double t_derivative(x)
{
    return a * cos(x) - b * sin(x) + c * (1 / pow(cos(x), 2));
}

void guass_seidel(int var)
{
    float matx [var][var+ 1];
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


int main()
{
    int w = 80; 
    string message = "\033[1mKUET Numerical Methods Console Application\033[0m"; 
    cout << setw((w + message.length() - 9) / 2) << message << endl;
    int s;
    cout<<"\033[1;32m-->\033[0m For \033[1;33mSolution of Linear Equations\033[0m select \033[1;33m1\033[0m"<<endl;
    cout<<"\033[1;32m-->\033[0m For \033[1;33mSolution of Non-Linear Equations\033[0m select \033[1;33m2\033[0m"<<endl;
    cout<<"\033[1;32m-->\033[0m For \033[1;33mSolution of Differential Equations\033[0m select \033[1;33m3\033[0m"<<endl;
    cout<<"\033[1;32m-->\033[0m For \033[1;33mMatrix Inversion\033[0m select \033[1;33m4\033[0m\n"<<endl;
    cout<<"\033[1;36m# Enter an option to select:\033[0m"<<endl;
    cin>>s;


    return 0;
}
