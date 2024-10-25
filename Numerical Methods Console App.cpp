#include<bits/stdc++.h>
using namespace std;

void guass_seidel()
{
    cout << " Enter No of variables: ";
    int var; 
    cin >> var;
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
