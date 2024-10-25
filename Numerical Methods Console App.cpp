#include<bits/stdc++.h>
using namespace std;



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
