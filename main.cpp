#include "Tool.h"

//范数
double norm(double* k,double num){
    double sum = 0.0;
    for (int i = 0; i < num; ++i) {
        sum += k[i] * k[i];
    }
    return sqrt(sum);
}

//输出
//打印
void ShowOne(double* x,double t,double num){
    cout<<"第"<<left<<setw(2)<<t<<"次迭代结果: ";
    for(int i=0;i<num;i++){
        cout<<"x"<<i+1<<":"<<left<<setw(12)<<x[i];
    }
    cout<<endl;
}

//雅可比迭代公式
void Jacobi(){
    cout<<"\n----------雅可比迭代公式---------"<<endl;
    cout<<"本环节将对增广矩阵为:\n"<<"10 -1 -2 7.2\n-1 10 -2 8.3\n-1 -1  5 4.2\n"<<"的线性方程组进行雅可比迭代公式求解"<<endl;
    double x[3] = {0,0,0};//解
    double tmpx[3] = {0,0,0};//临时
    double epsilon = 0.01;//精度
    cout<<"请输入精度:";
    epsilon = fabs(inputDouble());
    if (epsilon == 0) {
        throw myExpection("精度不能为0");
    }
    ShowOne(x,0,3);
    for (int i = 0; i < 100; ++i) {
        tmpx[0] = 0.1*x[1] + 0.2*x[2] + 0.72;
        tmpx[1] = 0.1*x[0] + 0.2*x[2] + 0.83;
        tmpx[2] = 0.2*x[0] + 0.2*x[1] + 0.84;
        ShowOne(tmpx,i+1,3);
        if (fabs(norm(x,3)- norm(tmpx,3))<epsilon){
            cout<<"达到精度要求，迭代结束"<<endl;
            return;
        }
        for (int j = 0; j < 3; ++j) {
            x[j] = tmpx[j];
        }
    }
    cout<<"迭代次数过多，仍未达到要求精度，迭代失败。"<<endl;
}

//高斯-塞德尔迭代
void Gauss_Seidel(){
    cout<<"\n----------高斯-塞德尔迭代---------"<<endl;
    cout<<"本环节将对增广矩阵为:\n"<<"10 -1 -2 7.2\n-1 10 -2 8.3\n-1 -1  5 4.2\n"<<"的线性方程组进行高斯-塞德尔迭代求解"<<endl;
    double x[3] = {0,0,0};//解
    double length = 0;
    double length2 = 0;
    double epsilon = 0.01;//精度
    cout<<"请输入精度:";
    epsilon = fabs(inputDouble());
    if (epsilon == 0) {
        throw myExpection("精度不能为0");
    }
    ShowOne(x,0,3);
    for (int i = 0; i < 100; ++i) {
        x[0] = 0.1*x[1] + 0.2*x[2] + 0.72;
        x[1] = 0.1*x[0] + 0.2*x[2] + 0.83;
        x[2] = 0.2*x[0] + 0.2*x[1] + 0.84;
        ShowOne(x,i+1,3);
        length2 = norm(x,3);
        if (fabs(length2- length)<epsilon){
            cout<<"达到精度要求，迭代结束"<<endl;
            return;
        }
        length = length2;
    }
    cout<<"迭代次数过多，仍未达到要求精度，迭代失败。"<<endl;
}

//超松弛法
void SOR(){
    cout<<"\n----------超松弛法---------"<<endl;
    cout<<"本环节将对增广矩阵为:\n"<<"10 -1 -2 7.2\n-1 10 -2 8.3\n-1 -1  5 4.2\n"<<"的线性方程组进行超松弛法求解"<<endl;
    double x[3] = {0,0,0};//解
    double length = 0;
    double length2 = 0;
    double epsilon = 0.01;//精度
    double omega = 0;//松弛因子
    cout<<"请输入精度:";
    epsilon = fabs(inputDouble());
    if (epsilon == 0) {
        throw myExpection("精度不能为0");
    }
    cout<<"请输入松弛因子(1-2):";
    omega = inputDouble();
    if (omega<1 || omega>2){
        throw myExpection("请输入范围为1-2的松弛因子。");
    }

    ShowOne(x,0,3);
    for (int i = 0; i < 100; ++i) {
        x[0] = (1-omega)*x[0]+omega*(0.1*x[1] + 0.2*x[2] + 0.72);
        x[1] = (1-omega)*x[1]+omega*(0.1*x[0] + 0.2*x[2] + 0.83);
        x[2] = (1-omega)*x[2]+omega*(0.2*x[0] + 0.2*x[1] + 0.84);
        ShowOne(x,i+1,3);
        length2 = norm(x,3);
        if (fabs(length2- length)<epsilon){
            cout<<"达到精度要求，迭代结束"<<endl;
            return;
        }
        length = length2;
    }
    cout<<"迭代次数过多，仍未达到要求精度，迭代失败。"<<endl;
}

bool jordan(int n, double a[][10]){
    //有n-1行,n列
    for (int i = 0; i < n; i++)
    {
        //判断是否非奇异矩阵
        if (a[i][i] == 0) {
            return false;
        }
        //处理第i行一行
        for (int j = i+1; j <= n; ++j) {
            a[i][j] /= a[i][i];
        }
        a[i][i] = 1.0;
        //处理其他
        for (int i2 = 0; i2 < n; ++i2) {
            if (i2 == i){
                continue;
            }
            double kl = a[i2][i];
            for (int j2 = i; j2 <= n; ++j2) {
                a[i2][j2] -= a[i][j2]*kl;
            }
        }
    }
    //判断Xn是否为有唯一解，否则返回错误
    if (a[n - 1][n - 1] == 0)
        return false;
    return true;
}
void Jordan(){
    cout << "\n----------约当消去法---------" << endl;
    int i, j, n;
    double a[10][10];
    cout << "请输入方程组的未知数的个数n:";
    n = inputInt();
    cout << "请逐个输入方程组的增广矩阵:";
    for (i = 0; i < n; i++) {
        //列数多一列
        for (j = 0; j <= n; j++)
            a[i][j] = inputDouble();
    }
    if (jordan(n, a)) {
        cout << "方程组的解是:" << endl;
        for (i = 0; i < n; i++)
            cout << setw(5) << a[i][n] << " ";
        cout << endl;
    } else {
        cout << "方程无解或有无穷解" << endl;
    }
}

//2 -1 3 1 4 2 5 4 1 2 0 7
//高斯消去法
bool gauss(int n, double a[][10])
{
    //有n-1行,n列
    for (int i = 0; i < n; i++)
    {
        //判断是否非奇异矩阵
        if (a[i][i] == 0) {
            return false;
        }
        //处理第i行
        for (int j = i+1; j <= n; ++j) {
            a[i][j] /= a[i][i];
        }
        a[i][i] = 1.0;
        //处理下面每一行
        for (int i2 = i+1; i2 < n; ++i2) {
            double kl = a[i2][i];
            for (int j2 = i; j2 <= n; ++j2) {
                a[i2][j2] -= a[i][j2]*kl;
            }
        }
    }
    //判断Xn是否为有唯一解，否则返回错误
    if (a[n - 1][n - 1] == 0)
        return false;
    //回代求解
    for (int j = n-1; j >= 0; --j) {
        for (int i = j-1; i>=0 ; --i) {
            //处理b矩阵
            a[i][n] -= a[j][n]*a[i][j]/a[j][j];
            a[i][j] = 0.0;
        }
    }
    return true;
}
void Gauss() {
    cout << "\n----------高斯消去法---------" << endl;
    int i, j, n;
    double a[10][10];
    cout << "请输入方程组的未知数的个数n:";
    n = inputInt();
    cout << "请逐个输入方程组的增广矩阵:";
    for (i = 0; i < n; i++) {
        //列数多一列
        for (j = 0; j <= n; j++)
            a[i][j] = inputDouble();
    }
    if (gauss(n, a)) {
        cout << "方程组的解是:" << endl;
        for (i = 0; i < n; i++)
            cout << setw(5) << a[i][n] << " ";
        cout << endl;
    } else {
        cout << "方程无解或有无穷解" << endl;
    }
}

//列主元高斯消去
bool mainGauss(int n, double a[][10]){
    //有n-1行,n列
    for (int i = 0; i < n; i++)
    {
        //判断是否非奇异矩阵
        if (a[i][i] == 0) {
            return false;
        }
        //处理第i行一行
        for (int j = i+1; j <= n; ++j) {
            a[i][j] /= a[i][i];
        }
        a[i][i] = 1.0;
        //交换
        int maxI = i+1;
        for (int i2 = maxI+1; i2 < n; ++i2) {
            if (a[i2][i]>a[maxI][i]){
                maxI = i2;
            }
        }
        if (maxI!=(i+1)){
            swap(a[maxI],a[i+1]);
        }
        //处理下面每一行
        for (int i2 = i+1; i2 < n; ++i2) {
            double kl = a[i2][i];
            for (int j2 = i; j2 <= n; ++j2) {
                a[i2][j2] -= a[i][j2]*kl;
            }
        }
    }
    //判断Xn是否为有唯一解，否则返回错误
    if (a[n - 1][n - 1] == 0)
        return false;
    //回代求解
    for (int j = n-1; j >= 0; --j) {
        for (int i = j-1; i>=0 ; --i) {
            //处理b矩阵
            a[i][n] -= a[j][n]*a[i][j]/a[j][j];
        }
    }
    return true;
}
void MainGauss(){
    cout << "\n----------列主元高斯消去法---------" << endl;
    int i, j, n;
    double a[10][10];
    cout << "请输入方程组的未知数的个数n:";
    n = inputInt();
    cout << "请逐个输入方程组的增广矩阵:";
    for (i = 0; i < n; i++) {
        //列数多一列
        for (j = 0; j <= n; j++)
            a[i][j] = inputDouble();
    }
    if (mainGauss(n, a)) {
        cout << "方程组的解是:" << endl;
        for (i = 0; i < n; i++)
            cout << setw(5) << a[i][n] << " ";
        cout << endl;
    } else {
        cout << "方程无解或有无穷解" << endl;
    }
}

int main() {
    try {
        cout << "这里是线性方程组求根问题的求解" << endl;
        cout << "您可以随时输入'q'来结束当前任务" << endl;
        cout << "----------请选择模式--------" << endl;
        cout << "----------迭代法--------" << endl;
        cout << "1.雅可比迭代公式" << endl;
        cout << "2.高斯-塞德尔迭代" << endl;
        cout << "3.超松弛法" << endl;
        cout << "----------直接法--------" << endl;
        cout << "4.约当消去法" << endl;
        cout << "5.高斯消去法" << endl;
        cout << "6.列主元高斯消去法" << endl;
        cout << "7.追赶法" << endl;
        cout << "8.平方根法" << endl;
        cout << "--------------------------" << endl;
        cout << "请选择:";
        int choice = inputInt();
        switch (choice) {
            case 1:
                Jacobi();
                break;
            case 2:
                Gauss_Seidel();
                break;
            case 3:
                SOR();
                break;
            case 4:
                Jordan();
                break;
            case 5:
                Gauss();
                break;
            case 6:
                MainGauss();
                break;
            case 7:
                break;
            default:
                throw myExpection("选择错误");
                break;
        }
    }catch(myExpection e){
        cout<<e.what()<<endl;
    }
    return 0;
}