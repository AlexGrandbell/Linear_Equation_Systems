//
// Created by AlexGrandbell on 2023.11.02.
//

#include "Tool.h"

double inputDouble(){
    return doubleNumberTest(inputString());
}

int inputInt(){
    return intNumberTest(inputString());
}

string inputString(){
    string tmp;
    cin>>tmp;
    exitTest(tmp);
    return tmp;
}

//double转换
double doubleNumberTest(string num){
    double tempDouble = 0;
    istringstream iss(num);
    if (iss >> tempDouble) {
        return tempDouble;
    } else {
        throw myExpection("在需要输入小数的地方输入的不是一个数字，操作失败！");
    }
};

int intNumberTest(string num){
    int tempDouble = 0;
    std::istringstream iss(num);
    if (iss >> tempDouble) {
        return tempDouble;
    } else {
        throw myExpection("在需要输入整数的地方输入的不是一个数字，操作失败！");
    }
}

void exitTest(string s){
    if (s == "q"){
        throw myExpection("已手动退出进程。");
    }
}