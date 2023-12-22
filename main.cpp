#include "Tool.h"

struct Node{
    double x;
    double y;
    Node(): x(0), y(0){}
    Node(double x, double y): x(x), y(y){}
};

//常微分方程
double func(double x,double y) {
    return y - 2 * x / y;
}

//打印
void Show(Node* n,int num){
    cout<<"求值结果为:"<<endl;
    for (int i = 0; i < num; ++i) {
        cout<<"x"<<setw(3)<<left<<i<<"= "<<setw(10)<<left<<fixed<<setprecision(4)<<n[i].x<<"y"<<setw(3)<<left<<i<<"= "<<left<<setw(10)<<fixed<<setprecision(4)<<n[i].y<<endl;
    }
}

//传值1
void inputInfo(double& x,double& y,double& h,int& num){
    cout<<"请输入初值的横坐标x:";
    x = inputDouble();//横坐标
    cout<<"请输入初值的纵坐标y:";
    y = inputDouble();//纵坐标
    cout<<"为了简化，请输入步长:";
    h = inputDouble();
    cout<<"请输入想要求的节点个数(小于100):";
    num = inputInt();
    string errorS = "输入数据过多,中止。";
    if (num > 100){
        throw myExpection("输入数据过多");
    }
}
//传值2
void inputInfo(double& x1,double& x2,double& y1,double& y2,double& h,int& num){
    cout<<"请输入初值的横坐标x:";
    x1 = inputDouble();//横坐标1
    cout<<"请输入初值的纵坐标y:";
    y1 = inputDouble();//纵坐标1
    cout<<"请输入第二个值的横坐标x:";
    x2 = inputDouble();//横坐标2
    cout<<"请输入第二个值的纵坐标y:";
    y2 = inputDouble();//纵坐标2
    h = x2-x1;
    if (h<0){
        throw myExpection("第二个横坐标要大于第一个横坐标");
    }
    cout<<"经计算，步长为"<<h<<endl;
    cout<<"请输入想要求的节点个数(小于100):";
    num = inputInt();
    if (num > 99 || num < 3){
        throw myExpection("输入数据数量错误");
    }
}
//传值4
void inputInfo(double* x,double* y,double& h,int& num){
    cout<<"为了简化，请输入步长:";
    h = inputDouble();
    cout<<"请输入第1个节点的横坐标x:";
    x[0] = inputDouble();//横坐标1
    for (int i = 1; i < 4; ++i) {
        x[i] =  x[i-1]+h;
    }
    for (int i = 0; i < 4; ++i) {
        cout<<"请输入第"<<i+1<<"个节点的纵坐标y:";
        y[i] = inputDouble();//纵坐标1
    }
    cout<<"请输入想要求的节点个数(小于100):";
    num = inputInt();
    if (num > 99 || num < 5){
        throw myExpection("输入数据数量错误");
    }
}

//欧拉格式
void Euler(){
    cout<<"\n--------欧拉格式法求值--------"<<endl;
    cout<<"本环节将对(y'=y-2x/y   y(0)=1)的常微分方程进行欧拉求值"<<endl;
    Node n[100];
    double x,y,h;
    int num;
    inputInfo(x,y,h,num);

    //欧拉格式
    n[0] = Node(x,y);
    for (int i = 1; i < num; ++i) {
        x += h;
        y = n[i - 1].y + h * func(n[i - 1].x, n[i - 1].y);
        n[i] = Node(x, y);
    }

    //输出
    Show(n,num);
}

//两步欧拉
void TwoSteps_Euler(){
    cout<<"\n-------两步欧拉格式法求值------"<<endl;
    cout<<"本环节将对(y'=y-2x/y   y(0)=1)的常微分方程进行两步欧拉求值"<<endl;
    Node n[100];
    double x1,x2,y1,y2,h;
    int num;
    inputInfo(x1,x2,y1,y2,h,num);

    //两步欧拉格式
    n[0] = Node(x1,y1);
    n[1] = Node(x2,y2);
    double x = x2,y;
    for (int i = 2; i < num; ++i) {
        x += h;
        y = n[i - 2].y + 2 * h * func(n[i - 1].x, n[i - 1].y);
        n[i] = Node(x, y);
    }

    //输出
    Show(n,num);
}

//改进欧拉
void Improved_Euler(){
    cout<<"\n------改进欧拉格式法求值------"<<endl;
    cout<<"本环节将对(y'=y-2x/y   y(0)=1)的常微分方程进行欧拉求值"<<endl;
    Node n[100];
    double x,y,h;
    int num;
    inputInfo(x,y,h,num);

    //改进欧拉格式
    n[0] = Node(x,y);
    for (int i = 1; i < num; ++i) {
        x += h;
        y = n[i - 1].y + h * func(n[i - 1].x, n[i - 1].y);
        y = n[i - 1].y + (h/2)*(func(n[i - 1].x, n[i - 1].y) + func(x, y));
        n[i] = Node(x, y);
    }

    //输出
    Show(n,num);
}

//四阶龙格-库塔
void Classical_Runge_Kutta(){
    cout<<"\n------四阶龙格-库塔经典格式------"<<endl;
    cout<<"本环节将对(y'=y-2x/y   y(0)=1)的常微分方程进行四阶龙格-库塔经典格式求值"<<endl;
    Node n[100];
    double x,y,h;
    int num;
    inputInfo(x,y,h,num);

    //四阶龙格-库塔经典格式
    double k1,k2,k3,k4;
    n[0] = Node(x,y);
    for (int i = 1; i < num; ++i) {
        k1 = func(x,y);
        k2 = func(x+0.5*h,y+h/2*k1);
        k3 = func(x+0.5*h,y+h/2*k2);
        k4 = func(x+h,y+h*k3);
        y = y + h/6*(k1+2*k2+2*k3+k4);
        x+=h;
        n[i] = Node(x,y);
    }

    //输出
    Show(n,num);
}

//亚当姆斯预报-校正系统
void Adams_Calibration_System(){
    cout<<"\n------亚当姆斯预报-校正系统------"<<endl;
    cout<<"本环节将对(y'=y-2x/y   y(0)=1)的常微分方程进行亚当姆斯预报-校正系统求值"<<endl;
    Node n[100];
    double x[4],y[4],h;
    int num;
    double tmpx,tmpy,tmpyfy;

    inputInfo(x,y,h,num);

    //测试
//    x[0]=0;x[1]=0.1;x[2]=0.2;x[3]=0.3;
//    y[0]=1.0;y[1]=1.0954;y[2]=1.1832;y[3]=1.2649;
//    h=0.1;
//    num=11;

    //亚当-姆斯校正系统
    for (int i = 0; i < 4; ++i) {
        n[i] = Node(x[i],y[i]);
    }

    for (int i = 4; i < num; ++i) {
        tmpx = n[i - 1].x + h;
        tmpy = (n[i - 1].y) + h * ((55 * func(n[i - 1].x, n[i - 1].y)) - (59 * func(n[i - 2].x, n[i - 2].y)) + (37 * func(n[i - 3].x, n[i - 3].y)) - (9 * func(n[i - 4].x, n[i - 4].y))) / 24.0;
        tmpyfy = func(tmpx, tmpy);
        tmpy = n[i - 1].y + (h / 24.0) * (9 * tmpyfy + 19 * func(n[i - 1].x, n[i - 1].y) - 5 * func(n[i - 2].x, n[i - 2].y) + func(n[i - 3].x, n[i - 3].y));
        n[i] = Node(tmpx, tmpy);
    }

    //输出
    Show(n,num);
}


int main() {
    try {
        cout << "这里是常微分方程初值问题的求解" << endl;
        cout << "您可以随时输入'q'来结束当前任务" << endl;
        cout << "----------请选择模式--------" << endl;
        cout << "1.欧拉格式" << endl;
        cout << "2.两步欧拉格式" << endl;
        cout << "3.改进欧拉格式" << endl;
        cout << "4.四阶龙格-库塔经典格式" << endl;
        cout << "5.亚当-姆斯校正系统" << endl;
        cout << "--------------------------" << endl;
        cout << "请选择:";
        int choice = inputInt();
        switch (choice) {
            case 1:
                Euler();
                break;
            case 2:
                TwoSteps_Euler();
                break;
            case 3:
                Improved_Euler();
                break;
            case 4:
                Classical_Runge_Kutta();
                break;
            case 5:
                Adams_Calibration_System();
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
