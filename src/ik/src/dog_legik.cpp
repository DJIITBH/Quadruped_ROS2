#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <eigen3/Eigen/Dense>
using namespace std;

#define PI 3.14159265358979323846


struct DHParameter {
    double a;
    double alpha;
    double d;
    double theta;
};

struct Cart_Pose {
    double x;
    double y;
    double z;
    double roll;
    double pitch;
    double yaw;
};

struct Cart_Vel {
    double x_dot;
    double y_dot;
    double z_dot;
    double roll_dot;
    double pitch_dot;
    double yaw_dot;
};

struct Joint_Pose
{
    double Theta_1;
    double Theta_2;
    double Theta_3;
    double Theta_4;
    double Theta_5;
    double Theta_6;

};

struct Joint_Vel
{

double q1_dot;
double q2_dot;
double q3_dot;
double q4_dot;
double q5_dot;
double q6_dot;
};

using Matrix4x4 = std::array<std::array<double, 4>, 4>;
using Vector3 = std::array<double, 3>;
using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

Matrix transposeMatrix(const Matrix& mat) {
    Matrix transposed(mat[0].size(), std::vector<double>(mat.size()));
    for (size_t i = 0; i < mat.size(); ++i) {
        for (size_t j = 0; j < mat[0].size(); ++j) {
            transposed[j][i] = mat[i][j];
        }
    }
    return transposed;
}



Matrix multiplyMatrices(const Matrix& A, const Matrix& B) {
    Matrix result(A.size(), std::vector<double>(B[0].size(), 0.0));
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < B[0].size(); ++j) {
            for (size_t k = 0; k < A[0].size(); ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}


Vector multiplyMatrixVector(const Matrix& mat, const Vector& vec) {
    Vector result(mat.size(), 0.0);
    for (size_t i = 0; i < mat.size(); ++i) {
        for (size_t j = 0; j < vec.size(); ++j) {
            result[i] += mat[i][j] * vec[j];
        }
    }
    return result;
}


Matrix dhTransform(double a, double alpha, double d, double theta) {
    Matrix T =  {{
        {cos(theta), -sin(theta) * cos(alpha),  sin(theta) * sin(alpha), a * cos(theta)},
        {sin(theta),  cos(theta) * cos(alpha), -cos(theta) * sin(alpha), a * sin(theta)},
        {0,           sin(alpha),               cos(alpha),               d},
        {0,           0,                        0,                        1}
    }};
    return T;
}

Vector3 crossProduct(const Vector3& a, const Vector3& b) {
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    };
}

Vector3 subtractVectors(const Vector3& a, const Vector3& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

Matrix forwardKinematicsTransform(const std::vector<DHParameter>& dhParams) {
    int count=0;
    Matrix T = {{
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    }};
    
    for (const auto& param : dhParams) {
        T = multiplyMatrices(T, dhTransform(param.a, param.alpha, param.d, param.theta));
    }
    return T; 
}

    
Cart_Pose forwardKinematics(const std::vector<DHParameter>& dhParams) {
    Matrix T = forwardKinematicsTransform(dhParams);

    Cart_Pose pose;
    pose.x = T[0][3];
    pose.y = T[1][3];
    pose.z = T[2][3];

    // Extracting roll, pitch, yaw from rotation matrix T
    pose.roll = atan2(T[2][1], T[2][2]);
    pose.pitch = atan2(-T[2][0], sqrt(T[2][1] * T[2][1] + T[2][2] * T[2][2]));//0
    pose.yaw = atan2(T[1][0], T[0][0]);
    return pose;
}

Matrix Inverse_Matrix(Matrix&T){
    Matrix T_inv = {{
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    }};
    Eigen::MatrixXd ans(4,4);
    for(int i =0; i<4;i++){
        for(int j=0;j<4;j++){
            ans(i,j)=T[i][j];
        }
    }
    // Code for pseudo-inverse calculation 
    Eigen::MatrixXd inverse = ans.inverse();
     for(int i =0; i<4;i++){
        for(int j=0;j<4;j++){
            T_inv[i][j]=inverse(i,j);
        }
    }
    return T_inv;
}


std::vector<std::vector<double>> calculateJacobian(const std::vector<DHParameter>& dhParams) {
    std::vector<Matrix> transforms;
    Matrix T = {{
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    }};

    for (const auto& param : dhParams) {
        T = multiplyMatrices(T, dhTransform(param.a, param.alpha, param.d, param.theta));
        transforms.push_back(T);
    }

    Vector3 p_e = {T[0][3], T[1][3], T[2][3]};
    std::vector<std::vector<double>> J(6, std::vector<double>(6, 0.0));

    for (int i = 0; i < 6; ++i) {
        Vector3 z_i = {transforms[i][0][2], transforms[i][1][2], transforms[i][2][2]};
        Vector3 p_i = {transforms[i][0][3], transforms[i][1][3], transforms[i][2][3]};
        Vector3 p_e_minus_p_i = subtractVectors(p_e, p_i);
        Vector3 J_v = crossProduct(z_i, p_e_minus_p_i);

        J[0][i] = J_v[0];
        J[1][i] = J_v[1];
        J[2][i] = J_v[2];

        J[3][i] = z_i[0];
        J[4][i] = z_i[1];
        J[5][i] = z_i[2];
    }

    return J;
}

Matrix pseudoInverse(const Matrix& mat) {
    // pseudo inverse method is computationally time consuming so for real time applications like otg or cartesian jogging it has to be avoided so a new method is required
    Matrix U, V;
    Vector S;
    Matrix transposed = transposeMatrix(mat);
    Matrix ans1(mat.size(),std::vector<double>(mat.size()));
    Matrix ans2(mat.size(),std::vector<double>(mat.size()));
    Matrix result(mat[0].size(), Vector(mat.size()));

    // Pseudo-inverse calculation
   ans1 =  multiplyMatrices(mat,transposed);
    Eigen::MatrixXd ans(ans1.size(),ans1.size());
    for(int i =0; i<ans1.size();i++){
        for(int j=0;j<ans1.size();j++){
            ans(i,j)=ans1[i][j];
        }
    }
    // Code for pseudo-inverse calculation 
    Eigen::MatrixXd inverse = ans.inverse();
     for(int i =0; i<ans2.size();i++){
        for(int j=0;j<ans2.size();j++){
            ans1[i][j]=inverse(i,j);
        }
    }
    result = multiplyMatrices(transposed,ans1);

    return result;
}


void Inverse_position_kinematics(const std::vector<DHParameter>& dhParams, Matrix T){
    //****************/ robot parameters**********
    double d1 = 0.0255;
    double d2 = 0.054;
    double a2 = 0.154;
    double d3 = 0;
    double a3 = 0.125;

    double alpha, beta = 0.0;



    double th11 =0;
    double th13,th14 =0;
    double th2 =0;
    double th3 =0;
    double th4 =0;
    double th51 =0;
    double th52 =0;
    double th53 =0;
    double th54 =0;
    double th61 =0;
    double th62 =0;
    double th63 =0;
    double th64 =0;


    Vector3 p3 = {T[0][3],T[1][3],T[2][3]}; // COORDINATES OF EE
    Vector3 z3 = {T[0][2],T[1][2],T[2][2]}; //DIRECTION OF AXIS Z6 OF EE
    Vector3 x3 = {T[0][0],T[1][0],T[2][0]}; //DIRECTION OF AXIS Z6 OF EE
    Vector3 p2 = {T[0][3]-(a3*x3[0]),T[1][3]-(a3*x3[1]),T[2][3]-(a3*x3[2])}; // COORDINATES OF FRAME 5

    // ************Calculating  th1**************
    // atan(val) -> has range (-pi/2,pi/2)
    // ************While atan2(y,x) has range (-pi,pi) and takes into account sign of x,y!******************
    if (p3[0]>0 && p3[1]>0){
        beta = atan(d2/a2);
        alpha = atan(p2[1]/p2[0]);
        th11 = alpha + beta;
        cout<<"th11 "<<th11<<endl;
    }
    
    else if(p3[0]<0 && p3[1]>0){
        beta = atan(d2/a2);
        alpha = atan(abs(p2[1]/p2[0]));
        th11 = alpha + beta;
        cout<<"th11s "<<th11<<endl;
    }
    else if(p3[0]<0 && p3[1]<0){
        beta = atan(d2/a2);
        alpha = PI - atan(abs(p2[1]/p2[0]));
        th13 = -(alpha - beta);      //check it again
        cout<<"th13 "<<th13<<endl;
    }
     else if(p3[0]>0 && p3[1]<0){
        beta = atan(d2/a2);
        alpha = atan(abs(p2[1]/p2[0]));
        th14 = -(alpha - beta);
        cout<<"th14 "<<th14<<endl;
    }
    // double th15 = acos(p2[0]/sqrt(p2[0]*p2[0] + p2[1]*p2[1]));
    // cout<<"th15 "<<th15<<endl;
//    ***************************th2*********************
    th2 = atan((p2[2]- d1)/sqrt(p2[0]*p2[0] + p2[1]*p2[1]));
    cout<<"th2 "<<th2<<endl;

    // ************************th3*********************
    th3 = atan((p3[2]-p2[2])/sqrt((p3[0]-p2[0])*(p3[0]-p2[0])+(p3[1]-p2[1])*(p3[1]-p2[1]))) - th2;
    cout<<"th3 "<<th3<<endl;
    
    

}

// class NumberPublisherNode : public rclcpp::Node
// {
// public:
//     NumberPublisherNode() : Node("number_publisher")
//     {
//         number_publisher_= this->create_publisher<example_interfaces::msg::Int64>("number", 10);
//         timer_=this->create_wall_timer(chrono::seconds(1), bind(&NumberPublisherNode::CounterCallback, this));
//         RCLCPP_INFO(this->get_logger(), "number_publisher has been started");
//     }
// private:
//     void CounterCallback(){
//         auto Int = example_interfaces::msg::Int64();
//         Int.data = counter_;
//         number_publisher_-> publish(Int);

//     }

//     rclcpp::Publisher<example_interfaces::msg::Int64>::SharedPtr number_publisher_;
//     rclcpp::TimerBase::SharedPtr timer_;
//     int counter_=2;

// };


int main() {
    Joint_Pose jpose;
    jpose = {0.23, 0.44, 1.3};
    std::vector<DHParameter> dhParams = {
        {0,     PI / 2,    0.0255,    jpose.Theta_1},
        {0.154,    0,         0.054,      jpose.Theta_2},
        {0.125,   0,    0,   jpose.Theta_3}
    };

    Cart_Pose pose = forwardKinematics(dhParams);
    Matrix T = forwardKinematicsTransform(dhParams);
    std::cout << "Position: (" << pose.x << ", " << pose.y << ", " << pose.z << ")\n";
    // std::cout << "Orientation: (roll: " << pose.roll << ", pitch: " << pose.pitch << ", yaw: " << pose.yaw << ")\n";
    Inverse_position_kinematics(dhParams,T);
    
    // rclcpp::init(argc, argv);
    // auto node = make_shared<NumberPublisherNode>();
    // rclcpp::spin(node);
    // rclcpp::shutdown();
    return 0;
}