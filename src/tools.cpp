#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    
    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if(estimations.size() != ground_truth.size() || estimations.size() == 0) {
        return rmse;
    }
    
    //accumulate squared diffs
    for(int i=0; i < estimations.size(); ++i) {
        VectorXd currDiff = estimations[i] - ground_truth[i];
        
        currDiff = currDiff.array()*currDiff.array();
        
        rmse += currDiff;
    }
    
    //calculate the mean
    rmse /= estimations.size();
    
    //calculate the squared root
    rmse = rmse.array().sqrt();
    
    //return the result
    return rmse;
}


VectorXd Tools::CartesianToPolarCTRV(const VectorXd& measurement) {
    
    VectorXd polar(3);
    
    const float px = measurement(0);
    const float py = measurement(1);
    const float v = measurement(2);
    const float yaw = measurement(3);
    
    const float v1 = cos(yaw)*v;
    const float v2 = sin(yaw)*v;
    
    const float rho = hypot(px, py);
    const float phi = atan2(py, px);
    const float rho_dot = (rho != 0) ? (px * v1 + py * v2) / rho : 0;
    
    polar << rho, phi, rho_dot;
    return polar;
}

VectorXd Tools::PolarToCartesian(const VectorXd& measurement) {
    
    VectorXd cartesian(4);
    
    const double rho = measurement(0);
    const double phi = measurement(1);
    const double rho_dot = measurement(2);
    
    const double px = rho * cos(phi);
    const double py = rho * sin(phi);
    const double vx = rho_dot * cos(phi);
    const double vy = rho_dot * sin(phi);
    
    cartesian << px, py, vx, vy;
    return cartesian;
}


float Tools::random_number(float from, float to)
{
    return ((to - from) * ((float)rand() / RAND_MAX)) + from;
}




