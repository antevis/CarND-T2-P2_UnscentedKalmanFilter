#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
    /**
    * Constructor.
    */
    Tools();

    /**
    * Destructor.
    */
    virtual ~Tools();

    /**
    * A helper method to calculate RMSE.
    */
    VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);
    
    VectorXd CartesianToPolarCTRV(const VectorXd& measurement);
    
    //Credits to Mithi Sevilla: https://github.com/mithi/fusion-ekf/blob/master/src/tools.cpp
    VectorXd PolarToCartesian(const VectorXd& measurement);
    
    float random_number(float from, float to);

};

#endif /* TOOLS_H_ */
