#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"


using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
    
public:
    
    UKF(); ///* Constructor
    virtual ~UKF(); ///* Destructor

    ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    VectorXd x_;
    
    /**
    * ProcessMeasurement
    * @param meas_package The latest measurement data of either radar or laser
    */
    void ProcessMeasurement(MeasurementPackage meas_package);
    
    ///* Radar NIS
    float nis_radar_;
    
    ///* Lidar NIS
    float nis_laser_;
    
private:
    
    ///* initially set to false, set to true in first call of ProcessMeasurement
    bool is_initialized_;
    
    ///* if this is false, laser measurements will be ignored (except for init)
    bool use_laser_;
    
    ///* if this is false, radar measurements will be ignored (except for init)
    bool use_radar_;
    
    ///* state covariance matrix
    MatrixXd P_;
    
    ///* predicted sigma points matrix
    MatrixXd Xsig_pred_;
    
    ///* time when the state is true, in us
    //Don't know what this is.
    //long long time_us_;
    
    /*
     MEASUREMENT MOISE. Since it's not suppused to be tuned, declared them as constants
     and initialized in header (which is acceptable).
     */
    ///* Laser measurement noise standard deviation position1 in m
    const double std_laspx_ = .15;
    ///* Laser measurement noise standard deviation position2 in m
    const double std_laspy_ = .15;
    ///* Radar measurement noise standard deviation radius in m
    const double std_radr_ = .3;
    ///* Radar measurement noise standard deviation angle in rad
    const double std_radphi_ = .03;
    ///* Radar measurement noise standard deviation radius change in m/s
    const double std_radrd_ = .3;
    
    ///* PROCESS NOISE (to be tuned)
    ///* Process noise standard deviation longitudinal acceleration in m/s^2
    double std_a_;
    ///* Process noise standard deviation yaw acceleration in rad/s^2
    double std_yawdd_;
    
    ///* Weights of sigma points
    VectorXd weights_;
    
    ///*Augmented state dimension
    const int n_aug_ = 7;
    
    ///* Sigma point spreading parameter
    double lambda_ = 3 - n_aug_;
    
    ///* State dimension
    const int n_x_ = 5;
    
    ///*Previous timestamp
    long long prev_t_;
    
    Tools tools;
    
    MatrixXd H_laser_;
    MatrixXd R_laser_; // Laser measurement noise covariance matrix
    MatrixXd R_radar_; // Radar measurement noise covariance matrix
    
    /**
     * Prediction Predicts sigma points, the state, and the state covariance
     * matrix
     * @param delta_t Time between k and k+1 in s
     */
    void Prediction(double delta_t);
    
    /**
     * Updates the state and the state covariance matrix using a laser measurement
     * @param meas_package The measurement at k+1
     */
    void UpdateLidar(MeasurementPackage meas_package);
    
    /**
     * Updates the state and the state covariance matrix using a radar measurement
     * @param meas_package The measurement at k+1
     */
    void UpdateRadar(MeasurementPackage meas_package);
    
    void PredictSigmaPoints(MatrixXd &Xsig_aug, double delta_t);
    
    void ComputeNewEstimate();
    
    void WeightSetup();
    
    MatrixXd AugSigmaPoints();
    
    struct SigmaPoint {
        double p_x;
        double p_y;
        double v;
        double yaw;
        double yawd;
        double nu_a;
        double nu_yawdd;
    };
    
    SigmaPoint GetSigmaPoint(MatrixXd& Xsig, char colum_idx);
    void SetSigmaPoint(SigmaPoint& point, char column_idx, MatrixXd& Xsig);
};

#endif /* UKF_H */


