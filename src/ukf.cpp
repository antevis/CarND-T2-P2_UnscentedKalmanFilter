#include "ukf.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(n_x_);

    // initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = .449;
    
    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = .77;
    
    is_initialized_ = false;
    
    
    H_laser_ = Eigen::MatrixXd(2, n_x_);
    H_laser_.fill(0.);
    for (int i = 0; i < H_laser_.rows(); ++i) {
        H_laser_(i, i) = 1.;
    }
    
    R_laser_ = Eigen::MatrixXd(2, 2);
    R_laser_ <<
        pow(std_laspx_, 2), 0.,
        0.,                  pow(std_laspy_, 2);
    
    
    R_radar_ = MatrixXd(3, 3);
    R_radar_ <<
        pow(std_radr_, 2),  0,                      0,
        0,                  pow(std_radphi_, 2),    0,
        0,                  0,                      pow(std_radrd_, 2);
    

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    
    if (!is_initialized_) {
        
        x_.fill(0.0);
        char x_size = x_.size();
        P_ = MatrixXd::Identity(x_size, x_size);
        
        
        if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            
            //Initialize state.
            x_(0) = meas_package.raw_measurements_[0];
            x_(1) = meas_package.raw_measurements_[1];
            
            
            
        } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            
            //Convert radar from polar to cartesian coordinates and initialize state.
            const Eigen::VectorXd init_state = tools.PolarToCartesian(meas_package.raw_measurements_);
            
            x_(0) = init_state(0);
            x_(1) = init_state(1);
            
            //Tweaking covariance of x and v improves performance for Dataset 2
            P_(0,0) = .1;
            P_(2,2) = .5;
        }
        
        prev_t_ = meas_package.timestamp_;
        
        WeightSetup();
        
        is_initialized_ = true;
        
        return;
    }
    
    const float dt = (meas_package.timestamp_ - prev_t_) / 1000000.0;	//dt - expressed in seconds
    
    Prediction(dt);
    
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        UpdateLidar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        UpdateRadar(meas_package);
    }
    
    prev_t_ = meas_package.timestamp_;
}


MatrixXd UKF::AugSigmaPoints() {
    VectorXd x_aug = VectorXd(7);
    
    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    
    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    Xsig_aug.fill(0.0);
    
    x_aug.fill(0.0);
    x_aug.head(n_x_) = x_;
    
    for (int i = n_x_; i< n_aug_; ++i) {
        x_aug(i) = 0;
    }
    
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    
    P_aug(n_x_, n_x_) = pow(std_a_, 2);
    P_aug(n_x_+1, n_x_+1) = pow(std_yawdd_, 2);
    
    MatrixXd A = P_aug.llt().matrixL();
    
    //create augmented sigma points
    const double sqrt_lbd_nx = sqrt(lambda_ + n_aug_);
    Xsig_aug.col(0) = x_aug;
    
    for (int i = 1; i <= n_aug_; ++i)
    {
        Xsig_aug.col(i) = x_aug + sqrt_lbd_nx * A.col(i-1);
        Xsig_aug.col(i+n_aug_) = x_aug - sqrt_lbd_nx * A.col(i-1);
    }
    return Xsig_aug;
}

UKF::SigmaPoint UKF::GetSigmaPoint(Eigen::MatrixXd& Xsig, char colum_idx) {
    //extract values for better readability
    const VectorXd sigma_point = Xsig.col(colum_idx);
    SigmaPoint pt = SigmaPoint();
    
    pt.p_x = sigma_point(0);;
    pt.p_y = sigma_point(1);;
    pt.v = sigma_point(2);;
    pt.yaw = sigma_point(3);;
    pt.yawd = sigma_point(4);
    pt.nu_a = sigma_point(5);
    pt.nu_yawdd = sigma_point(6);
    
    return pt;
}

void UKF::SetSigmaPoint(SigmaPoint& point, char column_idx, Eigen::MatrixXd& Xsig) {
    Xsig(0, column_idx) = point.p_x;
    Xsig(1, column_idx) = point.p_y;
    Xsig(2, column_idx) = point.v;
    Xsig(3, column_idx) = point.yaw;
    Xsig(4, column_idx) = point.yawd;
    
    if (Xsig.rows() > 5) {
        Xsig(5, column_idx) = point.nu_a;
        Xsig(6, column_idx) = point.nu_yawdd;
    }
    
    
}

void UKF::PredictSigmaPoints(MatrixXd &Xsig_aug, double delta_t) {
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    
    //predict sigma points
    for (int i = 0; i< 2*n_aug_+1; i++)
    {
        SigmaPoint source_point = GetSigmaPoint(Xsig_aug, i);
        SigmaPoint predicted_point;
        
        const float sin_yaw = sin(source_point.yaw);
        const float cos_yaw = cos(source_point.yaw);
        const float dt_squared = pow(delta_t, 2);
        
        //Noises
        const float v_noise = source_point.nu_a * delta_t;
        const float yaw_noise = 0.5 * source_point.nu_yawdd * dt_squared;
        const float yawd_noise = source_point.nu_yawdd * delta_t;
        const float px_noise = 0.5 * source_point.nu_a * dt_squared * cos_yaw;
        const float py_noise = 0.5 * source_point.nu_a * dt_squared * cos_yaw;
        
        //predicted state values
        //Adding noise inline
        predicted_point.v = source_point.v                             + v_noise;
        predicted_point.yaw = source_point.yaw + source_point.yawd * delta_t  + yaw_noise;
        predicted_point.yawd = source_point.yawd                       + yawd_noise;
        
        //avoid exploding
        if (fabs(source_point.yawd) > 0.001) {
            
            predicted_point.p_x = source_point.p_x + source_point.v / source_point.yawd * (sin(source_point.yaw + source_point.yawd * delta_t) - sin_yaw);
            predicted_point.p_y = source_point.p_y + source_point.v / source_point.yawd * (cos_yaw - cos(source_point.yaw + source_point.yawd * delta_t));
        }
        else {
            predicted_point.p_x = source_point.p_x + source_point.v * delta_t * cos_yaw;
            predicted_point.p_y = source_point.p_y + source_point.v * delta_t * sin_yaw;
        }
        
        //add noise
        predicted_point.p_x += px_noise;
        predicted_point.p_y += py_noise;
        
        SetSigmaPoint(predicted_point, i, Xsig_pred_);
    }
}


void UKF::WeightSetup() {
    weights_ = VectorXd(2*n_aug_+1);
    // set weights
    double weight_0 = lambda_/(lambda_+n_aug_);
    weights_(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
        double weight = 0.5/(n_aug_+lambda_);
        weights_(i) = weight;
    }
}

void UKF::ComputeNewEstimate() {
    //predicted state mean
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x_ += weights_(i) * Xsig_pred_.col(i);
    }
    
    //predicted state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        x_diff(3) = fmod(x_diff(3), M_PI);
        
        P_ += weights_(i) * x_diff * x_diff.transpose() ;
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    //Augmented sigma points
    MatrixXd Xsig_aug = AugSigmaPoints();
    
    //Sigma points predictions
    PredictSigmaPoints(Xsig_aug, delta_t);
    
    //Predict mean and covariance matrix
    ComputeNewEstimate();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Use lidar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the lidar NIS.
    */
    
    VectorXd y = meas_package.raw_measurements_ - H_laser_ * x_;
    MatrixXd PHt = P_ * H_laser_.transpose();
    MatrixXd S = H_laser_ * PHt + R_laser_;
    Eigen::MatrixXd Si = S.inverse();
    MatrixXd K = PHt * Si;
    
    //new state
    x_ = x_ + (K * y);
    long x_size = x_.size();
    P_ = (MatrixXd::Identity(x_size, x_size) - K * H_laser_) * P_;
    
    nis_laser_ = y.transpose() * Si * y;
}

/**
* Updates the state and the state covariance matrix using a radar measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    
    int n_z = 3; //set measurement dimension, radar can measure r, phi, and r_dot
    
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);  // Matrix for sigma points in measurement space
    VectorXd z_pred = VectorXd(n_z);                // Mean predicted measurement
    MatrixXd S = MatrixXd(n_z, n_z);                // Measurement covariance matrix S
    MatrixXd Tc = MatrixXd(n_x_, n_z);              // Matrix for cross correlation Tc
    
    z_pred.fill(0.0);
    S.fill(0.0);
    Tc.fill(0.0);
    
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        
        const Eigen::VectorXd sigma_point = Xsig_pred_.col(i);
        const Eigen::VectorXd polar_point = tools.CartesianToPolarCTRV(sigma_point);
        
        Zsig(0,i) = polar_point(0); //rho
        Zsig(1,i) = polar_point(1); //phi
        Zsig(2,i) = polar_point(2); //r_dot
        
        z_pred += weights_(i) * Zsig.col(i);
    }
    
    //calculate cross correlation matrix
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        
        VectorXd z_diff = Zsig.col(i) - z_pred;         //residual
        S += weights_(i) * z_diff * z_diff.transpose();
        
        VectorXd x_diff = Xsig_pred_.col(i) - x_;       // state difference
        Tc += weights_(i) * x_diff * z_diff.transpose();
    }
    
    S += R_radar_;
    Eigen::MatrixXd Si = S.inverse();
    MatrixXd K = Tc * Si;                                       //Kalman gain K;
    VectorXd z_diff = meas_package.raw_measurements_ - z_pred;  //residual
    z_diff(1) = fmod(z_diff(1), M_PI);                          //angle normalization
    
    //update state mean and covariance matrix
    x_ += K * z_diff;
    P_ -= K * S * K.transpose();
    
    nis_radar_ = z_diff.transpose() * Si * z_diff;
}
