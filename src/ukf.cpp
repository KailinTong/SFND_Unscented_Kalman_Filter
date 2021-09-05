#include <iostream>
#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */



    // State dimension
    n_x_ = 5;

    // Augmented state dimension
    n_aug_ = 7;

    // Sigma point spreading parameter
    lambda_ = 3 - n_aug_;

    // initialized
    is_initialized_ = false;

    // Weights of sigma points
    weights_ = VectorXd(2*n_aug_+1);
    double weight_0 = lambda_/(lambda_+n_aug_);
    double weight = 0.5/(lambda_+n_aug_);
    weights_.fill(weight);
    weights_(0) = weight_0;

    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
   if(is_initialized_){
       if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
           double dt = (meas_package.timestamp_ - time_us_)/1e6;
           time_us_ = meas_package.timestamp_;
           Prediction(dt);
           UpdateRadar(meas_package);
//           std::cout << x_ << std::endl;
//           std::cout << P_ << std::endl;

       }

       if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
           double dt = (meas_package.timestamp_ - time_us_)/1e6;
           time_us_ = meas_package.timestamp_;
           Prediction(dt);
           UpdateLidar(meas_package);
//           std::cout << x_ << std::endl;
//           std::cout << P_ << std::endl;
       }
   }
   else{
       if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
           P_ << std_laspx_ * std_laspx_, 0, 0, 0, 0,
                   0, std_laspy_ * std_laspy_, 0, 0, 0,
                   0, 0, 1, 0, 0,
                   0, 0, 0, 1, 0,
                   0, 0, 0, 0, 1;

           double radr =  meas_package.raw_measurements_(0);
           double radphi =  meas_package.raw_measurements_(1);
           double radrd =  meas_package.raw_measurements_(2);

           double px = radr * cos(radphi);
           double py = radr * sin(radphi);
           x_ << px, py, 0, 0, 0;

           time_us_ = meas_package.timestamp_; // in us

           is_initialized_ = true;
       }

       if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
           P_ << std_radr_* std_radr_, 0, 0, 0, 0,
                   0, std_radr_ * std_radr_, 0, 0, 0,
                   0, 0, std_radrd_ * std_radrd_, 0, 0,
                   0, 0, 0, std_radphi_ * std_radphi_, 0,
                   0, 0, 0, 0, std_radphi_ * std_radphi_;
           // following suggestions from https://knowledge.udacity.com/questions/654359

           double px = meas_package.raw_measurements_(0);
           double py = meas_package.raw_measurements_(1);
           x_ << px, py, 0, 0, 0;


           time_us_ = meas_package.timestamp_; // in us

           is_initialized_ = true;
       }
   }


}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */


    // create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);

    // create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

    // create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);


    // create augmented mean state
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    // create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5,5) = P_;
    P_aug(5,5) = std_a_*std_a_;
    P_aug(6,6) = std_yawdd_*std_yawdd_;

    // create square root matrix
    MatrixXd L = P_aug.llt().matrixL();

    // create augmented sigma points
    Xsig_aug.col(0)  = x_aug;
    for (int i = 0; i< n_aug_; ++i) {
        Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
    }


    // predict sigma points
    for (int i = 0; i< 2*n_aug_+1; ++i) {
        // extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);

        // predicted state values
        double px_p, py_p;

        // avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        } else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }

        // add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);

        // write predicted sigma point into right column
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;

    }

    // v
    Xsig_pred_.row(2) = Xsig_aug.row(2) + delta_t * Xsig_aug.row(5);
    // yaw
    Xsig_pred_.row(3) = Xsig_aug.row(3) + delta_t * Xsig_aug.row(4) + 0.5 * delta_t * delta_t * Xsig_aug.row(6);
    // yaw rate
    Xsig_pred_.row(4) = Xsig_aug.row(4) + delta_t * Xsig_aug.row(6);

    // predict state mean
    MatrixXd weighted_sum = (Xsig_pred_ * weights_);
    x_ = weighted_sum.rowwise().sum();
    // predict state covariance matrix

    P_.fill(0.0);

    for(int i = 0; i < 2 * n_aug_ + 1; i++){
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        // angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
    }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

    int n_z = 2;

    // create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    // mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);

    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);


    // transform sigma points into measurement space
    for(int i = 0; i < 2 * n_aug_ + 1; i++){
        double px =  Xsig_pred_(0, i);
        double py =  Xsig_pred_(1, i);
        double v =  Xsig_pred_(2, i);
        double yaw =  Xsig_pred_(3, i);
        double yawd =  Xsig_pred_(4, i);


        Zsig(0, i) = px;
        Zsig(1, i) = py;

    }
    // calculate mean predicted measurement
    z_pred = (Zsig * weights_).rowwise().mean();

    // calculate innovation covariance matrix S
    S.fill(0.0);
    for(int i = 0; i < n_aug_ * 2 + 1; i++){
        VectorXd z_diff =  Zsig.col(i) - z_pred;
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    MatrixXd R(n_z, n_z);
    R<<std_laspx_*std_laspx_, 0,
            0, std_laspy_*std_laspy_;
    S = S + R;


    // create z vector for incoming radar measurement
    VectorXd z = VectorXd(n_z);
    z <<    meas_package.raw_measurements_(0),   // px in m
            meas_package.raw_measurements_(1);   // py in m

    // create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    // calculate cross correlation matrix
    Tc.fill(0.0);
    for(int i = 0; i < 2 * n_aug_ + 1; i++){
        Tc = Tc + weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose();
    }
    // calculate Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    // update state mean and covariance matrix
    x_ = x_ + K * (z - z_pred);
    P_ = P_ - K * S * K.transpose();

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

    int n_z = 3;

    // create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    // mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);

    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);


    // transform sigma points into measurement space
    for(int i = 0; i < 2 * n_aug_ + 1; i++){
        double px =  Xsig_pred_(0, i);
        double py =  Xsig_pred_(1, i);
        double v =  Xsig_pred_(2, i);
        double yaw =  Xsig_pred_(3, i);
        double yawd =  Xsig_pred_(4, i);

        double radr = sqrt(px*px + py*py);
        double radphi = atan2(py, px);
        double radrd = (px*cos(yaw) + py*sin(yaw))*v/sqrt(px*px + py*py);

        Zsig(0, i) = radr;
        Zsig(1, i) = radphi;
        Zsig(2, i) = radrd;

    }
    // calculate mean predicted measurement
    z_pred = (Zsig * weights_).rowwise().mean();

    // calculate innovation covariance matrix S
    S.fill(0.0);
    for(int i = 0; i < n_aug_ * 2 + 1; i++){
        VectorXd z_diff =  Zsig.col(i) - z_pred;
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    MatrixXd R(n_z, n_z);
    R<<std_radr_*std_radr_, 0, 0,
            0, std_radphi_*std_radphi_, 0,
            0, 0, std_radrd_*std_radrd_;
    S = S + R;


    // create z vector for incoming radar measurement
    VectorXd z = VectorXd(n_z);
    double radr_meas =  meas_package.raw_measurements_(0);
    double radphi_meas =  meas_package.raw_measurements_(1);
    double radrd_meas =  meas_package.raw_measurements_(2);
    z <<    radr_meas,   // rho in m
            radphi_meas,   // phi in rad
            radrd_meas;   // rho_dot in m/s

    // create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    // calculate cross correlation matrix
    Tc.fill(0.0);
    for(int i = 0; i < 2 * n_aug_ + 1; i++){
        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        // angle normalization important!
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        // angle normalization important!
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    // calculate Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    // update state mean and covariance matrix
    x_ = x_ + K * (z - z_pred);
    P_ = P_ - K * S * K.transpose();

}