#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
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
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  is_initialized_ = false;

  P_ = MatrixXd::Identity(5,5);

  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  R_radar_ = MatrixXd(3,3);
  R_radar_ <<  std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  P_ = MatrixXd::Identity(5,5);
  R_laser_ = MatrixXd(2,2);
  R_laser_ <<  std_laspx_*std_laspx_, 0,
          0,std_laspy_*std_laspy_;


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (!is_initialized_){
    time_us_ = meas_package.timestamp_;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
	  float px = rho * cos(phi);
	  float py = rho * sin(phi);
      x_ << px, py, 0, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      x_ << meas_package.raw_measurements_[0],
            meas_package.raw_measurements_[1],
            0, 0, 0;
      // use smaller variance for lidar
      P_(0,0) = 0.3;
      P_(1,1) = 0.3;
    }
    is_initialized_ = true;
    return;
  }

  //compute the time elapsed between the current and previous measurements
  double dt = (double)(meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;
  while(dt > 0.1){
    const double dt_max = 0.05;
    Prediction(dt_max);
    dt -= dt_max;
  }
  Prediction(dt);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
    UpdateLidar(meas_package);
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */


  VectorXd x_aug = VectorXd(7);
  MatrixXd P_aug = MatrixXd(7, 7);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;


  MatrixXd L = P_aug.llt().matrixL();
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  for (int i = 0; i< 2*n_aug_+1; i++)
  {

    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);



    double px_p, py_p;
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin (yaw_p) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw_p) );
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    //add noise
    double dt2 = delta_t * delta_t;
    px_p += 0.5*nu_a*dt2 * cos(yaw);
    py_p += 0.5*nu_a*dt2 * sin(yaw);
    v_p  += nu_a*delta_t;

    yaw_p  += 0.5*nu_yawdd*dt2;
    yawd_p += nu_yawdd*delta_t;
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  
  VectorXd x_pred_ = Xsig_pred_ * weights_ ;
  MatrixXd P_pred_ = MatrixXd(n_x_, n_x_);
  P_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

  	VectorXd x_diff = Xsig_pred_.col(i) - x_pred_;
   // normalization
  	while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
  	while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

   P_pred_ = P_pred_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

  x_ = x_pred_;
  P_ = P_pred_;
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

  int n_z = 2; 
  VectorXd z = meas_package.raw_measurements_; 
  MatrixXd Z_sig = Xsig_pred_.block(0,0, n_z, 2 * n_aug_ + 1);

  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Z_sig.col(i);
  }

  // covariance
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  
    //residual
    VectorXd z_diff = Z_sig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd Tc = MatrixXd(5,2);
  Tc.fill(0.0);
  S = S + R_laser_;
  for (int i = 0; i < 2 * n_aug_ + 1; i++) { 
    VectorXd z_diff = Z_sig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman
  MatrixXd S_inv = S.inverse();
  MatrixXd K = Tc * S_inv;

  //residual
  VectorXd z_diff = z - z_pred;

  //update step
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  NIS_laser_ = z_diff.transpose() * S_inv * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  int n_z = 3; //dim of measurement vector
  VectorXd z = meas_package.raw_measurements_; // measurement vector
  // matrix for sigma points in measurement space
  MatrixXd Z_sig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;
    double radius = sqrt(p_x*p_x + p_y*p_y);
    double help_1 = cos(yaw)*v;
    double help_2 = sin(yaw)*v;
    double radius_ = sqrt(p_x*p_x + p_y*p_y);
    double small_num = 0.0001;
    double sq2 = 1.4142136;
    if (radius_ < small_num*sq2){
      p_x=small_num;
      p_y=small_num;
      radius_=small_num*sq2;
    }
    Z_sig(0,i) = radius;                         
    Z_sig(1,i) = atan2(p_y,p_x);            
    Z_sig(2,i) = (p_x*v1 + p_y*v2)/radius;       
  }

  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Z_sig.col(i);
  }

  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  
    //residual
    VectorXd z_diff = Z_sig.col(i) - z_pred;

    // normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S = S + R_radar_;

  MatrixXd T_c = MatrixXd(n_x_, n_z);
  T_c.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  
    VectorXd z_diff = Z_sig.col(i) - z_pred;
    // normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    T_c = T_c + weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd S_inv = S.inverse();
  MatrixXd K = T_c * S_inv;

  //residual
  VectorXd z_diff = z - z_pred;

  // normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  NIS_radar_ = z_diff.transpose() * S_inv * z_diff;
}
