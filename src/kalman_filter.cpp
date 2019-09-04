#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd y;
  MatrixXd S;
  MatrixXd K;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  
  y = z - H_ * x_;
  S = H_ * P_ * H_.transpose() + R_;
  K = P_ * H_.transpose() * S.inverse();

  // new state
  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  VectorXd y;
  MatrixXd S;
  MatrixXd K;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  // on this condition H_ is Hj
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);

  VectorXd h_x(3);
  h_x << sqrt(px*px + py*py),
         atan2(py, px),
         (px*vx + py*vy) / sqrt(px*px + py* py);
         
  y = z - h_x;
  
  while (y(1)> M_PI) y(1)-=2.*M_PI;
  while (y(1)<-M_PI) y(1)+=2.*M_PI;

  S = H_ * P_ * H_.transpose() + R_;
  K = P_ * H_.transpose() * S.inverse();

  // new state
  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
}
