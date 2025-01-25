#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8551265930993792021) {
   out_8551265930993792021[0] = delta_x[0] + nom_x[0];
   out_8551265930993792021[1] = delta_x[1] + nom_x[1];
   out_8551265930993792021[2] = delta_x[2] + nom_x[2];
   out_8551265930993792021[3] = delta_x[3] + nom_x[3];
   out_8551265930993792021[4] = delta_x[4] + nom_x[4];
   out_8551265930993792021[5] = delta_x[5] + nom_x[5];
   out_8551265930993792021[6] = delta_x[6] + nom_x[6];
   out_8551265930993792021[7] = delta_x[7] + nom_x[7];
   out_8551265930993792021[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_5164373455773049213) {
   out_5164373455773049213[0] = -nom_x[0] + true_x[0];
   out_5164373455773049213[1] = -nom_x[1] + true_x[1];
   out_5164373455773049213[2] = -nom_x[2] + true_x[2];
   out_5164373455773049213[3] = -nom_x[3] + true_x[3];
   out_5164373455773049213[4] = -nom_x[4] + true_x[4];
   out_5164373455773049213[5] = -nom_x[5] + true_x[5];
   out_5164373455773049213[6] = -nom_x[6] + true_x[6];
   out_5164373455773049213[7] = -nom_x[7] + true_x[7];
   out_5164373455773049213[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_6881006698762961884) {
   out_6881006698762961884[0] = 1.0;
   out_6881006698762961884[1] = 0;
   out_6881006698762961884[2] = 0;
   out_6881006698762961884[3] = 0;
   out_6881006698762961884[4] = 0;
   out_6881006698762961884[5] = 0;
   out_6881006698762961884[6] = 0;
   out_6881006698762961884[7] = 0;
   out_6881006698762961884[8] = 0;
   out_6881006698762961884[9] = 0;
   out_6881006698762961884[10] = 1.0;
   out_6881006698762961884[11] = 0;
   out_6881006698762961884[12] = 0;
   out_6881006698762961884[13] = 0;
   out_6881006698762961884[14] = 0;
   out_6881006698762961884[15] = 0;
   out_6881006698762961884[16] = 0;
   out_6881006698762961884[17] = 0;
   out_6881006698762961884[18] = 0;
   out_6881006698762961884[19] = 0;
   out_6881006698762961884[20] = 1.0;
   out_6881006698762961884[21] = 0;
   out_6881006698762961884[22] = 0;
   out_6881006698762961884[23] = 0;
   out_6881006698762961884[24] = 0;
   out_6881006698762961884[25] = 0;
   out_6881006698762961884[26] = 0;
   out_6881006698762961884[27] = 0;
   out_6881006698762961884[28] = 0;
   out_6881006698762961884[29] = 0;
   out_6881006698762961884[30] = 1.0;
   out_6881006698762961884[31] = 0;
   out_6881006698762961884[32] = 0;
   out_6881006698762961884[33] = 0;
   out_6881006698762961884[34] = 0;
   out_6881006698762961884[35] = 0;
   out_6881006698762961884[36] = 0;
   out_6881006698762961884[37] = 0;
   out_6881006698762961884[38] = 0;
   out_6881006698762961884[39] = 0;
   out_6881006698762961884[40] = 1.0;
   out_6881006698762961884[41] = 0;
   out_6881006698762961884[42] = 0;
   out_6881006698762961884[43] = 0;
   out_6881006698762961884[44] = 0;
   out_6881006698762961884[45] = 0;
   out_6881006698762961884[46] = 0;
   out_6881006698762961884[47] = 0;
   out_6881006698762961884[48] = 0;
   out_6881006698762961884[49] = 0;
   out_6881006698762961884[50] = 1.0;
   out_6881006698762961884[51] = 0;
   out_6881006698762961884[52] = 0;
   out_6881006698762961884[53] = 0;
   out_6881006698762961884[54] = 0;
   out_6881006698762961884[55] = 0;
   out_6881006698762961884[56] = 0;
   out_6881006698762961884[57] = 0;
   out_6881006698762961884[58] = 0;
   out_6881006698762961884[59] = 0;
   out_6881006698762961884[60] = 1.0;
   out_6881006698762961884[61] = 0;
   out_6881006698762961884[62] = 0;
   out_6881006698762961884[63] = 0;
   out_6881006698762961884[64] = 0;
   out_6881006698762961884[65] = 0;
   out_6881006698762961884[66] = 0;
   out_6881006698762961884[67] = 0;
   out_6881006698762961884[68] = 0;
   out_6881006698762961884[69] = 0;
   out_6881006698762961884[70] = 1.0;
   out_6881006698762961884[71] = 0;
   out_6881006698762961884[72] = 0;
   out_6881006698762961884[73] = 0;
   out_6881006698762961884[74] = 0;
   out_6881006698762961884[75] = 0;
   out_6881006698762961884[76] = 0;
   out_6881006698762961884[77] = 0;
   out_6881006698762961884[78] = 0;
   out_6881006698762961884[79] = 0;
   out_6881006698762961884[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_8735568865516118734) {
   out_8735568865516118734[0] = state[0];
   out_8735568865516118734[1] = state[1];
   out_8735568865516118734[2] = state[2];
   out_8735568865516118734[3] = state[3];
   out_8735568865516118734[4] = state[4];
   out_8735568865516118734[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_8735568865516118734[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_8735568865516118734[7] = state[7];
   out_8735568865516118734[8] = state[8];
}
void F_fun(double *state, double dt, double *out_579048534215247874) {
   out_579048534215247874[0] = 1;
   out_579048534215247874[1] = 0;
   out_579048534215247874[2] = 0;
   out_579048534215247874[3] = 0;
   out_579048534215247874[4] = 0;
   out_579048534215247874[5] = 0;
   out_579048534215247874[6] = 0;
   out_579048534215247874[7] = 0;
   out_579048534215247874[8] = 0;
   out_579048534215247874[9] = 0;
   out_579048534215247874[10] = 1;
   out_579048534215247874[11] = 0;
   out_579048534215247874[12] = 0;
   out_579048534215247874[13] = 0;
   out_579048534215247874[14] = 0;
   out_579048534215247874[15] = 0;
   out_579048534215247874[16] = 0;
   out_579048534215247874[17] = 0;
   out_579048534215247874[18] = 0;
   out_579048534215247874[19] = 0;
   out_579048534215247874[20] = 1;
   out_579048534215247874[21] = 0;
   out_579048534215247874[22] = 0;
   out_579048534215247874[23] = 0;
   out_579048534215247874[24] = 0;
   out_579048534215247874[25] = 0;
   out_579048534215247874[26] = 0;
   out_579048534215247874[27] = 0;
   out_579048534215247874[28] = 0;
   out_579048534215247874[29] = 0;
   out_579048534215247874[30] = 1;
   out_579048534215247874[31] = 0;
   out_579048534215247874[32] = 0;
   out_579048534215247874[33] = 0;
   out_579048534215247874[34] = 0;
   out_579048534215247874[35] = 0;
   out_579048534215247874[36] = 0;
   out_579048534215247874[37] = 0;
   out_579048534215247874[38] = 0;
   out_579048534215247874[39] = 0;
   out_579048534215247874[40] = 1;
   out_579048534215247874[41] = 0;
   out_579048534215247874[42] = 0;
   out_579048534215247874[43] = 0;
   out_579048534215247874[44] = 0;
   out_579048534215247874[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_579048534215247874[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_579048534215247874[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_579048534215247874[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_579048534215247874[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_579048534215247874[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_579048534215247874[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_579048534215247874[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_579048534215247874[53] = -9.8000000000000007*dt;
   out_579048534215247874[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_579048534215247874[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_579048534215247874[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_579048534215247874[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_579048534215247874[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_579048534215247874[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_579048534215247874[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_579048534215247874[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_579048534215247874[62] = 0;
   out_579048534215247874[63] = 0;
   out_579048534215247874[64] = 0;
   out_579048534215247874[65] = 0;
   out_579048534215247874[66] = 0;
   out_579048534215247874[67] = 0;
   out_579048534215247874[68] = 0;
   out_579048534215247874[69] = 0;
   out_579048534215247874[70] = 1;
   out_579048534215247874[71] = 0;
   out_579048534215247874[72] = 0;
   out_579048534215247874[73] = 0;
   out_579048534215247874[74] = 0;
   out_579048534215247874[75] = 0;
   out_579048534215247874[76] = 0;
   out_579048534215247874[77] = 0;
   out_579048534215247874[78] = 0;
   out_579048534215247874[79] = 0;
   out_579048534215247874[80] = 1;
}
void h_25(double *state, double *unused, double *out_6441091650644349478) {
   out_6441091650644349478[0] = state[6];
}
void H_25(double *state, double *unused, double *out_6867098003517625916) {
   out_6867098003517625916[0] = 0;
   out_6867098003517625916[1] = 0;
   out_6867098003517625916[2] = 0;
   out_6867098003517625916[3] = 0;
   out_6867098003517625916[4] = 0;
   out_6867098003517625916[5] = 0;
   out_6867098003517625916[6] = 1;
   out_6867098003517625916[7] = 0;
   out_6867098003517625916[8] = 0;
}
void h_24(double *state, double *unused, double *out_8791765496533727468) {
   out_8791765496533727468[0] = state[4];
   out_8791765496533727468[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4694448404512126350) {
   out_4694448404512126350[0] = 0;
   out_4694448404512126350[1] = 0;
   out_4694448404512126350[2] = 0;
   out_4694448404512126350[3] = 0;
   out_4694448404512126350[4] = 1;
   out_4694448404512126350[5] = 0;
   out_4694448404512126350[6] = 0;
   out_4694448404512126350[7] = 0;
   out_4694448404512126350[8] = 0;
   out_4694448404512126350[9] = 0;
   out_4694448404512126350[10] = 0;
   out_4694448404512126350[11] = 0;
   out_4694448404512126350[12] = 0;
   out_4694448404512126350[13] = 0;
   out_4694448404512126350[14] = 1;
   out_4694448404512126350[15] = 0;
   out_4694448404512126350[16] = 0;
   out_4694448404512126350[17] = 0;
}
void h_30(double *state, double *unused, double *out_3221843324619681176) {
   out_3221843324619681176[0] = state[4];
}
void H_30(double *state, double *unused, double *out_4662955728700308945) {
   out_4662955728700308945[0] = 0;
   out_4662955728700308945[1] = 0;
   out_4662955728700308945[2] = 0;
   out_4662955728700308945[3] = 0;
   out_4662955728700308945[4] = 1;
   out_4662955728700308945[5] = 0;
   out_4662955728700308945[6] = 0;
   out_4662955728700308945[7] = 0;
   out_4662955728700308945[8] = 0;
}
void h_26(double *state, double *unused, double *out_1635898321314185391) {
   out_1635898321314185391[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3125594684643569692) {
   out_3125594684643569692[0] = 0;
   out_3125594684643569692[1] = 0;
   out_3125594684643569692[2] = 0;
   out_3125594684643569692[3] = 0;
   out_3125594684643569692[4] = 0;
   out_3125594684643569692[5] = 0;
   out_3125594684643569692[6] = 0;
   out_3125594684643569692[7] = 1;
   out_3125594684643569692[8] = 0;
}
void h_27(double *state, double *unused, double *out_5459007557309051366) {
   out_5459007557309051366[0] = state[3];
}
void H_27(double *state, double *unused, double *out_6837719040500733856) {
   out_6837719040500733856[0] = 0;
   out_6837719040500733856[1] = 0;
   out_6837719040500733856[2] = 0;
   out_6837719040500733856[3] = 1;
   out_6837719040500733856[4] = 0;
   out_6837719040500733856[5] = 0;
   out_6837719040500733856[6] = 0;
   out_6837719040500733856[7] = 0;
   out_6837719040500733856[8] = 0;
}
void h_29(double *state, double *unused, double *out_5183813495024545477) {
   out_5183813495024545477[0] = state[1];
}
void H_29(double *state, double *unused, double *out_4152724384385916761) {
   out_4152724384385916761[0] = 0;
   out_4152724384385916761[1] = 1;
   out_4152724384385916761[2] = 0;
   out_4152724384385916761[3] = 0;
   out_4152724384385916761[4] = 0;
   out_4152724384385916761[5] = 0;
   out_4152724384385916761[6] = 0;
   out_4152724384385916761[7] = 0;
   out_4152724384385916761[8] = 0;
}
void h_28(double *state, double *unused, double *out_2190337616002282976) {
   out_2190337616002282976[0] = state[0];
}
void H_28(double *state, double *unused, double *out_4813263289269736153) {
   out_4813263289269736153[0] = 1;
   out_4813263289269736153[1] = 0;
   out_4813263289269736153[2] = 0;
   out_4813263289269736153[3] = 0;
   out_4813263289269736153[4] = 0;
   out_4813263289269736153[5] = 0;
   out_4813263289269736153[6] = 0;
   out_4813263289269736153[7] = 0;
   out_4813263289269736153[8] = 0;
}
void h_31(double *state, double *unused, double *out_4348613720485632255) {
   out_4348613720485632255[0] = state[8];
}
void H_31(double *state, double *unused, double *out_6897743965394586344) {
   out_6897743965394586344[0] = 0;
   out_6897743965394586344[1] = 0;
   out_6897743965394586344[2] = 0;
   out_6897743965394586344[3] = 0;
   out_6897743965394586344[4] = 0;
   out_6897743965394586344[5] = 0;
   out_6897743965394586344[6] = 0;
   out_6897743965394586344[7] = 0;
   out_6897743965394586344[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_8551265930993792021) {
  err_fun(nom_x, delta_x, out_8551265930993792021);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5164373455773049213) {
  inv_err_fun(nom_x, true_x, out_5164373455773049213);
}
void car_H_mod_fun(double *state, double *out_6881006698762961884) {
  H_mod_fun(state, out_6881006698762961884);
}
void car_f_fun(double *state, double dt, double *out_8735568865516118734) {
  f_fun(state,  dt, out_8735568865516118734);
}
void car_F_fun(double *state, double dt, double *out_579048534215247874) {
  F_fun(state,  dt, out_579048534215247874);
}
void car_h_25(double *state, double *unused, double *out_6441091650644349478) {
  h_25(state, unused, out_6441091650644349478);
}
void car_H_25(double *state, double *unused, double *out_6867098003517625916) {
  H_25(state, unused, out_6867098003517625916);
}
void car_h_24(double *state, double *unused, double *out_8791765496533727468) {
  h_24(state, unused, out_8791765496533727468);
}
void car_H_24(double *state, double *unused, double *out_4694448404512126350) {
  H_24(state, unused, out_4694448404512126350);
}
void car_h_30(double *state, double *unused, double *out_3221843324619681176) {
  h_30(state, unused, out_3221843324619681176);
}
void car_H_30(double *state, double *unused, double *out_4662955728700308945) {
  H_30(state, unused, out_4662955728700308945);
}
void car_h_26(double *state, double *unused, double *out_1635898321314185391) {
  h_26(state, unused, out_1635898321314185391);
}
void car_H_26(double *state, double *unused, double *out_3125594684643569692) {
  H_26(state, unused, out_3125594684643569692);
}
void car_h_27(double *state, double *unused, double *out_5459007557309051366) {
  h_27(state, unused, out_5459007557309051366);
}
void car_H_27(double *state, double *unused, double *out_6837719040500733856) {
  H_27(state, unused, out_6837719040500733856);
}
void car_h_29(double *state, double *unused, double *out_5183813495024545477) {
  h_29(state, unused, out_5183813495024545477);
}
void car_H_29(double *state, double *unused, double *out_4152724384385916761) {
  H_29(state, unused, out_4152724384385916761);
}
void car_h_28(double *state, double *unused, double *out_2190337616002282976) {
  h_28(state, unused, out_2190337616002282976);
}
void car_H_28(double *state, double *unused, double *out_4813263289269736153) {
  H_28(state, unused, out_4813263289269736153);
}
void car_h_31(double *state, double *unused, double *out_4348613720485632255) {
  h_31(state, unused, out_4348613720485632255);
}
void car_H_31(double *state, double *unused, double *out_6897743965394586344) {
  H_31(state, unused, out_6897743965394586344);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_lib_init(car)
