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
void err_fun(double *nom_x, double *delta_x, double *out_5464861526568687019) {
   out_5464861526568687019[0] = delta_x[0] + nom_x[0];
   out_5464861526568687019[1] = delta_x[1] + nom_x[1];
   out_5464861526568687019[2] = delta_x[2] + nom_x[2];
   out_5464861526568687019[3] = delta_x[3] + nom_x[3];
   out_5464861526568687019[4] = delta_x[4] + nom_x[4];
   out_5464861526568687019[5] = delta_x[5] + nom_x[5];
   out_5464861526568687019[6] = delta_x[6] + nom_x[6];
   out_5464861526568687019[7] = delta_x[7] + nom_x[7];
   out_5464861526568687019[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_7516906840074992970) {
   out_7516906840074992970[0] = -nom_x[0] + true_x[0];
   out_7516906840074992970[1] = -nom_x[1] + true_x[1];
   out_7516906840074992970[2] = -nom_x[2] + true_x[2];
   out_7516906840074992970[3] = -nom_x[3] + true_x[3];
   out_7516906840074992970[4] = -nom_x[4] + true_x[4];
   out_7516906840074992970[5] = -nom_x[5] + true_x[5];
   out_7516906840074992970[6] = -nom_x[6] + true_x[6];
   out_7516906840074992970[7] = -nom_x[7] + true_x[7];
   out_7516906840074992970[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_6431729039895127708) {
   out_6431729039895127708[0] = 1.0;
   out_6431729039895127708[1] = 0;
   out_6431729039895127708[2] = 0;
   out_6431729039895127708[3] = 0;
   out_6431729039895127708[4] = 0;
   out_6431729039895127708[5] = 0;
   out_6431729039895127708[6] = 0;
   out_6431729039895127708[7] = 0;
   out_6431729039895127708[8] = 0;
   out_6431729039895127708[9] = 0;
   out_6431729039895127708[10] = 1.0;
   out_6431729039895127708[11] = 0;
   out_6431729039895127708[12] = 0;
   out_6431729039895127708[13] = 0;
   out_6431729039895127708[14] = 0;
   out_6431729039895127708[15] = 0;
   out_6431729039895127708[16] = 0;
   out_6431729039895127708[17] = 0;
   out_6431729039895127708[18] = 0;
   out_6431729039895127708[19] = 0;
   out_6431729039895127708[20] = 1.0;
   out_6431729039895127708[21] = 0;
   out_6431729039895127708[22] = 0;
   out_6431729039895127708[23] = 0;
   out_6431729039895127708[24] = 0;
   out_6431729039895127708[25] = 0;
   out_6431729039895127708[26] = 0;
   out_6431729039895127708[27] = 0;
   out_6431729039895127708[28] = 0;
   out_6431729039895127708[29] = 0;
   out_6431729039895127708[30] = 1.0;
   out_6431729039895127708[31] = 0;
   out_6431729039895127708[32] = 0;
   out_6431729039895127708[33] = 0;
   out_6431729039895127708[34] = 0;
   out_6431729039895127708[35] = 0;
   out_6431729039895127708[36] = 0;
   out_6431729039895127708[37] = 0;
   out_6431729039895127708[38] = 0;
   out_6431729039895127708[39] = 0;
   out_6431729039895127708[40] = 1.0;
   out_6431729039895127708[41] = 0;
   out_6431729039895127708[42] = 0;
   out_6431729039895127708[43] = 0;
   out_6431729039895127708[44] = 0;
   out_6431729039895127708[45] = 0;
   out_6431729039895127708[46] = 0;
   out_6431729039895127708[47] = 0;
   out_6431729039895127708[48] = 0;
   out_6431729039895127708[49] = 0;
   out_6431729039895127708[50] = 1.0;
   out_6431729039895127708[51] = 0;
   out_6431729039895127708[52] = 0;
   out_6431729039895127708[53] = 0;
   out_6431729039895127708[54] = 0;
   out_6431729039895127708[55] = 0;
   out_6431729039895127708[56] = 0;
   out_6431729039895127708[57] = 0;
   out_6431729039895127708[58] = 0;
   out_6431729039895127708[59] = 0;
   out_6431729039895127708[60] = 1.0;
   out_6431729039895127708[61] = 0;
   out_6431729039895127708[62] = 0;
   out_6431729039895127708[63] = 0;
   out_6431729039895127708[64] = 0;
   out_6431729039895127708[65] = 0;
   out_6431729039895127708[66] = 0;
   out_6431729039895127708[67] = 0;
   out_6431729039895127708[68] = 0;
   out_6431729039895127708[69] = 0;
   out_6431729039895127708[70] = 1.0;
   out_6431729039895127708[71] = 0;
   out_6431729039895127708[72] = 0;
   out_6431729039895127708[73] = 0;
   out_6431729039895127708[74] = 0;
   out_6431729039895127708[75] = 0;
   out_6431729039895127708[76] = 0;
   out_6431729039895127708[77] = 0;
   out_6431729039895127708[78] = 0;
   out_6431729039895127708[79] = 0;
   out_6431729039895127708[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_57451106597728332) {
   out_57451106597728332[0] = state[0];
   out_57451106597728332[1] = state[1];
   out_57451106597728332[2] = state[2];
   out_57451106597728332[3] = state[3];
   out_57451106597728332[4] = state[4];
   out_57451106597728332[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_57451106597728332[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_57451106597728332[7] = state[7];
   out_57451106597728332[8] = state[8];
}
void F_fun(double *state, double dt, double *out_447865623327851108) {
   out_447865623327851108[0] = 1;
   out_447865623327851108[1] = 0;
   out_447865623327851108[2] = 0;
   out_447865623327851108[3] = 0;
   out_447865623327851108[4] = 0;
   out_447865623327851108[5] = 0;
   out_447865623327851108[6] = 0;
   out_447865623327851108[7] = 0;
   out_447865623327851108[8] = 0;
   out_447865623327851108[9] = 0;
   out_447865623327851108[10] = 1;
   out_447865623327851108[11] = 0;
   out_447865623327851108[12] = 0;
   out_447865623327851108[13] = 0;
   out_447865623327851108[14] = 0;
   out_447865623327851108[15] = 0;
   out_447865623327851108[16] = 0;
   out_447865623327851108[17] = 0;
   out_447865623327851108[18] = 0;
   out_447865623327851108[19] = 0;
   out_447865623327851108[20] = 1;
   out_447865623327851108[21] = 0;
   out_447865623327851108[22] = 0;
   out_447865623327851108[23] = 0;
   out_447865623327851108[24] = 0;
   out_447865623327851108[25] = 0;
   out_447865623327851108[26] = 0;
   out_447865623327851108[27] = 0;
   out_447865623327851108[28] = 0;
   out_447865623327851108[29] = 0;
   out_447865623327851108[30] = 1;
   out_447865623327851108[31] = 0;
   out_447865623327851108[32] = 0;
   out_447865623327851108[33] = 0;
   out_447865623327851108[34] = 0;
   out_447865623327851108[35] = 0;
   out_447865623327851108[36] = 0;
   out_447865623327851108[37] = 0;
   out_447865623327851108[38] = 0;
   out_447865623327851108[39] = 0;
   out_447865623327851108[40] = 1;
   out_447865623327851108[41] = 0;
   out_447865623327851108[42] = 0;
   out_447865623327851108[43] = 0;
   out_447865623327851108[44] = 0;
   out_447865623327851108[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_447865623327851108[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_447865623327851108[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_447865623327851108[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_447865623327851108[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_447865623327851108[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_447865623327851108[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_447865623327851108[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_447865623327851108[53] = -9.8000000000000007*dt;
   out_447865623327851108[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_447865623327851108[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_447865623327851108[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_447865623327851108[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_447865623327851108[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_447865623327851108[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_447865623327851108[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_447865623327851108[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_447865623327851108[62] = 0;
   out_447865623327851108[63] = 0;
   out_447865623327851108[64] = 0;
   out_447865623327851108[65] = 0;
   out_447865623327851108[66] = 0;
   out_447865623327851108[67] = 0;
   out_447865623327851108[68] = 0;
   out_447865623327851108[69] = 0;
   out_447865623327851108[70] = 1;
   out_447865623327851108[71] = 0;
   out_447865623327851108[72] = 0;
   out_447865623327851108[73] = 0;
   out_447865623327851108[74] = 0;
   out_447865623327851108[75] = 0;
   out_447865623327851108[76] = 0;
   out_447865623327851108[77] = 0;
   out_447865623327851108[78] = 0;
   out_447865623327851108[79] = 0;
   out_447865623327851108[80] = 1;
}
void h_25(double *state, double *unused, double *out_4364545812939371539) {
   out_4364545812939371539[0] = state[6];
}
void H_25(double *state, double *unused, double *out_2450797190468748385) {
   out_2450797190468748385[0] = 0;
   out_2450797190468748385[1] = 0;
   out_2450797190468748385[2] = 0;
   out_2450797190468748385[3] = 0;
   out_2450797190468748385[4] = 0;
   out_2450797190468748385[5] = 0;
   out_2450797190468748385[6] = 1;
   out_2450797190468748385[7] = 0;
   out_2450797190468748385[8] = 0;
}
void h_24(double *state, double *unused, double *out_4446211979141650588) {
   out_4446211979141650588[0] = state[4];
   out_4446211979141650588[1] = state[5];
}
void H_24(double *state, double *unused, double *out_3434863061310062053) {
   out_3434863061310062053[0] = 0;
   out_3434863061310062053[1] = 0;
   out_3434863061310062053[2] = 0;
   out_3434863061310062053[3] = 0;
   out_3434863061310062053[4] = 1;
   out_3434863061310062053[5] = 0;
   out_3434863061310062053[6] = 0;
   out_3434863061310062053[7] = 0;
   out_3434863061310062053[8] = 0;
   out_3434863061310062053[9] = 0;
   out_3434863061310062053[10] = 0;
   out_3434863061310062053[11] = 0;
   out_3434863061310062053[12] = 0;
   out_3434863061310062053[13] = 0;
   out_3434863061310062053[14] = 1;
   out_3434863061310062053[15] = 0;
   out_3434863061310062053[16] = 0;
   out_3434863061310062053[17] = 0;
}
void h_30(double *state, double *unused, double *out_3516726837814677646) {
   out_3516726837814677646[0] = state[4];
}
void H_30(double *state, double *unused, double *out_9079256541749186476) {
   out_9079256541749186476[0] = 0;
   out_9079256541749186476[1] = 0;
   out_9079256541749186476[2] = 0;
   out_9079256541749186476[3] = 0;
   out_9079256541749186476[4] = 1;
   out_9079256541749186476[5] = 0;
   out_9079256541749186476[6] = 0;
   out_9079256541749186476[7] = 0;
   out_9079256541749186476[8] = 0;
}
void h_26(double *state, double *unused, double *out_6797261298903552818) {
   out_6797261298903552818[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5755323160229548986) {
   out_5755323160229548986[0] = 0;
   out_5755323160229548986[1] = 0;
   out_5755323160229548986[2] = 0;
   out_5755323160229548986[3] = 0;
   out_5755323160229548986[4] = 0;
   out_5755323160229548986[5] = 0;
   out_5755323160229548986[6] = 0;
   out_5755323160229548986[7] = 1;
   out_5755323160229548986[8] = 0;
}
void h_27(double *state, double *unused, double *out_9119818263088886845) {
   out_9119818263088886845[0] = state[3];
}
void H_27(double *state, double *unused, double *out_7192724220159940229) {
   out_7192724220159940229[0] = 0;
   out_7192724220159940229[1] = 0;
   out_7192724220159940229[2] = 0;
   out_7192724220159940229[3] = 1;
   out_7192724220159940229[4] = 0;
   out_7192724220159940229[5] = 0;
   out_7192724220159940229[6] = 0;
   out_7192724220159940229[7] = 0;
   out_7192724220159940229[8] = 0;
}
void h_29(double *state, double *unused, double *out_6922248888625985153) {
   out_6922248888625985153[0] = state[1];
}
void H_29(double *state, double *unused, double *out_8569025197434794292) {
   out_8569025197434794292[0] = 0;
   out_8569025197434794292[1] = 1;
   out_8569025197434794292[2] = 0;
   out_8569025197434794292[3] = 0;
   out_8569025197434794292[4] = 0;
   out_8569025197434794292[5] = 0;
   out_8569025197434794292[6] = 0;
   out_8569025197434794292[7] = 0;
   out_8569025197434794292[8] = 0;
}
void h_28(double *state, double *unused, double *out_5703544973439816829) {
   out_5703544973439816829[0] = state[0];
}
void H_28(double *state, double *unused, double *out_396962476220858622) {
   out_396962476220858622[0] = 1;
   out_396962476220858622[1] = 0;
   out_396962476220858622[2] = 0;
   out_396962476220858622[3] = 0;
   out_396962476220858622[4] = 0;
   out_396962476220858622[5] = 0;
   out_396962476220858622[6] = 0;
   out_396962476220858622[7] = 0;
   out_396962476220858622[8] = 0;
}
void h_31(double *state, double *unused, double *out_266242514659999675) {
   out_266242514659999675[0] = state[8];
}
void H_31(double *state, double *unused, double *out_8919271632728985978) {
   out_8919271632728985978[0] = 0;
   out_8919271632728985978[1] = 0;
   out_8919271632728985978[2] = 0;
   out_8919271632728985978[3] = 0;
   out_8919271632728985978[4] = 0;
   out_8919271632728985978[5] = 0;
   out_8919271632728985978[6] = 0;
   out_8919271632728985978[7] = 0;
   out_8919271632728985978[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_5464861526568687019) {
  err_fun(nom_x, delta_x, out_5464861526568687019);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_7516906840074992970) {
  inv_err_fun(nom_x, true_x, out_7516906840074992970);
}
void car_H_mod_fun(double *state, double *out_6431729039895127708) {
  H_mod_fun(state, out_6431729039895127708);
}
void car_f_fun(double *state, double dt, double *out_57451106597728332) {
  f_fun(state,  dt, out_57451106597728332);
}
void car_F_fun(double *state, double dt, double *out_447865623327851108) {
  F_fun(state,  dt, out_447865623327851108);
}
void car_h_25(double *state, double *unused, double *out_4364545812939371539) {
  h_25(state, unused, out_4364545812939371539);
}
void car_H_25(double *state, double *unused, double *out_2450797190468748385) {
  H_25(state, unused, out_2450797190468748385);
}
void car_h_24(double *state, double *unused, double *out_4446211979141650588) {
  h_24(state, unused, out_4446211979141650588);
}
void car_H_24(double *state, double *unused, double *out_3434863061310062053) {
  H_24(state, unused, out_3434863061310062053);
}
void car_h_30(double *state, double *unused, double *out_3516726837814677646) {
  h_30(state, unused, out_3516726837814677646);
}
void car_H_30(double *state, double *unused, double *out_9079256541749186476) {
  H_30(state, unused, out_9079256541749186476);
}
void car_h_26(double *state, double *unused, double *out_6797261298903552818) {
  h_26(state, unused, out_6797261298903552818);
}
void car_H_26(double *state, double *unused, double *out_5755323160229548986) {
  H_26(state, unused, out_5755323160229548986);
}
void car_h_27(double *state, double *unused, double *out_9119818263088886845) {
  h_27(state, unused, out_9119818263088886845);
}
void car_H_27(double *state, double *unused, double *out_7192724220159940229) {
  H_27(state, unused, out_7192724220159940229);
}
void car_h_29(double *state, double *unused, double *out_6922248888625985153) {
  h_29(state, unused, out_6922248888625985153);
}
void car_H_29(double *state, double *unused, double *out_8569025197434794292) {
  H_29(state, unused, out_8569025197434794292);
}
void car_h_28(double *state, double *unused, double *out_5703544973439816829) {
  h_28(state, unused, out_5703544973439816829);
}
void car_H_28(double *state, double *unused, double *out_396962476220858622) {
  H_28(state, unused, out_396962476220858622);
}
void car_h_31(double *state, double *unused, double *out_266242514659999675) {
  h_31(state, unused, out_266242514659999675);
}
void car_H_31(double *state, double *unused, double *out_8919271632728985978) {
  H_31(state, unused, out_8919271632728985978);
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
