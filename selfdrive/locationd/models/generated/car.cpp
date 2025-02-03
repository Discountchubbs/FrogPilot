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
void err_fun(double *nom_x, double *delta_x, double *out_9171236554461587422) {
   out_9171236554461587422[0] = delta_x[0] + nom_x[0];
   out_9171236554461587422[1] = delta_x[1] + nom_x[1];
   out_9171236554461587422[2] = delta_x[2] + nom_x[2];
   out_9171236554461587422[3] = delta_x[3] + nom_x[3];
   out_9171236554461587422[4] = delta_x[4] + nom_x[4];
   out_9171236554461587422[5] = delta_x[5] + nom_x[5];
   out_9171236554461587422[6] = delta_x[6] + nom_x[6];
   out_9171236554461587422[7] = delta_x[7] + nom_x[7];
   out_9171236554461587422[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_6424056778957097950) {
   out_6424056778957097950[0] = -nom_x[0] + true_x[0];
   out_6424056778957097950[1] = -nom_x[1] + true_x[1];
   out_6424056778957097950[2] = -nom_x[2] + true_x[2];
   out_6424056778957097950[3] = -nom_x[3] + true_x[3];
   out_6424056778957097950[4] = -nom_x[4] + true_x[4];
   out_6424056778957097950[5] = -nom_x[5] + true_x[5];
   out_6424056778957097950[6] = -nom_x[6] + true_x[6];
   out_6424056778957097950[7] = -nom_x[7] + true_x[7];
   out_6424056778957097950[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_2838402224479150806) {
   out_2838402224479150806[0] = 1.0;
   out_2838402224479150806[1] = 0;
   out_2838402224479150806[2] = 0;
   out_2838402224479150806[3] = 0;
   out_2838402224479150806[4] = 0;
   out_2838402224479150806[5] = 0;
   out_2838402224479150806[6] = 0;
   out_2838402224479150806[7] = 0;
   out_2838402224479150806[8] = 0;
   out_2838402224479150806[9] = 0;
   out_2838402224479150806[10] = 1.0;
   out_2838402224479150806[11] = 0;
   out_2838402224479150806[12] = 0;
   out_2838402224479150806[13] = 0;
   out_2838402224479150806[14] = 0;
   out_2838402224479150806[15] = 0;
   out_2838402224479150806[16] = 0;
   out_2838402224479150806[17] = 0;
   out_2838402224479150806[18] = 0;
   out_2838402224479150806[19] = 0;
   out_2838402224479150806[20] = 1.0;
   out_2838402224479150806[21] = 0;
   out_2838402224479150806[22] = 0;
   out_2838402224479150806[23] = 0;
   out_2838402224479150806[24] = 0;
   out_2838402224479150806[25] = 0;
   out_2838402224479150806[26] = 0;
   out_2838402224479150806[27] = 0;
   out_2838402224479150806[28] = 0;
   out_2838402224479150806[29] = 0;
   out_2838402224479150806[30] = 1.0;
   out_2838402224479150806[31] = 0;
   out_2838402224479150806[32] = 0;
   out_2838402224479150806[33] = 0;
   out_2838402224479150806[34] = 0;
   out_2838402224479150806[35] = 0;
   out_2838402224479150806[36] = 0;
   out_2838402224479150806[37] = 0;
   out_2838402224479150806[38] = 0;
   out_2838402224479150806[39] = 0;
   out_2838402224479150806[40] = 1.0;
   out_2838402224479150806[41] = 0;
   out_2838402224479150806[42] = 0;
   out_2838402224479150806[43] = 0;
   out_2838402224479150806[44] = 0;
   out_2838402224479150806[45] = 0;
   out_2838402224479150806[46] = 0;
   out_2838402224479150806[47] = 0;
   out_2838402224479150806[48] = 0;
   out_2838402224479150806[49] = 0;
   out_2838402224479150806[50] = 1.0;
   out_2838402224479150806[51] = 0;
   out_2838402224479150806[52] = 0;
   out_2838402224479150806[53] = 0;
   out_2838402224479150806[54] = 0;
   out_2838402224479150806[55] = 0;
   out_2838402224479150806[56] = 0;
   out_2838402224479150806[57] = 0;
   out_2838402224479150806[58] = 0;
   out_2838402224479150806[59] = 0;
   out_2838402224479150806[60] = 1.0;
   out_2838402224479150806[61] = 0;
   out_2838402224479150806[62] = 0;
   out_2838402224479150806[63] = 0;
   out_2838402224479150806[64] = 0;
   out_2838402224479150806[65] = 0;
   out_2838402224479150806[66] = 0;
   out_2838402224479150806[67] = 0;
   out_2838402224479150806[68] = 0;
   out_2838402224479150806[69] = 0;
   out_2838402224479150806[70] = 1.0;
   out_2838402224479150806[71] = 0;
   out_2838402224479150806[72] = 0;
   out_2838402224479150806[73] = 0;
   out_2838402224479150806[74] = 0;
   out_2838402224479150806[75] = 0;
   out_2838402224479150806[76] = 0;
   out_2838402224479150806[77] = 0;
   out_2838402224479150806[78] = 0;
   out_2838402224479150806[79] = 0;
   out_2838402224479150806[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_7093649650184811885) {
   out_7093649650184811885[0] = state[0];
   out_7093649650184811885[1] = state[1];
   out_7093649650184811885[2] = state[2];
   out_7093649650184811885[3] = state[3];
   out_7093649650184811885[4] = state[4];
   out_7093649650184811885[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_7093649650184811885[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_7093649650184811885[7] = state[7];
   out_7093649650184811885[8] = state[8];
}
void F_fun(double *state, double dt, double *out_466002067398569504) {
   out_466002067398569504[0] = 1;
   out_466002067398569504[1] = 0;
   out_466002067398569504[2] = 0;
   out_466002067398569504[3] = 0;
   out_466002067398569504[4] = 0;
   out_466002067398569504[5] = 0;
   out_466002067398569504[6] = 0;
   out_466002067398569504[7] = 0;
   out_466002067398569504[8] = 0;
   out_466002067398569504[9] = 0;
   out_466002067398569504[10] = 1;
   out_466002067398569504[11] = 0;
   out_466002067398569504[12] = 0;
   out_466002067398569504[13] = 0;
   out_466002067398569504[14] = 0;
   out_466002067398569504[15] = 0;
   out_466002067398569504[16] = 0;
   out_466002067398569504[17] = 0;
   out_466002067398569504[18] = 0;
   out_466002067398569504[19] = 0;
   out_466002067398569504[20] = 1;
   out_466002067398569504[21] = 0;
   out_466002067398569504[22] = 0;
   out_466002067398569504[23] = 0;
   out_466002067398569504[24] = 0;
   out_466002067398569504[25] = 0;
   out_466002067398569504[26] = 0;
   out_466002067398569504[27] = 0;
   out_466002067398569504[28] = 0;
   out_466002067398569504[29] = 0;
   out_466002067398569504[30] = 1;
   out_466002067398569504[31] = 0;
   out_466002067398569504[32] = 0;
   out_466002067398569504[33] = 0;
   out_466002067398569504[34] = 0;
   out_466002067398569504[35] = 0;
   out_466002067398569504[36] = 0;
   out_466002067398569504[37] = 0;
   out_466002067398569504[38] = 0;
   out_466002067398569504[39] = 0;
   out_466002067398569504[40] = 1;
   out_466002067398569504[41] = 0;
   out_466002067398569504[42] = 0;
   out_466002067398569504[43] = 0;
   out_466002067398569504[44] = 0;
   out_466002067398569504[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_466002067398569504[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_466002067398569504[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_466002067398569504[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_466002067398569504[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_466002067398569504[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_466002067398569504[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_466002067398569504[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_466002067398569504[53] = -9.8000000000000007*dt;
   out_466002067398569504[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_466002067398569504[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_466002067398569504[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_466002067398569504[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_466002067398569504[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_466002067398569504[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_466002067398569504[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_466002067398569504[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_466002067398569504[62] = 0;
   out_466002067398569504[63] = 0;
   out_466002067398569504[64] = 0;
   out_466002067398569504[65] = 0;
   out_466002067398569504[66] = 0;
   out_466002067398569504[67] = 0;
   out_466002067398569504[68] = 0;
   out_466002067398569504[69] = 0;
   out_466002067398569504[70] = 1;
   out_466002067398569504[71] = 0;
   out_466002067398569504[72] = 0;
   out_466002067398569504[73] = 0;
   out_466002067398569504[74] = 0;
   out_466002067398569504[75] = 0;
   out_466002067398569504[76] = 0;
   out_466002067398569504[77] = 0;
   out_466002067398569504[78] = 0;
   out_466002067398569504[79] = 0;
   out_466002067398569504[80] = 1;
}
void h_25(double *state, double *unused, double *out_8685010582261463285) {
   out_8685010582261463285[0] = state[6];
}
void H_25(double *state, double *unused, double *out_1601500309678166422) {
   out_1601500309678166422[0] = 0;
   out_1601500309678166422[1] = 0;
   out_1601500309678166422[2] = 0;
   out_1601500309678166422[3] = 0;
   out_1601500309678166422[4] = 0;
   out_1601500309678166422[5] = 0;
   out_1601500309678166422[6] = 1;
   out_1601500309678166422[7] = 0;
   out_1601500309678166422[8] = 0;
}
void h_24(double *state, double *unused, double *out_4131104852282051912) {
   out_4131104852282051912[0] = state[4];
   out_4131104852282051912[1] = state[5];
}
void H_24(double *state, double *unused, double *out_5159425908480489558) {
   out_5159425908480489558[0] = 0;
   out_5159425908480489558[1] = 0;
   out_5159425908480489558[2] = 0;
   out_5159425908480489558[3] = 0;
   out_5159425908480489558[4] = 1;
   out_5159425908480489558[5] = 0;
   out_5159425908480489558[6] = 0;
   out_5159425908480489558[7] = 0;
   out_5159425908480489558[8] = 0;
   out_5159425908480489558[9] = 0;
   out_5159425908480489558[10] = 0;
   out_5159425908480489558[11] = 0;
   out_5159425908480489558[12] = 0;
   out_5159425908480489558[13] = 0;
   out_5159425908480489558[14] = 1;
   out_5159425908480489558[15] = 0;
   out_5159425908480489558[16] = 0;
   out_5159425908480489558[17] = 0;
}
void h_30(double *state, double *unused, double *out_1914175355911112349) {
   out_1914175355911112349[0] = state[4];
}
void H_30(double *state, double *unused, double *out_8518190651169783177) {
   out_8518190651169783177[0] = 0;
   out_8518190651169783177[1] = 0;
   out_8518190651169783177[2] = 0;
   out_8518190651169783177[3] = 0;
   out_8518190651169783177[4] = 1;
   out_8518190651169783177[5] = 0;
   out_8518190651169783177[6] = 0;
   out_8518190651169783177[7] = 0;
   out_8518190651169783177[8] = 0;
}
void h_26(double *state, double *unused, double *out_3087497911749733296) {
   out_3087497911749733296[0] = state[7];
}
void H_26(double *state, double *unused, double *out_2140003009195889802) {
   out_2140003009195889802[0] = 0;
   out_2140003009195889802[1] = 0;
   out_2140003009195889802[2] = 0;
   out_2140003009195889802[3] = 0;
   out_2140003009195889802[4] = 0;
   out_2140003009195889802[5] = 0;
   out_2140003009195889802[6] = 0;
   out_2140003009195889802[7] = 1;
   out_2140003009195889802[8] = 0;
}
void h_27(double *state, double *unused, double *out_5363781410166459319) {
   out_5363781410166459319[0] = state[3];
}
void H_27(double *state, double *unused, double *out_6343427339369358266) {
   out_6343427339369358266[0] = 0;
   out_6343427339369358266[1] = 0;
   out_6343427339369358266[2] = 0;
   out_6343427339369358266[3] = 1;
   out_6343427339369358266[4] = 0;
   out_6343427339369358266[5] = 0;
   out_6343427339369358266[6] = 0;
   out_6343427339369358266[7] = 0;
   out_6343427339369358266[8] = 0;
}
void h_29(double *state, double *unused, double *out_1502097933737957667) {
   out_1502097933737957667[0] = state[1];
}
void H_29(double *state, double *unused, double *out_9028421995484175361) {
   out_9028421995484175361[0] = 0;
   out_9028421995484175361[1] = 1;
   out_9028421995484175361[2] = 0;
   out_9028421995484175361[3] = 0;
   out_9028421995484175361[4] = 0;
   out_9028421995484175361[5] = 0;
   out_9028421995484175361[6] = 0;
   out_9028421995484175361[7] = 0;
   out_9028421995484175361[8] = 0;
}
void h_28(double *state, double *unused, double *out_6307291263068121754) {
   out_6307291263068121754[0] = state[0];
}
void H_28(double *state, double *unused, double *out_452334404569723341) {
   out_452334404569723341[0] = 1;
   out_452334404569723341[1] = 0;
   out_452334404569723341[2] = 0;
   out_452334404569723341[3] = 0;
   out_452334404569723341[4] = 0;
   out_452334404569723341[5] = 0;
   out_452334404569723341[6] = 0;
   out_452334404569723341[7] = 0;
   out_452334404569723341[8] = 0;
}
void h_31(double *state, double *unused, double *out_4443839867628224302) {
   out_4443839867628224302[0] = state[8];
}
void H_31(double *state, double *unused, double *out_1632146271555126850) {
   out_1632146271555126850[0] = 0;
   out_1632146271555126850[1] = 0;
   out_1632146271555126850[2] = 0;
   out_1632146271555126850[3] = 0;
   out_1632146271555126850[4] = 0;
   out_1632146271555126850[5] = 0;
   out_1632146271555126850[6] = 0;
   out_1632146271555126850[7] = 0;
   out_1632146271555126850[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_9171236554461587422) {
  err_fun(nom_x, delta_x, out_9171236554461587422);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_6424056778957097950) {
  inv_err_fun(nom_x, true_x, out_6424056778957097950);
}
void car_H_mod_fun(double *state, double *out_2838402224479150806) {
  H_mod_fun(state, out_2838402224479150806);
}
void car_f_fun(double *state, double dt, double *out_7093649650184811885) {
  f_fun(state,  dt, out_7093649650184811885);
}
void car_F_fun(double *state, double dt, double *out_466002067398569504) {
  F_fun(state,  dt, out_466002067398569504);
}
void car_h_25(double *state, double *unused, double *out_8685010582261463285) {
  h_25(state, unused, out_8685010582261463285);
}
void car_H_25(double *state, double *unused, double *out_1601500309678166422) {
  H_25(state, unused, out_1601500309678166422);
}
void car_h_24(double *state, double *unused, double *out_4131104852282051912) {
  h_24(state, unused, out_4131104852282051912);
}
void car_H_24(double *state, double *unused, double *out_5159425908480489558) {
  H_24(state, unused, out_5159425908480489558);
}
void car_h_30(double *state, double *unused, double *out_1914175355911112349) {
  h_30(state, unused, out_1914175355911112349);
}
void car_H_30(double *state, double *unused, double *out_8518190651169783177) {
  H_30(state, unused, out_8518190651169783177);
}
void car_h_26(double *state, double *unused, double *out_3087497911749733296) {
  h_26(state, unused, out_3087497911749733296);
}
void car_H_26(double *state, double *unused, double *out_2140003009195889802) {
  H_26(state, unused, out_2140003009195889802);
}
void car_h_27(double *state, double *unused, double *out_5363781410166459319) {
  h_27(state, unused, out_5363781410166459319);
}
void car_H_27(double *state, double *unused, double *out_6343427339369358266) {
  H_27(state, unused, out_6343427339369358266);
}
void car_h_29(double *state, double *unused, double *out_1502097933737957667) {
  h_29(state, unused, out_1502097933737957667);
}
void car_H_29(double *state, double *unused, double *out_9028421995484175361) {
  H_29(state, unused, out_9028421995484175361);
}
void car_h_28(double *state, double *unused, double *out_6307291263068121754) {
  h_28(state, unused, out_6307291263068121754);
}
void car_H_28(double *state, double *unused, double *out_452334404569723341) {
  H_28(state, unused, out_452334404569723341);
}
void car_h_31(double *state, double *unused, double *out_4443839867628224302) {
  h_31(state, unused, out_4443839867628224302);
}
void car_H_31(double *state, double *unused, double *out_1632146271555126850) {
  H_31(state, unused, out_1632146271555126850);
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
