#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_5464861526568687019);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_7516906840074992970);
void car_H_mod_fun(double *state, double *out_6431729039895127708);
void car_f_fun(double *state, double dt, double *out_57451106597728332);
void car_F_fun(double *state, double dt, double *out_447865623327851108);
void car_h_25(double *state, double *unused, double *out_4364545812939371539);
void car_H_25(double *state, double *unused, double *out_2450797190468748385);
void car_h_24(double *state, double *unused, double *out_4446211979141650588);
void car_H_24(double *state, double *unused, double *out_3434863061310062053);
void car_h_30(double *state, double *unused, double *out_3516726837814677646);
void car_H_30(double *state, double *unused, double *out_9079256541749186476);
void car_h_26(double *state, double *unused, double *out_6797261298903552818);
void car_H_26(double *state, double *unused, double *out_5755323160229548986);
void car_h_27(double *state, double *unused, double *out_9119818263088886845);
void car_H_27(double *state, double *unused, double *out_7192724220159940229);
void car_h_29(double *state, double *unused, double *out_6922248888625985153);
void car_H_29(double *state, double *unused, double *out_8569025197434794292);
void car_h_28(double *state, double *unused, double *out_5703544973439816829);
void car_H_28(double *state, double *unused, double *out_396962476220858622);
void car_h_31(double *state, double *unused, double *out_266242514659999675);
void car_H_31(double *state, double *unused, double *out_8919271632728985978);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}