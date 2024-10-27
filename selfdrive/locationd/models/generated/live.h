#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_35(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_4901497926982652881);
void live_err_fun(double *nom_x, double *delta_x, double *out_4218550092640272468);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_7473104647569890902);
void live_H_mod_fun(double *state, double *out_1172482510687413770);
void live_f_fun(double *state, double dt, double *out_6680059142450056333);
void live_F_fun(double *state, double dt, double *out_6093395981269628726);
void live_h_4(double *state, double *unused, double *out_4363039283608560512);
void live_H_4(double *state, double *unused, double *out_1208285306707785397);
void live_h_9(double *state, double *unused, double *out_114298876477406303);
void live_H_9(double *state, double *unused, double *out_6078933628556662073);
void live_h_10(double *state, double *unused, double *out_4017517240206656995);
void live_H_10(double *state, double *unused, double *out_974291533584954603);
void live_h_12(double *state, double *unused, double *out_4500566853436925760);
void live_H_12(double *state, double *unused, double *out_7589543683750518393);
void live_h_35(double *state, double *unused, double *out_4288276543319315052);
void live_H_35(double *state, double *unused, double *out_4843980651425504684);
void live_h_32(double *state, double *unused, double *out_5500753020028864922);
void live_H_32(double *state, double *unused, double *out_9184985427158611572);
void live_h_13(double *state, double *unused, double *out_9140390590474376911);
void live_H_13(double *state, double *unused, double *out_5368614765188541145);
void live_h_14(double *state, double *unused, double *out_114298876477406303);
void live_H_14(double *state, double *unused, double *out_6078933628556662073);
void live_h_33(double *state, double *unused, double *out_6487249515284105110);
void live_H_33(double *state, double *unused, double *out_1693423646786647080);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}