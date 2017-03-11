#define EIGEN_NO_DEBUG

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <iomanip>
#include "spline.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/KroneckerProduct>

using namespace std;
using namespace spline;
using namespace Eigen;
const double pi = 3.14159265358979323846264338327;
typedef Eigen::SparseVector<double> SpVec;
typedef Eigen::SparseMatrix<double> SpMat;

void GenerateKnot(int order, int refine, int insert, int repeat, double *insert_knot,
                  double *&knot, int &m);

void Geometry1(double xi, double eta, double &pxpxi, double &pxpeta, double &pypxi, double &pypeta,
               double &pxpxi_xi, double &pxpxi_eta, double &pxpeta_eta, double &pypxi_xi, double &pypxi_eta,
               double &pypeta_eta,
               double &x, double &y);

void Geometry2(double xi, double eta, double &pxpxi, double &pxpeta, double &pypxi, double &pypeta,
               double &pxpxi_xi, double &pxpxi_eta, double &pxpeta_eta, double &pypxi_xi, double &pypxi_eta,
               double &pypeta_eta,
               double &x, double &y);

void CombineKnot(double *knot1, int m1, double *knot2, int m2, double *&knot,
                 int &m) {
    vector<double> Knot;
    int i1 = 0, i2 = 0;
    while (i1 <= m1) {
        if (knot1[i1] == knot2[i2]) {
            Knot.push_back(knot1[i1]);
            i1++, i2++;
        } else if (knot1[i1] < knot2[i2]) {
            Knot.push_back(knot1[i1]);
            i1++;
        } else {
            Knot.push_back(knot2[i2]);
            i2++;
        }
    }

    m = Knot.end() - Knot.begin() - 1;
    knot = new double[m + 1];
    for (int i = 0; i < m + 1; i++)
        knot[i] = Knot[i];
}

double exactSolution(double x, double y) {
    return pow((x - 4) * (x - 2) * (y - 4) * (y - 2) * x * y, 2);
}
double forceTerm(double x, double y) {
    double s1, s2, s3, s4, s5, s6, s7;
    s1 = -36 * pow(x, 5) * (52 - 60 * y + 15 * y * y);
    s2 = 3 * pow(x, 6) * (52 - 60 * y + 15 * y * y);
    s3 = 9 * pow(x, 4) * (1008 - 1520 * y + 780 * y * y - 200 * pow(y, 3) + 25 * pow(y, 4));
    s4 = -24 * pow(x, 3) * (944 - 2160 * y + 1740 * y * y - 600 * pow(y, 3) + 75 * pow(y, 4));
    s5 = -36 * x * (512 - 2304 * y + 2816 * y * y - 1440 * pow(y, 3) + 380 * pow(y, 4) - 60 * pow(y, 5) + 5 * pow(y, 6));
    s6 = 9 * x * x * (3328 - 11264 * y + 11456 * y * y - 4640 * pow(y, 3) + 780 * pow(y, 4) - 60 * pow(y, 5) + 5 * pow(y, 6));
    s7 = 4 * (1024 - 4608 * y + 7488 * y * y - 5664 * pow(y, 3) + 2268 * pow(y, 4) - 468 * pow(y, 5) + 39 * pow(y, 6));
    return 8 * (s1 + s2 + s3 + s4 + s5 + s6 + s7);
}

int main() {
    Eigen::setNbThreads(8);
    double *gaussian = x7;
    double *weight = w7;
    int gaussian_points = 7;
    int order;
    int refine;
    cin >> order >> refine;
    int m_x_patch1, m_y_patch1;
    int m_x_patch2, m_y_patch2, m_y_patch2_projection;
    int p_x = order, p_y = order;
    double *knots_x_patch1, *knots_y_patch1;
    double *knots_x_patch2, *knots_y_patch2;
    double *knots_y_patch2_projection;
    double insertion_patch1[] = {1.0 / 2};
    double insertion_patch2[] = {1.0 / 3, 2.0 / 3};
    GenerateKnot(p_x, refine, 0, 1, insertion_patch1, knots_x_patch1, m_x_patch1);
    GenerateKnot(p_y, refine, 0, 1, insertion_patch1, knots_y_patch1, m_y_patch1);
    GenerateKnot(p_x, refine, 0, 1, insertion_patch2, knots_x_patch2, m_x_patch2);
    GenerateKnot(p_y, refine, 2, 1, insertion_patch2, knots_y_patch2, m_y_patch2);
    double *knots_y_coupling;
    int m_y_coupling;
    CombineKnot(knots_y_patch1, m_y_patch1, knots_y_patch2, m_y_patch2, knots_y_coupling, m_y_coupling);
    const int dof_x_patch1 = m_x_patch1 - p_x, dof_y_patch1 = m_y_patch1 - p_y,
            dof_x_patch2 = m_x_patch2 - p_x,
            dof_y_patch2 = m_y_patch2 - p_y;
    const int dof_patch1 = dof_x_patch1 * dof_y_patch1,
            dof_patch2 = dof_x_patch2 * dof_y_patch2;
    const int elements_x_patch1 = m_x_patch1 - 2 * p_x,
            elements_y_patch1 = m_y_patch1 - 2 * p_y,
            elements_x_patch2 = m_x_patch2 - 2 * p_x,
            elements_y_patch2 = m_y_patch2 - 2 * p_y;
    const int elements_y_coupling = m_y_coupling - 2 * p_y;
    MatrixXd M_C0 = MatrixXd::Zero(dof_y_patch2, dof_y_patch2), N2N1_C0 = MatrixXd::Zero(dof_y_patch2, dof_y_patch1);
    for (int ii_y = 0; ii_y < elements_y_coupling; ii_y++) {
        double J_y = (knots_y_coupling[ii_y + p_y + 1] - knots_y_coupling[ii_y + p_y]) / 2;
        double Middle_y = (knots_y_coupling[ii_y + p_y + 1] + knots_y_coupling[ii_y + p_y]) / 2;
        int i_y_patch1 = Findspan(m_y_patch1, p_y, knots_y_patch1, Middle_y);
        int i_y_patch2 = Findspan(m_y_patch2, p_y, knots_y_patch2, Middle_y);
        for (int jj_y = 0; jj_y < gaussian_points; jj_y++) {
            double eta = Middle_y + J_y * gaussian[jj_y];
            double xi = 0;
            double **ders_y_patch1, **ders_y_patch2;
            DersBasisFuns(i_y_patch1, eta, p_y, knots_y_patch1, 0, ders_y_patch1);
            DersBasisFuns(i_y_patch2, eta, p_y, knots_y_patch2, 0, ders_y_patch2);
            VectorXd Neta_patch1 = VectorXd::Zero(dof_y_patch1),
                    Neta_patch2 = VectorXd::Zero(dof_y_patch2);
            for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
                Neta_patch1(i_y_patch1 - p_y + kk_y) = ders_y_patch1[0][kk_y];
                Neta_patch2(i_y_patch2 - p_y + kk_y) = ders_y_patch2[0][kk_y];
            }
            for (int k = 0; k < 1; k++)
                delete ders_y_patch1[k];
            delete[] ders_y_patch1;
            for (int k = 0; k < 1; k++)
                delete ders_y_patch2[k];
            delete[] ders_y_patch2;
            double pxpxi, pxpeta, pypxi, pypeta, pxpxi_xi, pxpxi_eta, pxpeta_eta, pypxi_xi, pypxi_eta, pypeta_eta, x, y;
            Geometry2(xi, eta, pxpxi, pxpeta, pypxi, pypeta, pxpxi_xi, pxpxi_eta, pxpeta_eta, pypxi_xi, pypxi_eta,
                      pypeta_eta, x,
                      y);
            M_C0 += weight[jj_y] * Neta_patch2 * Neta_patch2.transpose() * J_y *
                    pow(pxpeta * pxpeta + pypeta * pypeta, .5);
            N2N1_C0 += weight[jj_y] * Neta_patch2 * Neta_patch1.transpose() * J_y *
                       pow(pxpeta * pxpeta + pypeta * pypeta, .5);
        }
    }
    MatrixXd C0_constrain(dof_y_patch2, dof_y_patch1);
    C0_constrain.setZero();
    C0_constrain.block(0, 0, dof_y_patch2 - 2, 1) = M_C0.block(0, 0, dof_y_patch2 - 2,
                                                               dof_y_patch2 - 2).partialPivLu().solve(
            N2N1_C0.block(0, 0, dof_y_patch2 - 2, 1));
    C0_constrain.block(1, 1, dof_y_patch2 - 3, 1) = M_C0.block(1, 1, dof_y_patch2 - 3,
                                                               dof_y_patch2 - 3).partialPivLu().solve(
            N2N1_C0.block(1, 1, dof_y_patch2 - 3, 1));
    C0_constrain.block(2, 2, dof_y_patch2 - 4, dof_y_patch1 - 4) = M_C0.block(2, 2, dof_y_patch2 - 4,
                                                                              dof_y_patch2 - 4).partialPivLu().solve(
            N2N1_C0.block(2, 2, dof_y_patch2 - 4, dof_y_patch1 - 4));
    C0_constrain.block(2, dof_y_patch1 - 2, dof_y_patch2 - 3, 1) = M_C0.block(2, 2, dof_y_patch2 - 3,
                                                                              dof_y_patch2 - 3).partialPivLu().solve(
            N2N1_C0.block(2, dof_y_patch1 - 2, dof_y_patch2 - 3, 1));
    C0_constrain.block(2, dof_y_patch1 - 1, dof_y_patch2 - 2, 1) = M_C0.block(2, 2, dof_y_patch2 - 2,
                                                                              dof_y_patch2 - 2).partialPivLu().solve(
            N2N1_C0.block(2, dof_y_patch1 - 1, dof_y_patch2 - 2, 1));
    cout << C0_constrain << endl;

    MatrixXd M_C1 = MatrixXd::Zero(dof_y_patch2, dof_y_patch2), N2N1_patch1 = MatrixXd::Zero(dof_y_patch2,
                                                                                             dof_y_patch1), N2N2_patch1 = MatrixXd::Zero(
            dof_y_patch2, dof_y_patch1),
            N1N1_patch1 = MatrixXd::Zero(dof_y_patch2, dof_y_patch1), N1N1_patch2 = MatrixXd::Zero(dof_y_patch2,
                                                                                                   dof_y_patch2);

    for (int ii_y = 0; ii_y < elements_y_coupling; ii_y++) {
        double J_y = (knots_y_coupling[ii_y + p_y + 1] - knots_y_coupling[ii_y + p_y]) / 2;
        double Middle_y = (knots_y_coupling[ii_y + p_y + 1] + knots_y_coupling[ii_y + p_y]) / 2;
        int i_y_patch1 = Findspan(m_y_patch1, p_y, knots_y_patch1, Middle_y);
        int i_y_patch2 = Findspan(m_y_patch2, p_y, knots_y_patch2, Middle_y);
        for (int jj_y = 0; jj_y < gaussian_points; jj_y++) {
            double eta_patch1, eta_patch2;
            eta_patch1 = eta_patch2 = Middle_y + J_y * gaussian[jj_y];
            double xi_patch1 = .999999999999999, xi_patch2 = 0;
            double **ders_y_patch1, **ders_y_patch2;
            double nxi_xi_second_last_patch1 = DersOneBasisFUn(p_x, m_x_patch1, knots_x_patch1, m_x_patch1 - p_x - 2,
                                                               xi_patch1, 1);
            double nxi_second_last_patch1 = DersOneBasisFUn(p_x, m_x_patch1, knots_x_patch1, m_x_patch1 - p_x - 2,
                                                            xi_patch1, 0);
            double nxi_xi_second_last_patch2 = DersOneBasisFUn(p_x, m_x_patch2, knots_x_patch2, 1, xi_patch2, 1);
            double nxi_second_last_patch2 = DersOneBasisFUn(p_x, m_x_patch2, knots_x_patch2, 1, xi_patch2, 0);
            double nxi_xi_last_patch1 = DersOneBasisFUn(p_x, m_x_patch1, knots_x_patch1, m_x_patch1 - p_x - 1,
                                                        xi_patch1, 1);
            double nxi_last_patch1 = DersOneBasisFUn(p_x, m_x_patch1, knots_x_patch1, m_x_patch1 - p_x - 1, xi_patch1,
                                                     0);
            double nxi_xi_last_patch2 = DersOneBasisFUn(p_x, m_x_patch2, knots_x_patch2, 0, xi_patch2, 1);
            double nxi_last_patch2 = DersOneBasisFUn(p_x, m_x_patch2, knots_x_patch2, 0, xi_patch2, 0);
            DersBasisFuns(i_y_patch1, eta_patch1, p_y, knots_y_patch1, 1, ders_y_patch1);
            DersBasisFuns(i_y_patch2, eta_patch2, p_y, knots_y_patch2, 1, ders_y_patch2);
            VectorXd Nxi_xiNeta_second_last_patch1 = VectorXd::Zero(dof_y_patch1),
                    NxiNeta_eta_second_last_patch1 = VectorXd::Zero(dof_y_patch1),
                    Nxi_xiNeta_second_last_patch2 = VectorXd::Zero(dof_y_patch2),
                    NxiNeta_eta_second_last_patch2 = VectorXd::Zero(dof_y_patch2),
                    NxiNeta_second_last_patch2 = VectorXd::Zero(dof_y_patch2), Nxi_xiNeta_last_patch1 = VectorXd::Zero(
                    dof_y_patch1),
                    NxiNeta_eta_last_patch1 = VectorXd::Zero(dof_y_patch1), Nxi_xiNeta_last_patch2 = VectorXd::Zero(
                    dof_y_patch2),
                    NxiNeta_eta_last_patch2 = VectorXd::Zero(dof_y_patch2), NxiNeta_last_patch2 = VectorXd::Zero(
                    dof_y_patch2), NxiNeta_last_patch1 = VectorXd::Zero(dof_y_patch1);
            for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
                Nxi_xiNeta_second_last_patch1(i_y_patch1 - p_y + kk_y) =
                        ders_y_patch1[0][kk_y] * nxi_xi_second_last_patch1;
                NxiNeta_eta_second_last_patch1(i_y_patch1 - p_y + kk_y) =
                        ders_y_patch1[1][kk_y] * nxi_second_last_patch1;
                Nxi_xiNeta_second_last_patch2(i_y_patch2 - p_y + kk_y) =
                        ders_y_patch2[0][kk_y] * nxi_xi_second_last_patch2;
                NxiNeta_eta_second_last_patch2(i_y_patch2 - p_y + kk_y) =
                        ders_y_patch2[1][kk_y] * nxi_second_last_patch2;
                NxiNeta_last_patch2(i_y_patch2 - p_y + kk_y) = ders_y_patch2[0][kk_y];
                Nxi_xiNeta_last_patch1(i_y_patch1 - p_y + kk_y) = ders_y_patch1[0][kk_y] * nxi_xi_last_patch1;
                NxiNeta_eta_last_patch1(i_y_patch1 - p_y + kk_y) = ders_y_patch1[1][kk_y] * nxi_last_patch1;
                Nxi_xiNeta_last_patch2(i_y_patch2 - p_y + kk_y) = ders_y_patch2[0][kk_y] * nxi_xi_last_patch2;
                NxiNeta_eta_last_patch2(i_y_patch2 - p_y + kk_y) = ders_y_patch2[1][kk_y] * nxi_last_patch2;
                NxiNeta_last_patch1(i_y_patch1 - p_y + kk_y) = ders_y_patch1[0][kk_y];
            }
            for (int k = 0; k < 2; k++)
                delete ders_y_patch1[k];
            delete[] ders_y_patch1;
            for (int k = 0; k < 2; k++)
                delete ders_y_patch2[k];
            delete[] ders_y_patch2;
            double pxpxi_patch1, pxpeta_patch1, pypxi_patch1, pypeta_patch1, x_patch1, y_patch1, pxpxi_patch2, pxpeta_patch2,
                    pypxi_patch2, pypeta_patch2, x_patch2, y_patch2, pxpxi_xi_patch1, pxpxi_eta_patch1, pxpeta_eta_patch1, pypxi_xi_patch1,
                    pypxi_eta_patch1, pypeta_eta_patch1, pxpxi_xi_patch2,
                    pxpxi_eta_patch2, pxpeta_eta_patch2, pypxi_xi_patch2, pypxi_eta_patch2, pypeta_eta_patch2;
            Geometry1(xi_patch1, eta_patch1, pxpxi_patch1, pxpeta_patch1, pypxi_patch1, pypeta_patch1, pxpxi_xi_patch1,
                      pxpxi_eta_patch1, pxpeta_eta_patch1, pypxi_xi_patch1, pypxi_eta_patch1, pypeta_eta_patch1,
                      x_patch1, y_patch1);
            Geometry2(xi_patch2, eta_patch2, pxpxi_patch2, pxpeta_patch2, pypxi_patch2, pypeta_patch2, pxpxi_xi_patch2,
                      pxpxi_eta_patch2, pxpeta_eta_patch2, pypxi_xi_patch2, pypxi_eta_patch2, pypeta_eta_patch2,
                      x_patch2, y_patch2);
            double alpha, beta, gamma;
            alpha = pxpxi_patch2 * pypxi_patch1 - pxpxi_patch1 * pypxi_patch2;
            beta = pxpeta_patch1 * pypxi_patch2 - pxpxi_patch2 * pypeta_patch1;
            gamma = pxpeta_patch1 * pypxi_patch1 - pxpxi_patch1 * pypeta_patch1;
            M_C1 += weight[jj_y] * NxiNeta_last_patch2 * Nxi_xiNeta_second_last_patch2.transpose() * J_y;
            N2N1_patch1 +=
                    beta / gamma * weight[jj_y] * NxiNeta_last_patch2 * Nxi_xiNeta_second_last_patch1.transpose() * J_y;
            N1N1_patch1 += beta / gamma * weight[jj_y] * NxiNeta_last_patch2 * Nxi_xiNeta_last_patch1.transpose() * J_y;
            N1N1_patch2 += weight[jj_y] * NxiNeta_last_patch2 * Nxi_xiNeta_last_patch2.transpose() * J_y;
            N2N2_patch1 +=
                    alpha / gamma * weight[jj_y] * NxiNeta_last_patch2 * NxiNeta_eta_last_patch1.transpose() * J_y;
        }
    }
    MatrixXd C1_constrain_ker(dof_y_patch2, dof_y_patch1);
    C1_constrain_ker.setZero();
    C1_constrain_ker.block(0, 0, dof_y_patch2 - 2, 1) = M_C1.block(0, 0, dof_y_patch2 - 2,
                                                                   dof_y_patch2 - 2).partialPivLu().solve(
            N2N1_patch1.block(0, 0, dof_y_patch2 - 2, 1));
    C1_constrain_ker.block(1, 1, dof_y_patch2 - 3, 1) = M_C1.block(1, 1, dof_y_patch2 - 3,
                                                                   dof_y_patch2 - 3).partialPivLu().solve(
            N2N1_patch1.block(1, 1, dof_y_patch2 - 3, 1));
    C1_constrain_ker.block(2, 2, dof_y_patch2 - 4, dof_y_patch1 - 4) = M_C1.block(2, 2, dof_y_patch2 - 4,
                                                                                  dof_y_patch2 -
                                                                                  4).partialPivLu().solve(
            N2N1_patch1.block(2, 2, dof_y_patch2 - 4, dof_y_patch1 - 4));

    C1_constrain_ker.block(2, dof_y_patch1 - 2, dof_y_patch2 - 3, 1) = M_C1.block(2, 2, dof_y_patch2 - 3,
                                                                                  dof_y_patch2 -
                                                                                  3).partialPivLu().solve(
            N2N1_patch1.block(2, dof_y_patch1 - 2, dof_y_patch2 - 3, 1));
    C1_constrain_ker.block(2, dof_y_patch1 - 1, dof_y_patch2 - 2, 1) = M_C1.block(2, 2, dof_y_patch2 - 2,
                                                                                  dof_y_patch2 -
                                                                                  2).partialPivLu().solve(
            N2N1_patch1.block(2, dof_y_patch1 - 1, dof_y_patch2 - 2, 1));
    cout << C1_constrain_ker << endl << endl;
    MatrixXd C1_constrain_img(dof_y_patch2, dof_y_patch1);
    C1_constrain_img.setZero();
    MatrixXd C1_constrain_img_rhs = N1N1_patch1 + N2N2_patch1 - N1N1_patch2 * C0_constrain;
    C1_constrain_img.block(0, 0, dof_y_patch2 - 2, 1) = M_C1.block(0, 0, dof_y_patch2 - 2,
                                                                   dof_y_patch2 - 2).partialPivLu().solve(
            C1_constrain_img_rhs.block(0, 0, dof_y_patch2 - 2, 1));
    C1_constrain_img.block(1, 1, dof_y_patch2 - 3, 1) = M_C1.block(1, 1, dof_y_patch2 - 3,
                                                                   dof_y_patch2 - 3).partialPivLu().solve(
            C1_constrain_img_rhs.block(1, 1, dof_y_patch2 - 3, 1));
    C1_constrain_img.block(2, 2, dof_y_patch2 - 4, dof_y_patch1 - 4) = M_C1.block(2, 2, dof_y_patch2 - 4,
                                                                                  dof_y_patch2 -
                                                                                  4).partialPivLu().solve(
            C1_constrain_img_rhs.block(2, 2, dof_y_patch2 - 4, dof_y_patch1 - 4));


    C1_constrain_img.block(2, dof_y_patch1 - 2, dof_y_patch2 - 3, 1) = M_C1.block(2, 2, dof_y_patch2 - 3,
                                                                                  dof_y_patch2 -
                                                                                  3).partialPivLu().solve(
            C1_constrain_img_rhs.block(2, dof_y_patch1 - 2, dof_y_patch2 - 3, 1));
    C1_constrain_img.block(2, dof_y_patch1 - 1, dof_y_patch2 - 2, 1) = M_C1.block(2, 2, dof_y_patch2 - 2, dof_y_patch2 -
                                                                                                          2).partialPivLu().solve(
            C1_constrain_img_rhs.block(2, dof_y_patch1 - 1, dof_y_patch2 - 2, 1));
    cout << C1_constrain_img << endl;
    vector<Eigen::Triplet<double> > coefficients;
    SpMat assemble_C0(dof_patch1 + dof_patch2, dof_patch1 + dof_patch2 - dof_y_patch2);
    for (int i = 0; i < dof_patch1; i++)
        coefficients.push_back(Eigen::Triplet<double>(i, i, 1));
    for (int i = 0; i < dof_y_patch1; i++) {
        for (int j = 0; j < dof_y_patch2; j++) {
            coefficients.push_back(
                    Eigen::Triplet<double>(j + dof_patch1, i + dof_patch1 - dof_y_patch1, C0_constrain(j, i)));
        }
    }
    for (int i = 0; i < dof_patch2 - dof_y_patch2; i++)
        coefficients.push_back(Eigen::Triplet<double>(i + dof_patch1 + dof_y_patch2, i + dof_patch1, 1));
    assemble_C0.setFromTriplets(coefficients.begin(), coefficients.end());

    SpMat assemble_C1(dof_patch1 + dof_patch2 - dof_y_patch2, dof_patch1 + dof_patch2 - 2 * dof_y_patch2);
    coefficients.clear();
    for (int i = 0; i < dof_patch1; i++)
        coefficients.push_back(Eigen::Triplet<double>(i, i, 1));
    for (int i = 0; i < dof_y_patch1; i++) {
        for (int j = 0; j < dof_y_patch2; j++) {
            coefficients.push_back(
                    Eigen::Triplet<double>(j + dof_patch1, i + dof_patch1 - 2 * dof_y_patch1, C1_constrain_ker(j,
                                                                                                               i)));
            coefficients.push_back(
                    Eigen::Triplet<double>(j + dof_patch1, i + dof_patch1 - dof_y_patch1, C1_constrain_img(j, i)));
        }
    }
    for (int i = 0; i < dof_patch2 - 2 * dof_y_patch2; i++)
        coefficients.push_back(Eigen::Triplet<double>(i + dof_patch1 + dof_y_patch2, i + dof_patch1, 1));
    assemble_C1.setFromTriplets(coefficients.begin(), coefficients.end());
    coefficients.clear();


    SpMat K_patch1(dof_patch1, dof_patch1);
    VectorXd F_patch1 = VectorXd::Zero(dof_patch1);

    for (int ii_x = 0; ii_x < elements_x_patch1; ii_x++) {
        double J_x = (knots_x_patch1[ii_x + p_x + 1] - knots_x_patch1[ii_x + p_x]) / 2;
        double Middle_x = (knots_x_patch1[ii_x + p_x + 1] + knots_x_patch1[ii_x +
                                                                           p_x]) / 2;
        int i_x = Findspan(m_x_patch1, p_x, knots_x_patch1, Middle_x);
        for (int ii_y = 0; ii_y < elements_y_patch1; ii_y++) {
            double J_y = (knots_y_patch1[ii_y + p_y + 1] - knots_y_patch1[ii_y + p_y]) / 2;
            double Middle_y = (knots_y_patch1[ii_y + p_y + 1] + knots_y_patch1[ii_y +
                                                                               p_y]) / 2;
            int i_y = Findspan(m_y_patch1, p_y, knots_y_patch1, Middle_y);
            for (int jj_x = 0; jj_x < gaussian_points; jj_x++) {
                for (int jj_y = 0; jj_y < gaussian_points; jj_y++) {
                    double xi = Middle_x + J_x * gaussian[jj_x];
                    double eta = Middle_y + J_y * gaussian[jj_y];
                    double **ders_x, **ders_y;
                    DersBasisFuns(i_x, xi, p_x, knots_x_patch1, 2, ders_x);
                    DersBasisFuns(i_y, eta, p_y, knots_y_patch1, 2, ders_y);
                    SpVec Nxi(dof_x_patch1), Nxi_xi(dof_x_patch1), Nxi_xi_xi(dof_x_patch1), Neta(dof_y_patch1),
                            Neta_eta(dof_y_patch1), Neta_eta_eta(dof_y_patch1);
                    for (int kk_x = 0; kk_x < p_x + 1; kk_x++) {
                        Nxi.insert(i_x - p_x + kk_x) = ders_x[0][kk_x];
                        Nxi_xi.insert(i_x - p_x + kk_x) = ders_x[1][kk_x];
                        Nxi_xi_xi.insert(i_x - p_x + kk_x) = ders_x[2][kk_x];
                    }
                    for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
                        Neta.insert(i_y - p_y + kk_y) = ders_y[0][kk_y];
                        Neta_eta.insert(i_y - p_y + kk_y) = ders_y[1][kk_y];
                        Neta_eta_eta.insert(i_y - p_y + kk_y) = ders_y[2][kk_y];
                    }
                    for (int k = 0; k < 3; k++)
                        delete ders_x[k];
                    delete[] ders_x;
                    for (int k = 0; k < 3; k++)
                        delete ders_y[k];
                    delete[] ders_y;
                    SpVec Nxi_xiNeta, NxiNeta_eta, Nxi_xi_xiNeta, NxiNeta_eta_eta, Nxi_xiNeta_eta, NxiNeta;
                    Nxi_xiNeta = (kroneckerProduct(Nxi_xi, Neta)).sparseView();
                    NxiNeta_eta = (kroneckerProduct(Nxi, Neta_eta)).sparseView();
                    Nxi_xi_xiNeta = (kroneckerProduct(Nxi_xi_xi, Neta)).sparseView();
                    NxiNeta_eta_eta = (kroneckerProduct(Nxi, Neta_eta_eta)).sparseView();
                    Nxi_xiNeta_eta = (kroneckerProduct(Nxi_xi, Neta_eta)).sparseView();
                    NxiNeta = (kroneckerProduct(Nxi, Neta)).sparseView();
                    double pxpxi, pxpeta, pypxi, pypeta, pxpxi_xi, pxpxi_eta, pxpeta_eta, pypxi_xi, pypxi_eta, pypeta_eta, x, y;
                    Geometry1(xi, eta, pxpxi, pxpeta, pypxi, pypeta, pxpxi_xi, pxpxi_eta, pxpeta_eta, pypxi_xi,
                              pypxi_eta, pypeta_eta, x,
                              y);
                    double force = forceTerm(x, y);
                    double Jacobian = pxpxi * pypeta - pxpeta * pypxi;
                    MatrixXd Hessian(5, 5);
                    Hessian << pxpxi, pypxi, 0, 0, 0, pxpeta, pypeta, 0, 0, 0, pxpxi_xi, pypxi_xi, pxpxi * pxpxi, 2 *
                                                                                                                  pxpxi *
                                                                                                                  pypxi,
                            pypxi * pypxi, pxpxi_eta, pypxi_eta, pxpxi * pxpeta, pxpxi * pypeta + pxpeta * pypxi,
                            pypxi * pypeta, pxpeta_eta, pypeta_eta,
                            pxpeta * pxpeta, 2 * pxpeta * pypeta, pypeta * pypeta;
                    MatrixXd Hessian_inv = Hessian.inverse();
                    SpVec Nx_x_xNy, NxNy_y_y;
                    Nx_x_xNy = Hessian_inv(2, 0) * Nxi_xiNeta + Hessian_inv(2, 1) * NxiNeta_eta +
                               Hessian_inv(2, 2) * Nxi_xi_xiNeta + Hessian_inv(2, 3) * Nxi_xiNeta_eta +
                               Hessian_inv(2, 4) * NxiNeta_eta_eta;
                    NxNy_y_y = Hessian_inv(4, 0) * Nxi_xiNeta + Hessian_inv(4, 1) * NxiNeta_eta +
                               Hessian_inv(4, 2) * Nxi_xi_xiNeta + Hessian_inv(4, 3) * Nxi_xiNeta_eta +
                               Hessian_inv(4, 4) * NxiNeta_eta_eta;
                    SpVec B = Nx_x_xNy + NxNy_y_y;
                    K_patch1 += weight[jj_x] * weight[jj_y] * Jacobian * B * B.transpose() * J_x * J_y;
                    F_patch1 += weight[jj_x] * weight[jj_y] * Jacobian * NxiNeta * force * J_x * J_y;
                }
            }
        }
    }

    for (int k = 0; k < K_patch1.outerSize(); ++k)
        for (SpMat::InnerIterator it(K_patch1, k); it; ++it) {
            coefficients.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
        }

    SpMat K_patch2(dof_patch2, dof_patch2);
    VectorXd F_patch2 = VectorXd::Zero(dof_patch2);
    MatrixXd constraints = (assemble_C0 * assemble_C1).block(dof_patch1, dof_patch1 - 2 * dof_y_patch1, dof_patch2,
                                                             dof_patch2 + 2 * dof_y_patch1 - 2 * dof_y_patch2);
    for (int ii_x = 0; ii_x < elements_x_patch2; ii_x++) {
        double J_x = (knots_x_patch2[ii_x + p_x + 1] - knots_x_patch2[ii_x + p_x]) / 2;
        double Middle_x = (knots_x_patch2[ii_x + p_x + 1] + knots_x_patch2[ii_x +
                                                                           p_x]) / 2;
        int i_x = Findspan(m_x_patch2, p_x, knots_x_patch2, Middle_x);
        for (int ii_y = 0; ii_y < elements_y_patch2; ii_y++) {
            double J_y = (knots_y_patch2[ii_y + p_y + 1] - knots_y_patch2[ii_y + p_y]) / 2;
            double Middle_y = (knots_y_patch2[ii_y + p_y + 1] + knots_y_patch2[ii_y +
                                                                               p_y]) / 2;
            int i_y = Findspan(m_y_patch2, p_y, knots_y_patch2, Middle_y);
            for (int jj_x = 0; jj_x < gaussian_points; jj_x++) {
                for (int jj_y = 0; jj_y < gaussian_points; jj_y++) {
                    double xi = Middle_x + J_x * gaussian[jj_x];
                    double eta = Middle_y + J_y * gaussian[jj_y];
                    double **ders_x, **ders_y;
                    DersBasisFuns(i_x, xi, p_x, knots_x_patch2, 2, ders_x);
                    DersBasisFuns(i_y, eta, p_y, knots_y_patch2, 2, ders_y);
                    SpVec Nxi(dof_x_patch2), Nxi_xi(dof_x_patch2), Nxi_xi_xi(dof_x_patch2), Neta(dof_y_patch2),
                            Neta_eta(dof_y_patch2), Neta_eta_eta(dof_y_patch2);
                    for (int kk_x = 0; kk_x < p_x + 1; kk_x++) {
                        Nxi.insert(i_x - p_x + kk_x) = ders_x[0][kk_x];
                        Nxi_xi.insert(i_x - p_x + kk_x) = ders_x[1][kk_x];
                        Nxi_xi_xi.insert(i_x - p_x + kk_x) = ders_x[2][kk_x];
                    }
                    for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
                        Neta.insert(i_y - p_y + kk_y) = ders_y[0][kk_y];
                        Neta_eta.insert(i_y - p_y + kk_y) = ders_y[1][kk_y];
                        Neta_eta_eta.insert(i_y - p_y + kk_y) = ders_y[2][kk_y];
                    }
                    for (int k = 0; k < 3; k++)
                        delete ders_x[k];
                    delete[] ders_x;
                    for (int k = 0; k < 3; k++)
                        delete ders_y[k];
                    delete[] ders_y;

                    double WNxi = (constraints.transpose() * Nxi).sum();
                    double WNxi_xi = (constraints.transpose() * Nxi_xi).sum();
                    double WNxi_xi_xi = (constraints.transpose() * Nxi_xi_xi).sum();
                    double WNeta = (constraints.transpose() * Neta).sum();
                    double WNeta_eta = (constraints.transpose() * Neta_eta).sum();
                    double WNeta_eta_eta = (constraints.transpose() * Neta_eta_eta).sum();
                    SpVec Rxi(dof_x_patch2), Rxi_xi(dof_x_patch2), Rxi_xi_xi(dof_x_patch2), Reta(dof_y_patch2),
                            Reta_eta(dof_y_patch2), Reta_eta_eta(dof_y_patch2);
                    Rxi = 1.0 / WNxi * Nxi;
                    Rxi_xi = 1.0 / WNxi * Nxi_xi - WNxi_xi / pow(WNxi, 2) * Nxi;
                    Rxi_xi_xi = -2 * WNxi_xi / pow(WNxi, 2) * Nxi_xi + 2 * pow(WNxi_xi, 2) / pow(WNxi, 3) * Nxi +
                                1.0 / WNxi * Nxi_xi_xi - WNxi_xi_xi / pow(WNxi, 2) * Nxi;
                    Reta = 1.0 / WNeta * Neta;
                    Reta_eta = 1.0 / WNeta * Neta_eta - WNeta_eta / pow(WNeta, 2) * Neta;
                    Reta_eta_eta =
                            -2 * WNeta_eta / pow(WNeta, 2) * Neta_eta + 2 * pow(WNeta_eta, 2) / pow(WNeta, 3) * Neta +
                            1.0 / WNeta * Neta_eta_eta - WNeta_eta_eta / pow(WNeta, 2) * Neta;

                    SpVec Nxi_xiNeta, NxiNeta_eta, Nxi_xi_xiNeta, NxiNeta_eta_eta, Nxi_xiNeta_eta, NxiNeta;
                    Nxi_xiNeta = (kroneckerProduct(Rxi_xi, Reta)).sparseView();
                    NxiNeta_eta = (kroneckerProduct(Rxi, Reta_eta)).sparseView();
                    Nxi_xi_xiNeta = (kroneckerProduct(Rxi_xi_xi, Reta)).sparseView();
                    NxiNeta_eta_eta = (kroneckerProduct(Rxi, Reta_eta_eta)).sparseView();
                    Nxi_xiNeta_eta = (kroneckerProduct(Rxi_xi, Reta_eta)).sparseView();
                    NxiNeta = (kroneckerProduct(Rxi, Reta)).sparseView();
                    double pxpxi, pxpeta, pypxi, pypeta, pxpxi_xi, pxpxi_eta, pxpeta_eta, pypxi_xi, pypxi_eta, pypeta_eta, x, y;
                    Geometry2(xi, eta, pxpxi, pxpeta, pypxi, pypeta, pxpxi_xi, pxpxi_eta, pxpeta_eta, pypxi_xi,
                              pypxi_eta, pypeta_eta, x,
                              y);
                    double force = forceTerm(x, y);
                    double Jacobian = pxpxi * pypeta - pxpeta * pypxi;
                    MatrixXd Hessian(5, 5);
                    Hessian << pxpxi, pypxi, 0, 0, 0, pxpeta, pypeta, 0, 0, 0, pxpxi_xi, pypxi_xi, pxpxi * pxpxi, 2 *
                                                                                                                  pxpxi *
                                                                                                                  pypxi,
                            pypxi * pypxi, pxpxi_eta, pypxi_eta, pxpxi * pxpeta, pxpxi * pypeta + pxpeta * pypxi,
                            pypxi * pypeta, pxpeta_eta, pypeta_eta,
                            pxpeta * pxpeta, 2 * pxpeta * pypeta, pypeta * pypeta;
                    MatrixXd Hessian_inv = Hessian.inverse();
                    SpVec Nx_x_xNy, NxNy_y_y;
                    Nx_x_xNy = Hessian_inv(2, 0) * Nxi_xiNeta + Hessian_inv(2, 1) * NxiNeta_eta +
                               Hessian_inv(2, 2) * Nxi_xi_xiNeta + Hessian_inv(2, 3) * Nxi_xiNeta_eta +
                               Hessian_inv(2, 4) * NxiNeta_eta_eta;
                    NxNy_y_y = Hessian_inv(4, 0) * Nxi_xiNeta + Hessian_inv(4, 1) * NxiNeta_eta +
                               Hessian_inv(4, 2) * Nxi_xi_xiNeta + Hessian_inv(4, 3) * Nxi_xiNeta_eta +
                               Hessian_inv(4, 4) * NxiNeta_eta_eta;
                    SpVec B = Nx_x_xNy + NxNy_y_y;
                    K_patch2 += weight[jj_x] * weight[jj_y] * Jacobian * B * B.transpose() * J_x * J_y;
                    F_patch2 += weight[jj_x] * weight[jj_y] * Jacobian * NxiNeta * force * J_x * J_y;
                }
            }
        }
    }
    for (int k = 0; k < K_patch2.outerSize(); ++k)
        for (SpMat::InnerIterator it(K_patch2, k); it; ++it) {
            coefficients.push_back(Eigen::Triplet<double>(it.row() + dof_patch1, it.col() + dof_patch1, it.value()));
        }
    SpMat K(dof_patch1 + dof_patch2, dof_patch1 + dof_patch2);
    K.setFromTriplets(coefficients.begin(), coefficients.end());
    coefficients.clear();
    VectorXd F(dof_patch1 + dof_patch2);
    F.segment(0, dof_patch1) = F_patch1;
    F.segment(dof_patch1, dof_patch2) = F_patch2;


    SpMat K_assembled = assemble_C1.transpose() * assemble_C0.transpose() * K * assemble_C0 * assemble_C1;
    VectorXd F_assembled = assemble_C1.transpose() * assemble_C0.transpose() * F;
    coefficients.clear();
    int x_it = 0, y_it = 0;
    for (int i = 0; i < dof_x_patch1; i++) {
        for (int j = 0; j < dof_y_patch1; j++) {
            if ((i == 0) || (i == 1) || (j == 0) || (j == 1) || (j == dof_y_patch1 - 2) || (j == dof_y_patch1 - 1)) {
                y_it++;
            } else {
                coefficients.push_back(Eigen::Triplet<double>(x_it, y_it, 1));
                x_it++;
                y_it++;
            }
        }
    }
    for (int i = 0; i < dof_x_patch2 - 2; i++) {
        for (int j = 0; j < dof_y_patch2; j++) {
            if ((i == dof_x_patch2 - 4) || (i == dof_x_patch2 - 3) || (j == 0) || (j == 1) || (j == dof_y_patch2 - 2)
                || (j == dof_y_patch2 - 1)) {
                y_it++;
            } else {
                coefficients.push_back(Eigen::Triplet<double>(x_it, y_it, 1));
                x_it++;
                y_it++;
            }
        }
    }
    SpMat transform(dof_patch1 + dof_patch2 - 2 * dof_y_patch2 - 2 * (dof_y_patch1 + dof_y_patch2) - 4 *
                                                                                                     (dof_x_patch1 +
                                                                                                      dof_x_patch2 -
                                                                                                      2) + 16,
                    dof_patch1 + dof_patch2 - 2 * dof_y_patch2);
    transform.setFromTriplets(coefficients.begin(), coefficients.end());

    SpMat K_cal = transform * K_assembled * transform.transpose();
    VectorXd F_cal = transform * F_assembled;


    ConjugateGradient<SparseMatrix<double>, Lower | Upper> cg;
    cg.compute(K_cal);
    VectorXd U_cal = cg.solve(F_cal);
    VectorXd U_patch1 = VectorXd::Zero(dof_patch1);
    x_it = 0;
    for (int i = 0; i < dof_x_patch1; i++) {
        for (int j = 0; j < dof_y_patch1; j++) {
            if ((i == 0) || (i == 1) || (j == 0) || (j == 1) || (j == dof_y_patch1 - 2) || (j == dof_y_patch1 - 1)) {
                U_patch1[j + i * dof_y_patch1] = 0;
            } else {
                U_patch1[j + i * dof_y_patch1] = U_cal[x_it];
                x_it++;
            }
        }
    }
    double L2_norm = 0, L2_norm_error = 0;
    for (int ii_x = 0; ii_x < elements_x_patch1; ii_x++) {
        double J_x = (knots_x_patch1[ii_x + p_x + 1] - knots_x_patch1[ii_x + p_x]) / 2;
        double Middle_x = (knots_x_patch1[ii_x + p_x + 1] + knots_x_patch1[ii_x +
                                                                           p_x]) / 2;
        int i_x = Findspan(m_x_patch1, p_x, knots_x_patch1, Middle_x);
        for (int ii_y = 0; ii_y < elements_y_patch1; ii_y++) {
            double J_y = (knots_y_patch1[ii_y + p_y + 1] - knots_y_patch1[ii_y + p_y]) / 2;
            double Middle_y = (knots_y_patch1[ii_y + p_y + 1] + knots_y_patch1[ii_y +
                                                                               p_y]) / 2;
            int i_y = Findspan(m_y_patch1, p_y, knots_y_patch1, Middle_y);
            for (int jj_x = 0; jj_x < gaussian_points; jj_x++) {
                for (int jj_y = 0; jj_y < gaussian_points; jj_y++) {
                    double xi = Middle_x + J_x * gaussian[jj_x];
                    double eta = Middle_y + J_y * gaussian[jj_y];
                    double **ders_x, **ders_y;
                    DersBasisFuns(i_x, xi, p_x, knots_x_patch1, 1, ders_x);
                    DersBasisFuns(i_y, eta, p_y, knots_y_patch1, 1, ders_y);
                    VectorXd Nxi = VectorXd::Zero(dof_x_patch1),
                            Neta = VectorXd::Zero(dof_y_patch1);
                    for (int kk_x = 0; kk_x < p_x + 1; kk_x++) {
                        Nxi(i_x - p_x + kk_x) = ders_x[0][kk_x];
                    }
                    for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
                        Neta(i_y - p_y + kk_y) = ders_y[0][kk_y];
                    }
                    for (int k = 0; k < 2; k++)
                        delete ders_x[k];
                    delete[] ders_x;
                    for (int k = 0; k < 2; k++)
                        delete ders_y[k];
                    delete[] ders_y;
                    VectorXd NxiNeta;
                    NxiNeta = kroneckerProduct(Nxi, Neta);
                    double pxpxi, pxpeta, pypxi, pypeta, pxpxi_xi, pxpxi_eta, pxpeta_eta, pypxi_xi, pypxi_eta, pypeta_eta, x, y;
                    Geometry1(xi, eta, pxpxi, pxpeta, pypxi, pypeta, pxpxi_xi, pxpxi_eta, pxpeta_eta, pypxi_xi,
                              pypxi_eta, pypeta_eta, x,
                              y);
                    double Jacobian = pxpxi * pypeta - pxpeta * pypxi;
                    MatrixXd U = (NxiNeta.transpose()) * U_patch1;
                    L2_norm_error += weight[jj_x] * weight[jj_y] * (pow(U(0, 0) - exactSolution(x,
                                                                                                y), 2)) * J_x * J_y
                                     * Jacobian;
                    L2_norm += weight[jj_x] * weight[jj_y] * (pow(exactSolution(x, y),
                                                                  2)) * J_x * J_y * Jacobian;
                }
            }
        }
    }
    cout << sqrt(L2_norm_error) / sqrt(L2_norm) << endl;
    return 0;
}

void GenerateKnot(int order, int refine, int insert, int repeat, double *insert_knot,
                  double *&knot, int &m) {
    vector<double> Knot;
    for (int i = 0; i <= order; i++) {
        Knot.push_back(0);
    }
    for (int i = 0; i < insert; i++) {
        for (int j = 0; j < repeat; j++) {
            Knot.push_back(insert_knot[i]);
        }
    }
    for (int i = 0; i <= order; i++) {
        Knot.push_back(1);
    }
    for (int i = 0; i < refine; i++) {
        for (int j = 0; j < Knot.end() - Knot.begin() - 1; j++) {
            double insertion;
            if (Knot[j] != Knot[j + 1]) {
                insertion = (Knot[j] + Knot[j + 1]) / 2;
                for (int k = 0; k < repeat; k++) {
                    Knot.insert(Knot.begin() + j + 1, insertion);
                    j++;
                }
            }
        }
    }
    m = Knot.end() - Knot.begin() - 1;
    knot = new double[m + 1];
    for (int i = 0; i < m + 1; i++)
        knot[i] = Knot[i];
}

void Geometry1(double xi, double eta, double &pxpxi, double &pxpeta, double &pypxi, double &pypeta,
               double &pxpxi_xi, double &pxpxi_eta, double &pxpeta_eta, double &pypxi_xi, double &pypxi_eta,
               double &pypeta_eta,
               double &x, double &y) {
    double knot_x[] = {0, 0, 1, 1};
    double knot_y[] = {0, 0, 1, 1};
    MatrixXd B_x(2, 2);
    MatrixXd B_y(2, 2);

    B_x << 0, 2, 0, 2;
    B_y << 0, 0, 4, 4;

    int p_x = 1, p_y = 1;
    int m_x = 3, m_y = 3;
    int dof_x = m_x - p_x, dof_y = m_y - p_y;
    int i_x = Findspan(m_x, p_x, knot_x, xi);
    int i_y = Findspan(m_y, p_y, knot_y, eta);
    double **ders_x, **ders_y;
    DersBasisFuns(i_x, xi, p_x, knot_x, 2, ders_x);
    DersBasisFuns(i_y, eta, p_y, knot_y, 2, ders_y);
    SpVec Nxi(dof_x), Nxi_xi(dof_x), Nxi_xi_xi(dof_x), Neta(dof_y), Neta_eta(dof_y), Neta_eta_eta(dof_y);
    for (int kk_x = 0; kk_x < p_x + 1; kk_x++) {
        Nxi.insert(i_x - p_x + kk_x) = ders_x[0][kk_x];
        Nxi_xi.insert(i_x - p_x + kk_x) = ders_x[1][kk_x];
        Nxi_xi_xi.insert(i_x - p_x + kk_x) = ders_x[2][kk_x];
    }
    for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
        Neta.insert(i_y - p_y + kk_y) = ders_y[0][kk_y];
        Neta_eta.insert(i_y - p_y + kk_y) = ders_y[1][kk_y];
        Neta_eta_eta.insert(i_y - p_y + kk_y) = ders_y[2][kk_y];
    }
    for (int k = 0; k < 3; k++)
        delete ders_x[k];
    delete[] ders_x;
    for (int k = 0; k < 3; k++)
        delete ders_y[k];
    delete[] ders_y;
    MatrixXd pxpxi_temp, pxpeta_temp, pypxi_temp, pypeta_temp, pxpxi_xi_temp, pxpxi_eta_temp, pxpeta_eta_temp,
            pypxi_xi_temp, pypxi_eta_temp, pypeta_eta_temp;
    pxpxi_temp = Neta.transpose() * B_x * Nxi_xi;
    pypxi_temp = Neta.transpose() * B_y * Nxi_xi;
    pxpeta_temp = Neta_eta.transpose() * B_x * Nxi;
    pypeta_temp = Neta_eta.transpose() * B_y * Nxi;
    pxpxi_xi_temp = Neta.transpose() * B_x * Nxi_xi_xi;
    pxpxi_eta_temp = Neta_eta.transpose() * B_x * Nxi_xi;
    pxpeta_eta_temp = Neta_eta_eta.transpose() * B_x * Nxi;
    pypxi_xi_temp = Neta.transpose() * B_y * Nxi_xi_xi;
    pypxi_eta_temp = Neta_eta.transpose() * B_y * Nxi_xi;
    pypeta_eta_temp = Neta_eta_eta.transpose() * B_y * Nxi;
    pxpxi = pxpxi_temp(0, 0);
    pxpeta = pxpeta_temp(0, 0);
    pypxi = pypxi_temp(0, 0);
    pypeta = pypeta_temp(0, 0);
    pxpxi_xi = pxpxi_xi_temp(0, 0);
    pxpxi_eta = pxpxi_eta_temp(0, 0);
    pxpeta_eta = pxpeta_eta_temp(0, 0);
    pypxi_xi = pypxi_xi_temp(0, 0);
    pypxi_eta = pypxi_eta_temp(0, 0);
    pypeta_eta = pypeta_eta_temp(0, 0);
    MatrixXd x0 = Neta.transpose() * B_x * Nxi;
    MatrixXd y0 = Neta.transpose() * B_y * Nxi;
    x = x0(0, 0);
    y = y0(0, 0);
}

void Geometry2(double xi, double eta, double &pxpxi, double &pxpeta, double &pypxi, double &pypeta,
               double &pxpxi_xi, double &pxpxi_eta, double &pxpeta_eta, double &pypxi_xi, double &pypxi_eta,
               double &pypeta_eta,
               double &x, double &y) {
    double knot_x[] = {0, 0, 1, 1};
    double knot_y[] = {0, 0, 1, 1};
    MatrixXd B_x(2, 2);
    MatrixXd B_y(2, 2);

    B_x << 2, 4, 2, 4;
    B_y << 0, 0, 4, 4;

    int p_x = 1, p_y = 1;
    int m_x = 3, m_y = 3;
    int dof_x = m_x - p_x, dof_y = m_y - p_y;
    int i_x = Findspan(m_x, p_x, knot_x, xi);
    int i_y = Findspan(m_y, p_y, knot_y, eta);
    double **ders_x, **ders_y;
    DersBasisFuns(i_x, xi, p_x, knot_x, 2, ders_x);
    DersBasisFuns(i_y, eta, p_y, knot_y, 2, ders_y);
    SpVec Nxi(dof_x), Nxi_xi(dof_x), Nxi_xi_xi(dof_x), Neta(dof_y), Neta_eta(dof_y), Neta_eta_eta(dof_y);
    for (int kk_x = 0; kk_x < p_x + 1; kk_x++) {
        Nxi.insert(i_x - p_x + kk_x) = ders_x[0][kk_x];
        Nxi_xi.insert(i_x - p_x + kk_x) = ders_x[1][kk_x];
        Nxi_xi_xi.insert(i_x - p_x + kk_x) = ders_x[2][kk_x];
    }
    for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
        Neta.insert(i_y - p_y + kk_y) = ders_y[0][kk_y];
        Neta_eta.insert(i_y - p_y + kk_y) = ders_y[1][kk_y];
        Neta_eta_eta.insert(i_y - p_y + kk_y) = ders_y[2][kk_y];
    }
    for (int k = 0; k < 3; k++)
        delete ders_x[k];
    delete[] ders_x;
    for (int k = 0; k < 3; k++)
        delete ders_y[k];
    delete[] ders_y;
    MatrixXd pxpxi_temp, pxpeta_temp, pypxi_temp, pypeta_temp, pxpxi_xi_temp, pxpxi_eta_temp, pxpeta_eta_temp,
            pypxi_xi_temp, pypxi_eta_temp, pypeta_eta_temp;
    pxpxi_temp = Neta.transpose() * B_x * Nxi_xi;
    pypxi_temp = Neta.transpose() * B_y * Nxi_xi;
    pxpeta_temp = Neta_eta.transpose() * B_x * Nxi;
    pypeta_temp = Neta_eta.transpose() * B_y * Nxi;
    pxpxi_xi_temp = Neta.transpose() * B_x * Nxi_xi_xi;
    pxpxi_eta_temp = Neta_eta.transpose() * B_x * Nxi_xi;
    pxpeta_eta_temp = Neta_eta_eta.transpose() * B_x * Nxi;
    pypxi_xi_temp = Neta.transpose() * B_y * Nxi_xi_xi;
    pypxi_eta_temp = Neta_eta.transpose() * B_y * Nxi_xi;
    pypeta_eta_temp = Neta_eta_eta.transpose() * B_y * Nxi;
    pxpxi = pxpxi_temp(0, 0);
    pxpeta = pxpeta_temp(0, 0);
    pypxi = pypxi_temp(0, 0);
    pypeta = pypeta_temp(0, 0);
    pxpxi_xi = pxpxi_xi_temp(0, 0);
    pxpxi_eta = pxpxi_eta_temp(0, 0);
    pxpeta_eta = pxpeta_eta_temp(0, 0);
    pypxi_xi = pypxi_xi_temp(0, 0);
    pypxi_eta = pypxi_eta_temp(0, 0);
    pypeta_eta = pypeta_eta_temp(0, 0);
    MatrixXd x0 = Neta.transpose() * B_x * Nxi;
    MatrixXd y0 = Neta.transpose() * B_y * Nxi;
    x = x0(0, 0);
    y = y0(0, 0);
}