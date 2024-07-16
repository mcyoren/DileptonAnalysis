#include "Run14AuAuLeptonCombyReco.h"

void Run14AuAuLeptonCombyReco::InitParams()
{
    BBC_VERTEX_CUT = 10;
    TRIGGERBIT = 4; // 4 - BBCLL1(>1 tubes) narrowvtx, 5 - BBCLL1(>1 tubes), 6 - BBCLL1(>1 tubes) novertex
    E_PT = 0.2;     // 0.15 in RDanalyzer for simulation to study smearing effect
    MAX_PT = 5;
    QUALITY[0] = 31;
    QUALITY[1] = 51;
    QUALITY[2] = 63;
    Z_GLOBAL = 75;
    DC_DEADMAP = 1; // by default, DC deadmap is on
    N0 = 2;
    DISP = 5;
    CHI2_NPE0 = 10;
    EOVERP = -999;
    DEP[0] = -2; // -2, can be applied offline
    DEP[1] =  5;  // 5
    PROB = 0.01; // 0.01
    TEMC = 5;
    EMCDPHI = 0.05; // ~5 sigma 0.05
    EMCDZ = 25;     // ~5 sigma 25
    EMCSDPHI = 5;   // remove sigmalized cut for Run14AuAu because the recalibrator doesn't work for simulation
    EMCSDZ = 5;
    RICH_GHOST = 10; // by default, RICH ghost cut is off

    std::cout << "Run14AuAuCut loaded... " << std::endl;
    std::cout << "***Event selection*** " << std::endl;
    std::cout <<  "BBC VERTEX CUT = "<<BBC_VERTEX_CUT<<std::endl;
    if (TRIGGERBIT > -99)
        std::cout << "trigger bit = " << TRIGGERBIT << std::endl;
    std::cout << "ZDC conincidence required" << std::endl;

    std::cout << "***Single track selection***" << std::endl;
    if (QUALITY[0] > -99)
        std::cout << "quality = " << QUALITY[0] << std::endl;
    if (QUALITY[1] > -99)
        std::cout << "quality = " << QUALITY[1] << std::endl;
    if (QUALITY[2] > -99)
        std::cout << "quality = " << QUALITY[2] << std::endl;
    if (Z_GLOBAL > -99)
        std::cout << "|z_global| < " << Z_GLOBAL << std::endl;
    if (DC_DEADMAP)
    {
        std::cout << "dc dead map in" << std::endl;
        init_rungroup();
        init_dead_area();
    }

    std::cout << "***e+/- track selection***" << std::endl;
    if (E_PT > -99)
        std::cout << "pT > " << E_PT << std::endl;
    if (MAX_PT > -99)
        std::cout << "pT < " << MAX_PT << std::endl;
    if (N0 > -99)
        std::cout << "n0 > " << N0 << std::endl;
    if (DISP > -99)
        std::cout << "disp < " << DISP << std::endl;
    if (CHI2_NPE0 > -99)
        std::cout << "chi2/npe0 < " << CHI2_NPE0 << std::endl;
    if (EMCDPHI > -99)
        std::cout << "emc matching cut in: emcdphi < " << EMCDPHI << "rad" << std::endl;
    if (EMCDZ > -99)
        std::cout << "emc matching cut in: emcdz < " << EMCDZ << "cm" << std::endl;
    if (EMCSDPHI > -99)
        std::cout << "emc matching cut (sigmalized) in: emcsdphi_e < " << EMCSDPHI << "sigma" << std::endl;
    if (EMCSDZ > -99)
        std::cout << "emc matching cut (sigmalized) in: emcsdz_e < " << EMCSDZ << "sigma" << std::endl;
    if (PROB > -99)
        std::cout << "prob > " << PROB << std::endl;
    if (TEMC > -99)
        std::cout << "TEMC < " << TEMC << std::endl;
    if (EOVERP > -99)
        std::cout << "E/p > " << EOVERP << std::endl;
    if (DEP[0] > -99)
        std::cout << "dep > " << DEP[0] << std::endl;
    if (DEP[1] > -99)
        std::cout << "dep < " << DEP[1] << std::endl;
    if (RICH_GHOST > -99)
    {
        std::cout << "RICH ghost cut in: decenter < " << RICH_GHOST << "cm" << std::endl;
    }

    if(!DC_DEADMAP)
    {
        std::cout << "setting init_dead_area to -9999" << std::endl;
        for (int i = 0; i < N_RUN_GRP; ++i)
        {
            for (int j = 0; j < N_SIDE; ++j)
            {
                for (int k = 0; k < N_ARM; ++k)
                {
                    for (int l = 0; l < MAX_DEAD_AREA; ++l)
                    {
                        dcmap_xx1[i][j][k][l] = -9999;
                        dcmap_yy1[i][j][k][l] = -9999;
                        dcmap_xx2[i][j][k][l] = -9999;
                        dcmap_yy2[i][j][k][l] = -9999;
                        dcmap_xx3[i][j][k][l] = -9999;
                        dcmap_yy3[i][j][k][l] = -9999;
                        dcmap_xx4[i][j][k][l] = -9999;
                        dcmap_yy4[i][j][k][l] = -9999;
                    }
                }
            }
        }
    }
}

float Run14AuAuLeptonCombyReco::get_board(float phi_dc)
{
    if (phi_dc > 1.57)
        return ((3.72402 - phi_dc + 0.008047 * cos(phi_dc + 0.87851)) / 0.01963496);
    return ((0.573231 + phi_dc - 0.0046 * cos(phi_dc + 0.05721)) / 0.01963496);
}

void Run14AuAuLeptonCombyReco::init_rungroup()
{
    for (int igrp = 0; igrp < N_RUN_GRP; ++igrp)
    {
        for (int irun = 0; irun < MAX_RUN; ++irun)
        {
            dcmap_runs[igrp][irun] = -9999;
        }
    }

    // run group 1
    int igrp = 1;
    for (int irun = 0; irun < N_dcmap_runs_rg1; ++irun)
    {
        dcmap_runs[igrp][irun] = dcmap_runs_rg1[irun];
    }

    // run group 2
    igrp = 2;
    for (int irun = 0; irun < N_dcmap_runs_rg2; ++irun)
    {
        dcmap_runs[igrp][irun] = dcmap_runs_rg2[irun];
    }

    // run group 3
    igrp = 3;
    for (int irun = 0; irun < N_dcmap_runs_rg3; ++irun)
    {
        dcmap_runs[igrp][irun] = dcmap_runs_rg3[irun];
    }

    // run group 4
    igrp = 4;
    for (int irun = 0; irun < N_dcmap_runs_rg4; ++irun)
    {
        dcmap_runs[igrp][irun] = dcmap_runs_rg4[irun];
    }

    // run group 5
    igrp = 5;
    for (int irun = 0; irun < N_dcmap_runs_rg5; ++irun)
    {
        dcmap_runs[igrp][irun] = dcmap_runs_rg5[irun];
    }

    // run group 6
    igrp = 6;
    for (int irun = 0; irun < N_dcmap_runs_rg6; ++irun)
    {
        dcmap_runs[igrp][irun] = dcmap_runs_rg6[irun];
    }

    // run group 7
    igrp = 7;
    for (int irun = 0; irun < N_dcmap_runs_rg7; ++irun)
    {
        dcmap_runs[igrp][irun] = dcmap_runs_rg7[irun];
    }

    // run group 8
    igrp = 8;
    for (int irun = 0; irun < N_dcmap_runs_rg8; ++irun)
    {
        dcmap_runs[igrp][irun] = dcmap_runs_rg8[irun];
    }
}

int Run14AuAuLeptonCombyReco::get_rungroup(int run_num)
{
    for (int igrp = 0; igrp < N_RUN_GRP; ++igrp)
    {
        for (int irun = 0; irun < MAX_RUN; ++irun)
        {
            if (run_num == dcmap_runs[igrp][irun])
                return igrp;
        }
    }
    return 0; // no match (bad runs)
}

void Run14AuAuLeptonCombyReco::init_dead_area()
{

    // run group 1
    int igrp = 1;
    for (int iarea = 0; iarea < (int)(N_ES_rg1 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][0][iarea] = dcmap_x_ES_rg1[4 * iarea];
        dcmap_yy1[igrp][0][0][iarea] = dcmap_y_ES_rg1[4 * iarea];
        dcmap_xx2[igrp][0][0][iarea] = dcmap_x_ES_rg1[4 * iarea + 1];
        dcmap_yy2[igrp][0][0][iarea] = dcmap_y_ES_rg1[4 * iarea + 1];
        dcmap_xx3[igrp][0][0][iarea] = dcmap_x_ES_rg1[4 * iarea + 2];
        dcmap_yy3[igrp][0][0][iarea] = dcmap_y_ES_rg1[4 * iarea + 2];
        dcmap_xx4[igrp][0][0][iarea] = dcmap_x_ES_rg1[4 * iarea + 3];
        dcmap_yy4[igrp][0][0][iarea] = dcmap_y_ES_rg1[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WS_rg1 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][1][iarea] = dcmap_x_WS_rg1[4 * iarea];
        dcmap_yy1[igrp][0][1][iarea] = dcmap_y_WS_rg1[4 * iarea];
        dcmap_xx2[igrp][0][1][iarea] = dcmap_x_WS_rg1[4 * iarea + 1];
        dcmap_yy2[igrp][0][1][iarea] = dcmap_y_WS_rg1[4 * iarea + 1];
        dcmap_xx3[igrp][0][1][iarea] = dcmap_x_WS_rg1[4 * iarea + 2];
        dcmap_yy3[igrp][0][1][iarea] = dcmap_y_WS_rg1[4 * iarea + 2];
        dcmap_xx4[igrp][0][1][iarea] = dcmap_x_WS_rg1[4 * iarea + 3];
        dcmap_yy4[igrp][0][1][iarea] = dcmap_y_WS_rg1[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_EN_rg1 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][0][iarea] = dcmap_x_EN_rg1[4 * iarea];
        dcmap_yy1[igrp][1][0][iarea] = dcmap_y_EN_rg1[4 * iarea];
        dcmap_xx2[igrp][1][0][iarea] = dcmap_x_EN_rg1[4 * iarea + 1];
        dcmap_yy2[igrp][1][0][iarea] = dcmap_y_EN_rg1[4 * iarea + 1];
        dcmap_xx3[igrp][1][0][iarea] = dcmap_x_EN_rg1[4 * iarea + 2];
        dcmap_yy3[igrp][1][0][iarea] = dcmap_y_EN_rg1[4 * iarea + 2];
        dcmap_xx4[igrp][1][0][iarea] = dcmap_x_EN_rg1[4 * iarea + 3];
        dcmap_yy4[igrp][1][0][iarea] = dcmap_y_EN_rg1[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WN_rg1 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][1][iarea] = dcmap_x_WN_rg1[4 * iarea];
        dcmap_yy1[igrp][1][1][iarea] = dcmap_y_WN_rg1[4 * iarea];
        dcmap_xx2[igrp][1][1][iarea] = dcmap_x_WN_rg1[4 * iarea + 1];
        dcmap_yy2[igrp][1][1][iarea] = dcmap_y_WN_rg1[4 * iarea + 1];
        dcmap_xx3[igrp][1][1][iarea] = dcmap_x_WN_rg1[4 * iarea + 2];
        dcmap_yy3[igrp][1][1][iarea] = dcmap_y_WN_rg1[4 * iarea + 2];
        dcmap_xx4[igrp][1][1][iarea] = dcmap_x_WN_rg1[4 * iarea + 3];
        dcmap_yy4[igrp][1][1][iarea] = dcmap_y_WN_rg1[4 * iarea + 3];
    }

    // run group 2
    igrp = 2;
    for (int iarea = 0; iarea < (int)(N_ES_rg2 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][0][iarea] = dcmap_x_ES_rg2[4 * iarea];
        dcmap_yy1[igrp][0][0][iarea] = dcmap_y_ES_rg2[4 * iarea];
        dcmap_xx2[igrp][0][0][iarea] = dcmap_x_ES_rg2[4 * iarea + 1];
        dcmap_yy2[igrp][0][0][iarea] = dcmap_y_ES_rg2[4 * iarea + 1];
        dcmap_xx3[igrp][0][0][iarea] = dcmap_x_ES_rg2[4 * iarea + 2];
        dcmap_yy3[igrp][0][0][iarea] = dcmap_y_ES_rg2[4 * iarea + 2];
        dcmap_xx4[igrp][0][0][iarea] = dcmap_x_ES_rg2[4 * iarea + 3];
        dcmap_yy4[igrp][0][0][iarea] = dcmap_y_ES_rg2[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WS_rg2 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][1][iarea] = dcmap_x_WS_rg2[4 * iarea];
        dcmap_yy1[igrp][0][1][iarea] = dcmap_y_WS_rg2[4 * iarea];
        dcmap_xx2[igrp][0][1][iarea] = dcmap_x_WS_rg2[4 * iarea + 1];
        dcmap_yy2[igrp][0][1][iarea] = dcmap_y_WS_rg2[4 * iarea + 1];
        dcmap_xx3[igrp][0][1][iarea] = dcmap_x_WS_rg2[4 * iarea + 2];
        dcmap_yy3[igrp][0][1][iarea] = dcmap_y_WS_rg2[4 * iarea + 2];
        dcmap_xx4[igrp][0][1][iarea] = dcmap_x_WS_rg2[4 * iarea + 3];
        dcmap_yy4[igrp][0][1][iarea] = dcmap_y_WS_rg2[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_EN_rg2 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][0][iarea] = dcmap_x_EN_rg2[4 * iarea];
        dcmap_yy1[igrp][1][0][iarea] = dcmap_y_EN_rg2[4 * iarea];
        dcmap_xx2[igrp][1][0][iarea] = dcmap_x_EN_rg2[4 * iarea + 1];
        dcmap_yy2[igrp][1][0][iarea] = dcmap_y_EN_rg2[4 * iarea + 1];
        dcmap_xx3[igrp][1][0][iarea] = dcmap_x_EN_rg2[4 * iarea + 2];
        dcmap_yy3[igrp][1][0][iarea] = dcmap_y_EN_rg2[4 * iarea + 2];
        dcmap_xx4[igrp][1][0][iarea] = dcmap_x_EN_rg2[4 * iarea + 3];
        dcmap_yy4[igrp][1][0][iarea] = dcmap_y_EN_rg2[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WN_rg2 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][1][iarea] = dcmap_x_WN_rg2[4 * iarea];
        dcmap_yy1[igrp][1][1][iarea] = dcmap_y_WN_rg2[4 * iarea];
        dcmap_xx2[igrp][1][1][iarea] = dcmap_x_WN_rg2[4 * iarea + 1];
        dcmap_yy2[igrp][1][1][iarea] = dcmap_y_WN_rg2[4 * iarea + 1];
        dcmap_xx3[igrp][1][1][iarea] = dcmap_x_WN_rg2[4 * iarea + 2];
        dcmap_yy3[igrp][1][1][iarea] = dcmap_y_WN_rg2[4 * iarea + 2];
        dcmap_xx4[igrp][1][1][iarea] = dcmap_x_WN_rg2[4 * iarea + 3];
        dcmap_yy4[igrp][1][1][iarea] = dcmap_y_WN_rg2[4 * iarea + 3];
    }

    // run group 3
    igrp = 3;
    for (int iarea = 0; iarea < (int)(N_ES_rg3 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][0][iarea] = dcmap_x_ES_rg3[4 * iarea];
        dcmap_yy1[igrp][0][0][iarea] = dcmap_y_ES_rg3[4 * iarea];
        dcmap_xx2[igrp][0][0][iarea] = dcmap_x_ES_rg3[4 * iarea + 1];
        dcmap_yy2[igrp][0][0][iarea] = dcmap_y_ES_rg3[4 * iarea + 1];
        dcmap_xx3[igrp][0][0][iarea] = dcmap_x_ES_rg3[4 * iarea + 2];
        dcmap_yy3[igrp][0][0][iarea] = dcmap_y_ES_rg3[4 * iarea + 2];
        dcmap_xx4[igrp][0][0][iarea] = dcmap_x_ES_rg3[4 * iarea + 3];
        dcmap_yy4[igrp][0][0][iarea] = dcmap_y_ES_rg3[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WS_rg3 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][1][iarea] = dcmap_x_WS_rg3[4 * iarea];
        dcmap_yy1[igrp][0][1][iarea] = dcmap_y_WS_rg3[4 * iarea];
        dcmap_xx2[igrp][0][1][iarea] = dcmap_x_WS_rg3[4 * iarea + 1];
        dcmap_yy2[igrp][0][1][iarea] = dcmap_y_WS_rg3[4 * iarea + 1];
        dcmap_xx3[igrp][0][1][iarea] = dcmap_x_WS_rg3[4 * iarea + 2];
        dcmap_yy3[igrp][0][1][iarea] = dcmap_y_WS_rg3[4 * iarea + 2];
        dcmap_xx4[igrp][0][1][iarea] = dcmap_x_WS_rg3[4 * iarea + 3];
        dcmap_yy4[igrp][0][1][iarea] = dcmap_y_WS_rg3[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_EN_rg3 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][0][iarea] = dcmap_x_EN_rg3[4 * iarea];
        dcmap_yy1[igrp][1][0][iarea] = dcmap_y_EN_rg3[4 * iarea];
        dcmap_xx2[igrp][1][0][iarea] = dcmap_x_EN_rg3[4 * iarea + 1];
        dcmap_yy2[igrp][1][0][iarea] = dcmap_y_EN_rg3[4 * iarea + 1];
        dcmap_xx3[igrp][1][0][iarea] = dcmap_x_EN_rg3[4 * iarea + 2];
        dcmap_yy3[igrp][1][0][iarea] = dcmap_y_EN_rg3[4 * iarea + 2];
        dcmap_xx4[igrp][1][0][iarea] = dcmap_x_EN_rg3[4 * iarea + 3];
        dcmap_yy4[igrp][1][0][iarea] = dcmap_y_EN_rg3[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WN_rg3 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][1][iarea] = dcmap_x_WN_rg3[4 * iarea];
        dcmap_yy1[igrp][1][1][iarea] = dcmap_y_WN_rg3[4 * iarea];
        dcmap_xx2[igrp][1][1][iarea] = dcmap_x_WN_rg3[4 * iarea + 1];
        dcmap_yy2[igrp][1][1][iarea] = dcmap_y_WN_rg3[4 * iarea + 1];
        dcmap_xx3[igrp][1][1][iarea] = dcmap_x_WN_rg3[4 * iarea + 2];
        dcmap_yy3[igrp][1][1][iarea] = dcmap_y_WN_rg3[4 * iarea + 2];
        dcmap_xx4[igrp][1][1][iarea] = dcmap_x_WN_rg3[4 * iarea + 3];
        dcmap_yy4[igrp][1][1][iarea] = dcmap_y_WN_rg3[4 * iarea + 3];
    }

    // run group 4
    igrp = 4;
    for (int iarea = 0; iarea < (int)(N_ES_rg4 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][0][iarea] = dcmap_x_ES_rg4[4 * iarea];
        dcmap_yy1[igrp][0][0][iarea] = dcmap_y_ES_rg4[4 * iarea];
        dcmap_xx2[igrp][0][0][iarea] = dcmap_x_ES_rg4[4 * iarea + 1];
        dcmap_yy2[igrp][0][0][iarea] = dcmap_y_ES_rg4[4 * iarea + 1];
        dcmap_xx3[igrp][0][0][iarea] = dcmap_x_ES_rg4[4 * iarea + 2];
        dcmap_yy3[igrp][0][0][iarea] = dcmap_y_ES_rg4[4 * iarea + 2];
        dcmap_xx4[igrp][0][0][iarea] = dcmap_x_ES_rg4[4 * iarea + 3];
        dcmap_yy4[igrp][0][0][iarea] = dcmap_y_ES_rg4[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WS_rg4 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][1][iarea] = dcmap_x_WS_rg4[4 * iarea];
        dcmap_yy1[igrp][0][1][iarea] = dcmap_y_WS_rg4[4 * iarea];
        dcmap_xx2[igrp][0][1][iarea] = dcmap_x_WS_rg4[4 * iarea + 1];
        dcmap_yy2[igrp][0][1][iarea] = dcmap_y_WS_rg4[4 * iarea + 1];
        dcmap_xx3[igrp][0][1][iarea] = dcmap_x_WS_rg4[4 * iarea + 2];
        dcmap_yy3[igrp][0][1][iarea] = dcmap_y_WS_rg4[4 * iarea + 2];
        dcmap_xx4[igrp][0][1][iarea] = dcmap_x_WS_rg4[4 * iarea + 3];
        dcmap_yy4[igrp][0][1][iarea] = dcmap_y_WS_rg4[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_EN_rg4 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][0][iarea] = dcmap_x_EN_rg4[4 * iarea];
        dcmap_yy1[igrp][1][0][iarea] = dcmap_y_EN_rg4[4 * iarea];
        dcmap_xx2[igrp][1][0][iarea] = dcmap_x_EN_rg4[4 * iarea + 1];
        dcmap_yy2[igrp][1][0][iarea] = dcmap_y_EN_rg4[4 * iarea + 1];
        dcmap_xx3[igrp][1][0][iarea] = dcmap_x_EN_rg4[4 * iarea + 2];
        dcmap_yy3[igrp][1][0][iarea] = dcmap_y_EN_rg4[4 * iarea + 2];
        dcmap_xx4[igrp][1][0][iarea] = dcmap_x_EN_rg4[4 * iarea + 3];
        dcmap_yy4[igrp][1][0][iarea] = dcmap_y_EN_rg4[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WN_rg4 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][1][iarea] = dcmap_x_WN_rg4[4 * iarea];
        dcmap_yy1[igrp][1][1][iarea] = dcmap_y_WN_rg4[4 * iarea];
        dcmap_xx2[igrp][1][1][iarea] = dcmap_x_WN_rg4[4 * iarea + 1];
        dcmap_yy2[igrp][1][1][iarea] = dcmap_y_WN_rg4[4 * iarea + 1];
        dcmap_xx3[igrp][1][1][iarea] = dcmap_x_WN_rg4[4 * iarea + 2];
        dcmap_yy3[igrp][1][1][iarea] = dcmap_y_WN_rg4[4 * iarea + 2];
        dcmap_xx4[igrp][1][1][iarea] = dcmap_x_WN_rg4[4 * iarea + 3];
        dcmap_yy4[igrp][1][1][iarea] = dcmap_y_WN_rg4[4 * iarea + 3];
    }

    // run group 5
    igrp = 5;
    for (int iarea = 0; iarea < (int)(N_ES_rg5 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][0][iarea] = dcmap_x_ES_rg5[4 * iarea];
        dcmap_yy1[igrp][0][0][iarea] = dcmap_y_ES_rg5[4 * iarea];
        dcmap_xx2[igrp][0][0][iarea] = dcmap_x_ES_rg5[4 * iarea + 1];
        dcmap_yy2[igrp][0][0][iarea] = dcmap_y_ES_rg5[4 * iarea + 1];
        dcmap_xx3[igrp][0][0][iarea] = dcmap_x_ES_rg5[4 * iarea + 2];
        dcmap_yy3[igrp][0][0][iarea] = dcmap_y_ES_rg5[4 * iarea + 2];
        dcmap_xx4[igrp][0][0][iarea] = dcmap_x_ES_rg5[4 * iarea + 3];
        dcmap_yy4[igrp][0][0][iarea] = dcmap_y_ES_rg5[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WS_rg5 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][1][iarea] = dcmap_x_WS_rg5[4 * iarea];
        dcmap_yy1[igrp][0][1][iarea] = dcmap_y_WS_rg5[4 * iarea];
        dcmap_xx2[igrp][0][1][iarea] = dcmap_x_WS_rg5[4 * iarea + 1];
        dcmap_yy2[igrp][0][1][iarea] = dcmap_y_WS_rg5[4 * iarea + 1];
        dcmap_xx3[igrp][0][1][iarea] = dcmap_x_WS_rg5[4 * iarea + 2];
        dcmap_yy3[igrp][0][1][iarea] = dcmap_y_WS_rg5[4 * iarea + 2];
        dcmap_xx4[igrp][0][1][iarea] = dcmap_x_WS_rg5[4 * iarea + 3];
        dcmap_yy4[igrp][0][1][iarea] = dcmap_y_WS_rg5[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_EN_rg5 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][0][iarea] = dcmap_x_EN_rg5[4 * iarea];
        dcmap_yy1[igrp][1][0][iarea] = dcmap_y_EN_rg5[4 * iarea];
        dcmap_xx2[igrp][1][0][iarea] = dcmap_x_EN_rg5[4 * iarea + 1];
        dcmap_yy2[igrp][1][0][iarea] = dcmap_y_EN_rg5[4 * iarea + 1];
        dcmap_xx3[igrp][1][0][iarea] = dcmap_x_EN_rg5[4 * iarea + 2];
        dcmap_yy3[igrp][1][0][iarea] = dcmap_y_EN_rg5[4 * iarea + 2];
        dcmap_xx4[igrp][1][0][iarea] = dcmap_x_EN_rg5[4 * iarea + 3];
        dcmap_yy4[igrp][1][0][iarea] = dcmap_y_EN_rg5[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WN_rg5 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][1][iarea] = dcmap_x_WN_rg5[4 * iarea];
        dcmap_yy1[igrp][1][1][iarea] = dcmap_y_WN_rg5[4 * iarea];
        dcmap_xx2[igrp][1][1][iarea] = dcmap_x_WN_rg5[4 * iarea + 1];
        dcmap_yy2[igrp][1][1][iarea] = dcmap_y_WN_rg5[4 * iarea + 1];
        dcmap_xx3[igrp][1][1][iarea] = dcmap_x_WN_rg5[4 * iarea + 2];
        dcmap_yy3[igrp][1][1][iarea] = dcmap_y_WN_rg5[4 * iarea + 2];
        dcmap_xx4[igrp][1][1][iarea] = dcmap_x_WN_rg5[4 * iarea + 3];
        dcmap_yy4[igrp][1][1][iarea] = dcmap_y_WN_rg5[4 * iarea + 3];
    }

    // run group 6
    igrp = 6;
    for (int iarea = 0; iarea < (int)(N_ES_rg6 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][0][iarea] = dcmap_x_ES_rg6[4 * iarea];
        dcmap_yy1[igrp][0][0][iarea] = dcmap_y_ES_rg6[4 * iarea];
        dcmap_xx2[igrp][0][0][iarea] = dcmap_x_ES_rg6[4 * iarea + 1];
        dcmap_yy2[igrp][0][0][iarea] = dcmap_y_ES_rg6[4 * iarea + 1];
        dcmap_xx3[igrp][0][0][iarea] = dcmap_x_ES_rg6[4 * iarea + 2];
        dcmap_yy3[igrp][0][0][iarea] = dcmap_y_ES_rg6[4 * iarea + 2];
        dcmap_xx4[igrp][0][0][iarea] = dcmap_x_ES_rg6[4 * iarea + 3];
        dcmap_yy4[igrp][0][0][iarea] = dcmap_y_ES_rg6[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WS_rg6 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][1][iarea] = dcmap_x_WS_rg6[4 * iarea];
        dcmap_yy1[igrp][0][1][iarea] = dcmap_y_WS_rg6[4 * iarea];
        dcmap_xx2[igrp][0][1][iarea] = dcmap_x_WS_rg6[4 * iarea + 1];
        dcmap_yy2[igrp][0][1][iarea] = dcmap_y_WS_rg6[4 * iarea + 1];
        dcmap_xx3[igrp][0][1][iarea] = dcmap_x_WS_rg6[4 * iarea + 2];
        dcmap_yy3[igrp][0][1][iarea] = dcmap_y_WS_rg6[4 * iarea + 2];
        dcmap_xx4[igrp][0][1][iarea] = dcmap_x_WS_rg6[4 * iarea + 3];
        dcmap_yy4[igrp][0][1][iarea] = dcmap_y_WS_rg6[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_EN_rg6 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][0][iarea] = dcmap_x_EN_rg6[4 * iarea];
        dcmap_yy1[igrp][1][0][iarea] = dcmap_y_EN_rg6[4 * iarea];
        dcmap_xx2[igrp][1][0][iarea] = dcmap_x_EN_rg6[4 * iarea + 1];
        dcmap_yy2[igrp][1][0][iarea] = dcmap_y_EN_rg6[4 * iarea + 1];
        dcmap_xx3[igrp][1][0][iarea] = dcmap_x_EN_rg6[4 * iarea + 2];
        dcmap_yy3[igrp][1][0][iarea] = dcmap_y_EN_rg6[4 * iarea + 2];
        dcmap_xx4[igrp][1][0][iarea] = dcmap_x_EN_rg6[4 * iarea + 3];
        dcmap_yy4[igrp][1][0][iarea] = dcmap_y_EN_rg6[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WN_rg6 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][1][iarea] = dcmap_x_WN_rg6[4 * iarea];
        dcmap_yy1[igrp][1][1][iarea] = dcmap_y_WN_rg6[4 * iarea];
        dcmap_xx2[igrp][1][1][iarea] = dcmap_x_WN_rg6[4 * iarea + 1];
        dcmap_yy2[igrp][1][1][iarea] = dcmap_y_WN_rg6[4 * iarea + 1];
        dcmap_xx3[igrp][1][1][iarea] = dcmap_x_WN_rg6[4 * iarea + 2];
        dcmap_yy3[igrp][1][1][iarea] = dcmap_y_WN_rg6[4 * iarea + 2];
        dcmap_xx4[igrp][1][1][iarea] = dcmap_x_WN_rg6[4 * iarea + 3];
        dcmap_yy4[igrp][1][1][iarea] = dcmap_y_WN_rg6[4 * iarea + 3];
    }

    // run group 7
    igrp = 7;
    for (int iarea = 0; iarea < (int)(N_ES_rg7 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][0][iarea] = dcmap_x_ES_rg7[4 * iarea];
        dcmap_yy1[igrp][0][0][iarea] = dcmap_y_ES_rg7[4 * iarea];
        dcmap_xx2[igrp][0][0][iarea] = dcmap_x_ES_rg7[4 * iarea + 1];
        dcmap_yy2[igrp][0][0][iarea] = dcmap_y_ES_rg7[4 * iarea + 1];
        dcmap_xx3[igrp][0][0][iarea] = dcmap_x_ES_rg7[4 * iarea + 2];
        dcmap_yy3[igrp][0][0][iarea] = dcmap_y_ES_rg7[4 * iarea + 2];
        dcmap_xx4[igrp][0][0][iarea] = dcmap_x_ES_rg7[4 * iarea + 3];
        dcmap_yy4[igrp][0][0][iarea] = dcmap_y_ES_rg7[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WS_rg7 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][1][iarea] = dcmap_x_WS_rg7[4 * iarea];
        dcmap_yy1[igrp][0][1][iarea] = dcmap_y_WS_rg7[4 * iarea];
        dcmap_xx2[igrp][0][1][iarea] = dcmap_x_WS_rg7[4 * iarea + 1];
        dcmap_yy2[igrp][0][1][iarea] = dcmap_y_WS_rg7[4 * iarea + 1];
        dcmap_xx3[igrp][0][1][iarea] = dcmap_x_WS_rg7[4 * iarea + 2];
        dcmap_yy3[igrp][0][1][iarea] = dcmap_y_WS_rg7[4 * iarea + 2];
        dcmap_xx4[igrp][0][1][iarea] = dcmap_x_WS_rg7[4 * iarea + 3];
        dcmap_yy4[igrp][0][1][iarea] = dcmap_y_WS_rg7[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_EN_rg7 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][0][iarea] = dcmap_x_EN_rg7[4 * iarea];
        dcmap_yy1[igrp][1][0][iarea] = dcmap_y_EN_rg7[4 * iarea];
        dcmap_xx2[igrp][1][0][iarea] = dcmap_x_EN_rg7[4 * iarea + 1];
        dcmap_yy2[igrp][1][0][iarea] = dcmap_y_EN_rg7[4 * iarea + 1];
        dcmap_xx3[igrp][1][0][iarea] = dcmap_x_EN_rg7[4 * iarea + 2];
        dcmap_yy3[igrp][1][0][iarea] = dcmap_y_EN_rg7[4 * iarea + 2];
        dcmap_xx4[igrp][1][0][iarea] = dcmap_x_EN_rg7[4 * iarea + 3];
        dcmap_yy4[igrp][1][0][iarea] = dcmap_y_EN_rg7[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WN_rg7 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][1][iarea] = dcmap_x_WN_rg7[4 * iarea];
        dcmap_yy1[igrp][1][1][iarea] = dcmap_y_WN_rg7[4 * iarea];
        dcmap_xx2[igrp][1][1][iarea] = dcmap_x_WN_rg7[4 * iarea + 1];
        dcmap_yy2[igrp][1][1][iarea] = dcmap_y_WN_rg7[4 * iarea + 1];
        dcmap_xx3[igrp][1][1][iarea] = dcmap_x_WN_rg7[4 * iarea + 2];
        dcmap_yy3[igrp][1][1][iarea] = dcmap_y_WN_rg7[4 * iarea + 2];
        dcmap_xx4[igrp][1][1][iarea] = dcmap_x_WN_rg7[4 * iarea + 3];
        dcmap_yy4[igrp][1][1][iarea] = dcmap_y_WN_rg7[4 * iarea + 3];
    }

    // run group 8
    igrp = 8;
    for (int iarea = 0; iarea < (int)(N_ES_rg8 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][0][iarea] = dcmap_x_ES_rg8[4 * iarea];
        dcmap_yy1[igrp][0][0][iarea] = dcmap_y_ES_rg8[4 * iarea];
        dcmap_xx2[igrp][0][0][iarea] = dcmap_x_ES_rg8[4 * iarea + 1];
        dcmap_yy2[igrp][0][0][iarea] = dcmap_y_ES_rg8[4 * iarea + 1];
        dcmap_xx3[igrp][0][0][iarea] = dcmap_x_ES_rg8[4 * iarea + 2];
        dcmap_yy3[igrp][0][0][iarea] = dcmap_y_ES_rg8[4 * iarea + 2];
        dcmap_xx4[igrp][0][0][iarea] = dcmap_x_ES_rg8[4 * iarea + 3];
        dcmap_yy4[igrp][0][0][iarea] = dcmap_y_ES_rg8[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WS_rg8 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][1][iarea] = dcmap_x_WS_rg8[4 * iarea];
        dcmap_yy1[igrp][0][1][iarea] = dcmap_y_WS_rg8[4 * iarea];
        dcmap_xx2[igrp][0][1][iarea] = dcmap_x_WS_rg8[4 * iarea + 1];
        dcmap_yy2[igrp][0][1][iarea] = dcmap_y_WS_rg8[4 * iarea + 1];
        dcmap_xx3[igrp][0][1][iarea] = dcmap_x_WS_rg8[4 * iarea + 2];
        dcmap_yy3[igrp][0][1][iarea] = dcmap_y_WS_rg8[4 * iarea + 2];
        dcmap_xx4[igrp][0][1][iarea] = dcmap_x_WS_rg8[4 * iarea + 3];
        dcmap_yy4[igrp][0][1][iarea] = dcmap_y_WS_rg8[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_EN_rg8 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][0][iarea] = dcmap_x_EN_rg8[4 * iarea];
        dcmap_yy1[igrp][1][0][iarea] = dcmap_y_EN_rg8[4 * iarea];
        dcmap_xx2[igrp][1][0][iarea] = dcmap_x_EN_rg8[4 * iarea + 1];
        dcmap_yy2[igrp][1][0][iarea] = dcmap_y_EN_rg8[4 * iarea + 1];
        dcmap_xx3[igrp][1][0][iarea] = dcmap_x_EN_rg8[4 * iarea + 2];
        dcmap_yy3[igrp][1][0][iarea] = dcmap_y_EN_rg8[4 * iarea + 2];
        dcmap_xx4[igrp][1][0][iarea] = dcmap_x_EN_rg8[4 * iarea + 3];
        dcmap_yy4[igrp][1][0][iarea] = dcmap_y_EN_rg8[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WN_rg8 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][1][iarea] = dcmap_x_WN_rg8[4 * iarea];
        dcmap_yy1[igrp][1][1][iarea] = dcmap_y_WN_rg8[4 * iarea];
        dcmap_xx2[igrp][1][1][iarea] = dcmap_x_WN_rg8[4 * iarea + 1];
        dcmap_yy2[igrp][1][1][iarea] = dcmap_y_WN_rg8[4 * iarea + 1];
        dcmap_xx3[igrp][1][1][iarea] = dcmap_x_WN_rg8[4 * iarea + 2];
        dcmap_yy3[igrp][1][1][iarea] = dcmap_y_WN_rg8[4 * iarea + 2];
        dcmap_xx4[igrp][1][1][iarea] = dcmap_x_WN_rg8[4 * iarea + 3];
        dcmap_yy4[igrp][1][1][iarea] = dcmap_y_WN_rg8[4 * iarea + 3];
    }
}

bool Run14AuAuLeptonCombyReco::left_bound(float x, float y, float xx1, float yy1, float xx2, float yy2)
{ // bound line defined by (xx1, yy1), (xx2, yy2), check if (x, y) is on the right-hand side of the bound line
    if (xx1 == xx2 && x > xx1)
        return true;
    else
    {
        const float k = (yy2 - yy1) / (xx2 - xx1); // slope
        if (k > 0 && y < k * (x - xx2) + yy2)
            return true;
        if (k < 0 && y > k * (x - xx2) + yy2)
            return true;
    }
    return false;
}

bool Run14AuAuLeptonCombyReco::right_bound(float x, float y, float xx1, float yy1, float xx2, float yy2)
{ // need to also consider slope infinity case (xx1==xx2)
    if (xx1 == xx2 && x < xx1)
        return true;
    else
    {
        const float k = (yy2 - yy1) / (xx2 - xx1); // slope
        if (k > 0 && y > k * (x - xx2) + yy2)
            return true;
        if (k < 0 && y < k * (x - xx2) + yy2)
            return true;
    }
    return false;
}

bool Run14AuAuLeptonCombyReco::dead_region(float x, float y, float xx1, float yy1, float xx2, float yy2, float xx3, float yy3, float xx4, float yy4)
{ // if in dead region, return true
    // alpha cut is equivalent to pt, so the middle region better be kept to keep the large pt
    if (xx1 < -999)
        return false; // not a dead region
    if (xx1 > -999 && xx1 < 0 && right_bound(x, y, xx3, yy3, xx4, yy4))
        return true; // only right bound
    if (xx3 > -999 && xx3 < 0 && left_bound(x, y, xx1, yy1, xx2, yy2))
        return true; // only left bound
    if (left_bound(x, y, xx1, yy1, xx2, yy2) && right_bound(x, y, xx3, yy3, xx4, yy4))
        return true;
    return false;
}

void Run14AuAuLeptonCombyReco::MoonWalk()
{
    fCDH = nullptr;
    se = nullptr;
    for (int itow = 0; itow < 24768; itow++)
    {
        Walk[itow] = 0;
        Walk2[itow] = 0;
        Walk3[itow] = 0;
        Walk4[itow] = 0;
        Walk5[itow] = 0;
        T0Offset[itow] = 0;
        T0OffsetSigma[itow] = 0;
    }
    for (int isec = 0; isec < 8; isec++)
    {
        SectorOffset[isec] = 0;
    }
    fafter = new TF1("f", "[0]*exp([1]/x)*pow(x,[2])", 0.2, 20);
    fafter->SetParameters(-8.25403, -5.4072, -0.33457);

    phi_conv_hist=nullptr; the_conv_hist=nullptr; r_conv_hist=nullptr; dzed_conv_hist=nullptr,phiv_conv_hist=nullptr;

    return;
}

void Run14AuAuLeptonCombyReco::InitWalk(PHCompositeNode *topNode)
{
    std::cout << __FILE__ << ":" << __LINE__ << " in InitRun" << std::endl;

    int runnumber = 0;

    const RunHeader *runheader =
        findNode::getClass<RunHeader>(topNode, "RunHeader");
    if (!runheader)
    {
        std::cout << PHWHERE << "Failed to find RunHeader Node" << std::endl;
    }
    runnumber = runheader->get_RunNumber();
    std::cout << "RecalEMCalTOF::InitRun: Run Number = " << runnumber << std::endl;
    fCDH = new emcCalibrationDataHelper(runnumber, false);

    char dummy[100];
    int run;
    int sector;
    int itowerid;
    float max = 0, peak = 0, sigma = 0;
    float walkconst = 0, walkconst2 = 0, woffset = 0, secondExp=0, normSecondExp=0;
    TOAD toad("Run14AuAuLeptonComby");
    
    std::string file_location0 = toad.location("WalkCorrection.txt");
    std::ifstream file_walk(file_location0.c_str());
    std::string file_location1 = toad.location("TowerByTower.txt");
    std::ifstream file_tower(file_location1.c_str());
    std::string file_location2 = toad.location("SectorBySector.txt");
    std::ifstream file_sector(file_location2.c_str());

    if (file_walk.is_open())
    {
        std::cout << "Recal Open '" << file_location0.c_str() << std::endl;
    }
    else
    {
        std::cout << "File " << file_location0.c_str() << " doesn't exist" << std::endl;
        exit(0);
    }
    if (file_tower.is_open())
    {
        std::cout << "Recal Open '" << file_location1.c_str() << std::endl;
    }
    else
    {
        std::cout << "File " << file_location1.c_str() << " doesn't exist" << std::endl;
        exit(0);
    }
    if (file_sector.is_open())
    {
        std::cout << "Recal Open '" << file_location2.c_str() << std::endl;
    }
    else
    {
        std::cout << "File " << file_location2.c_str() << " doesn't exist" << std::endl;
        exit(0);
    }

    while(file_walk>>dummy>>itowerid>>walkconst>>walkconst2>>woffset>>secondExp>>normSecondExp){

      Walk[itowerid] =walkconst;
      Walk2[itowerid]=walkconst2;
      Walk3[itowerid]=woffset;
      Walk4[itowerid]=secondExp;
      Walk5[itowerid]=normSecondExp;
    }
    file_walk.close();

    while (file_tower >> dummy >> itowerid >> max >> peak >> sigma)
    {
        T0Offset[itowerid] = peak;
        T0OffsetSigma[itowerid] = sigma;
    }
    file_tower.close();

    while (file_sector >> run >> dummy >> sector >> max >> peak >> sigma)
    {
        if (run == runnumber)
        {
            SectorOffset[sector] = peak;
        }
    }
    file_sector.close();

    se = Fun4AllServer::instance();
    std::cout << "RecalEMCalTOF::Init: " << "Book tdc_adc_PbGl histogram" << std::endl;
    phi_conv_hist = new TH2D("phi_conv_hist","phi_conv_hist",500,-0.5,0.5,50,0,5);
    se->registerHisto("phi_conv_hist", phi_conv_hist );
    the_conv_hist = new TH2D("the_conv_hist","the_conv_hist",500,-0.5,0.5,50,0,5);
    se->registerHisto("the_conv_hist", the_conv_hist );
    r_conv_hist = new TH2D("r_conv_hist","r_conv_hist",500,-1,49,50,0,5);
    se->registerHisto("r_conv_hist", r_conv_hist );
    dzed_conv_hist = new TH2D("dzed_conv_hist","dzed_conv_hist",500,-20,20,50,0,5);
    se->registerHisto("dzed_conv_hist", dzed_conv_hist );
    phiv_conv_hist = new TH2D("phiv_conv_hist","phiv_conv_hist",500,0,3.14159,50,0,5);
    se->registerHisto("phiv_conv_hist", phiv_conv_hist );

    return ;
}

void Run14AuAuLeptonCombyReco::Walking(PHCompositeNode *topNode)
{
    const PHGlobal *_phglobal_ptr =
        findNode::getClass<PHGlobal>(topNode, "PHGlobal");
    const emcClusterContainer* _emcClusterContainer_ptr =
        findNode::getClass<emcClusterContainer>(topNode, "emcClusterContainer");
    const emcTowerContainer* _emcTowerContainer_ptr =
        findNode::getClass<emcTowerContainer>(topNode, "emcHitContainer");
    if(!_phglobal_ptr||!_emcClusterContainer_ptr||!_emcTowerContainer_ptr)
    {
        std::cout<<"bad-vad nodes"<<std::endl;
    }
    float fVtx = _phglobal_ptr->getBbcZVertex();
    float bbct0 = _phglobal_ptr->getBbcTimeZero();
    
    int nclusters = _emcClusterContainer_ptr->size();
    int ntowers = _emcTowerContainer_ptr->size();
    
    for (int i = 0; i < nclusters; i++)
    { 
        emcClusterContent *cluster = _emcClusterContainer_ptr->getCluster(i);

        int clustercent = cluster->towerid(0);
        emcTowerContent *tower = nullptr;
        for (int itow = 0; itow < ntowers; itow++)
        {
            emcTowerContent *towertemp = _emcTowerContainer_ptr->getTower(itow);
            if (towertemp->towerid() == clustercent)
            {
                tower = towertemp;
            }
        }

        if (tower == nullptr)
            continue;
        
        int ifem, channel, isec;
        EmcIndexer::PXPXSM144CH(clustercent, ifem, channel);
        const emcCalFEM* LC = fCDH->getCalibration(ifem, "LCTofs");
        float lc = LC->getValueFast(channel, 0);
        lc = ((lc > 25. && lc < 65.) ? lc : 40.0) / 1000.;

        const int TDC = tower->TDC();
        const int ADC = tower->ADC();
        const float x = cluster->x();
        const float y = cluster->y();
        const float z = cluster->z() - fVtx;
        while(clustercent>24767) clustercent -=24768;
        if(clustercent<0) continue;
        if (clustercent < 15552)
        {
            isec = clustercent / (72 * 36); 
            //iy = (clustercent - isec * 2592) / 72;
            //iz = (clustercent - isec * 2592) % 72;
        }
        else
        {
            isec = 6 + (clustercent - 6 * 72 * 36) / (96 * 48);
            //iy = (clustercent - 15552 - 4608 * (isec - 6)) / 96;
            //iz = (clustercent - 15552 - 4608 * (isec - 6)) % 96;
        }
        int sector_clockwise = isec;
        if(isec<4) sector_clockwise = isec;
        else if(isec==4) sector_clockwise = 5; // sector E2
        else if(isec==5) sector_clockwise = 4;  // sector E3
        else if(isec==6) sector_clockwise = 7; // sector E0
        else if(isec==7) sector_clockwise = 6;  // sector E1
        else continue;
        
        const float d = sqrt(x * x + y * y + z * z);
        const float c = 29.979245829979; //[cm/ns]
        const float t_flash = d / c;
        //const float t0_offset = T0Offset[clustercent]; ///// i have no idea why??
        const float t0_offset = bbct0;
        const float sec_offset = SectorOffset[sector_clockwise];
        const float walk = Walk[clustercent]*(pow(ADC,Walk2[clustercent]))+Walk3[clustercent]*(pow(ADC,Walk4[clustercent]))+Walk5[clustercent]; //  new fitting function by balazs 5 parameters,  walk3 is tower offset in TDC
        //float fTime = -lc * (TDC - walk) - t0_offset - sec_offset - t_flash;
        float fTime =  lc * (4095 -TDC)- walk - t0_offset -sec_offset - t_flash; // Edited by Attia
        if (TDC < 0)
            fTime = -9999;
        // Afterburner
        //fTime = fTime - fafter->Eval(cluster->ecent());
        // std::cout << fTime << std::endl;
        
        cluster->set_tofcorr(fTime - 0*bbct0);
        
    }
    
    return;
}

void Run14AuAuLeptonCombyReco::StopWalking()
{
    delete fCDH;
}