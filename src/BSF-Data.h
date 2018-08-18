/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: BSF-MR Skeleton
Module: BSF-Data.h (BSF Data)
Prefix: BD
Authors: Nadezhda A. Ezhova, Leonid B. Sokolinsky
Creation Date: 09.04.2017
==============================================================================*/

#pragma once

#include "BSF-Include.h"

// ========================= Constants ========================
#define BD_EXIT	true			// End of Work
// ========================= Variables ========================
static int	BD_listSize;		// Size of the list
static int	BD_rank;			// Rank of the current process MPI 
static int	BD_size;			// Number of MPI processes
static int	BD_masterRank;		// Rank of Master
static int	BD_numOfWorkers;	// Number of Workers = Number of MPI processes - 1
static int	BD_elemsPerWorker;	// Number of elements per worker
static int	BD_tailLength;		// Number of elements for last worker
static bool BD_exit;			// Indicator of ending the process
// ========================= Time Variables ========================
static double BD_t;				// Total time
static double BD_t_s_L;			// Sending time per one worker (including latency)
static double BD_t_S;			// Total sending time
static double BD_t_r_L;			// Receiving time per one worker (including latency)
static double BD_t_L;			// Latency
static double BD_t_W;			// Map time
static double BD_t_A_w;			// Total worker reduce time
static double BD_t_A_m;			// Total master reduce time
static double BD_t_p;
static int BD_iterCount;
// ========================= Data stuctures ========================
static PT_bsf_data_T BD_data;
// ========================= Lists ========================
PT_bsf_mapElem_T* BD_mapList;
BT_extendedReduceElem_T* BD_extendedReduceList;
BT_extendedReduceElem_T* BD_extendedReduceResult_P;
BT_order_T* BD_order;		// Orders for Workers
MPI_Status*  BD_status;
MPI_Request* BD_request;
