/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: BSF-MR Skeleton
Module: BSF-Types.h (BSF Types)
Prefix: BT
Authors: Nadezhda A. Ezhova, Leonid B. Sokolinsky
Creation Date: 09.04.2017
==============================================================================*/
#pragma once
#include "BSF-Include.h"

struct BT_order_T {  
	char exit;		// true, if worker must stop
	int index;		// index of the first MAP sublist element
	int length;		// length of MAP sublist
	PT_bsf_data_T data;
};

struct BT_extendedReduceElem_T { // Extended element type of reduce list
	PT_bsf_reduceElem_T elem;	// Element of reduce list
	int counter;			// Counter
	double t_W;
	double t_A_w;
};
