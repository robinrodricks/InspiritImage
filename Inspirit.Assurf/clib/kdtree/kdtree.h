/** @file   kdtree.h
 ** @brief  KD-tree
 ** @author Andrea Vedaldi
 **/

/* AUTORIGHTS
Copyright (C) 2007-10 Andrea Vedaldi and Brian Fulkerson

This file is part of VLFeat, available under the terms of the
GNU GPLv2, or (at your option) any later version.
*/

#ifndef VL_KDTREE_H
#define VL_KDTREE_H

#include "types.h"
#include "random.h"

#define VL_KDTREE_SPLIT_HEALP_SIZE 5
#define VL_KDTREE_MEDIAN 0
#define VL_KDTREE_MEAN 1

typedef struct _VlKDTreeNode VlKDTreeNode ;
typedef struct _VlKDTreeSplitDimension VlKDTreeSplitDimension ;
typedef struct _VlKDTreeDataIndexEntry VlKDTreeDataIndexEntry ;
typedef struct _VlKDForestSearchState VlKDForestSearchState ;

struct _VlKDTreeNode
{
  vl_uindex parent ;
  vl_index lowerChild ;
  vl_index upperChild ;
  unsigned int splitDimension ;
  double splitThreshold ;
  double lowerBound ;
  double upperBound ;
} ;

struct _VlKDTreeSplitDimension
{
  unsigned int dimension ;
  double mean ;
  double variance ;
} ;

struct _VlKDTreeDataIndexEntry
{
  vl_index index ;
  double value ;
} ;

/** @brief Neighbor of a query point */
typedef struct _VlKDForestNeighbor {
  double distance ;   /**< distance to the query point */
  vl_uindex index ;   /**< index of the neighbor in the KDTree data */
} VlKDForestNeighbor ;

typedef struct _VlKDTree
{
  VlKDTreeNode * nodes ;
  vl_size numUsedNodes ;
  vl_size numAllocatedNodes ;
  VlKDTreeDataIndexEntry * dataIndex ;
  unsigned int depth ;
} VlKDTree ;

struct _VlKDForestSearchState
{
  VlKDTree * tree ;
  vl_uindex nodeIndex ;
  double distanceLowerBound ;
} ;

/** @brief KDForest object */
typedef struct _VlKDForest
{
  vl_size dimension ;

  /* random number generator */
  VlRand * rand ;

  /* indexed data */
  vl_type dataType ;
  void const * data ;
  vl_size numData ;
  void (*distanceFunction)(void) ;

  /* tree structure */
  VlKDTree ** trees ;
  vl_size numTrees ;

  /* build */
  int thresholdingMethod ;
  VlKDTreeSplitDimension splitHeapArray [VL_KDTREE_SPLIT_HEALP_SIZE] ;
  vl_size splitHeapNumNodes ;
  vl_size splitHeapSize ;

  /* querying */
  VlKDForestSearchState * searchHeapArray ;
  vl_size searchHeapNumNodes ;
  vl_uindex searchId ;
  vl_uindex * searchIdBook ;

  vl_size searchMaxNumComparisons ;
  vl_size searchNumComparisons;
  vl_size searchNumRecursions ;
  vl_size searchNumSimplifications ;
} VlKDForest ;

/** @name Creatind and disposing
 ** @{ */
VlKDForest * vl_kdforest_new (vl_size dimension, vl_size numTrees) ;
void vl_kdforest_delete (VlKDForest * self) ;
/** @} */

/** @name Building and querying
 ** @{ */
void vl_kdforest_build (VlKDForest * self,
                                  vl_size numData,
                                  void const * data) ;
vl_size vl_kdforest_query (VlKDForest * self,
                                     VlKDForestNeighbor * neighbors,
                                     vl_size numNeighbors,
                                     void const * query) ;
/** @} */

/** @name Retrieving and setting parameters
 ** @{ */
inline vl_size vl_kdforest_get_depth_of_tree (VlKDForest const * self, vl_uindex treeIndex) ;
inline vl_size vl_kdforest_get_num_nodes_of_tree (VlKDForest const * self, vl_uindex treeIndex) ;
inline vl_size vl_kdforest_get_num_trees (VlKDForest const * self) ;
inline vl_size vl_kdforest_get_data_dimension (VlKDForest const * self) ;
inline vl_type vl_kdforest_get_data_type (VlKDForest const * self) ;
inline void vl_kdforest_set_max_num_comparisons (VlKDForest * self, vl_size n) ;
inline vl_size vl_kdforest_get_max_num_comparisons (VlKDForest * self) ;
inline void vl_kdforest_set_thresholding_method (VlKDForest * self, int method) ;
inline int vl_kdforest_get_thresholding_method (VlKDForest const * self) ;
/** @} */

inline double vec_compare_function(vl_size dimension, double const * X, double const * Y);


/* VL_KDTREE_H */
#endif
