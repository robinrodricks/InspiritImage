//#include "kdtree/kdtree.h"

VlKDForest * refTree;

int getMatchesByTree(RefObject *obj, IPoint *currP, IPoint *refP, const int num1, const int num2, IPointMatch *matches, int prevMatchedCount, const double point_match_factor)
{
	int i, matchedCount = prevMatchedCount;
	
	double pmap[obj->pointsCount];
	int imap[obj->pointsCount];
	
	for(i=0; i<obj->pointsCount; i++)
	{
		pmap[i] = 100.0;
		imap[i] = -1;
	}
	for(i=0; i<prevMatchedCount; i++)
	{
		pmap[ matches[i].second->localIndex ] = matches[i].confidence;
		imap[ matches[i].second->localIndex ] = i;
	}
	
	VlKDForestNeighbor neighbors[2];
	
	refTree = obj->kdf;
	IPointMatch *match;
	
	for(i = 0; i < num1; i++)
	{
		IPoint *curr = &currP[i];
		
		if(curr->matched) continue;
		
		vl_kdforest_query (refTree, &*neighbors, 2, curr->descriptor);
		
		const double ratio = neighbors[0].distance / neighbors[1].distance;
		const int match_idx = neighbors[0].index;
		const double curr_dist = neighbors[0].distance;
		
		if(ratio < point_match_factor && pmap[match_idx] > curr_dist)
		{
			
			if(imap[match_idx] > -1)
			{
				match = &matches[ imap[match_idx] ];
				
				match->first = curr;
				match->second = &obj->points[match_idx];
				match->confidence = curr_dist;
				match->lk_tracked = 0;
				match->tracked = 0;
				match->prev = 0;
				match->wasGood = 0;
			}
			else
			{
				match = &matches[matchedCount];
				
				match->first = curr;
				match->second = &obj->points[match_idx];
				match->confidence = curr_dist;
				match->lk_tracked = 0;
				match->tracked = 0;
				match->prev = 0;
				match->wasGood = 0;
				
				imap[match_idx] = matchedCount;
				
				matchedCount++;
			}
			
			pmap[match_idx] = curr_dist;
		}
	}
	
	double max_dist = 0.0;
	for(i = 0; i < matchedCount; i++) 
	{
		match = &matches[i];
		if(max_dist < match->confidence) 
		{
			max_dist = match->confidence;
		}
	}
	for(i = 0; i < matchedCount; i++)
	{
		match = &matches[i];
		match->normConfidence = match->confidence / max_dist;
	}
	
	return matchedCount;
}

int getMatchesByMultiTree(IPoint *currP, IPoint *refP, const int num1, const int num2, IPointMatch *matches, int prevMatchedCount, const double point_match_factor)
{
	int i, t, matchedCount = prevMatchedCount;
	double screen_dist[num1];
	int prev_match_ind[num1];
	int prev_obj_ind[num1];
	int real_match_idx;
	
	VlKDForestNeighbor neighbors[2];
	
	for(i = 0; i < num1; i++) screen_dist[i] = 100.0, prev_match_ind[i]=0, prev_obj_ind[i]=-1;
	
	for(t = 0; t < referenceCount; t++)
	{
		RefObject *obj = &refObjectsMap[t];
		refTree = obj->kdf;
		
		double pmap[obj->pointsCount];
		int imap[obj->pointsCount];
		
		obj->prevMatchedPointsCount = obj->matchedPointsCount;
		obj->matchedPointsCount = 0;
		
		for(i=0; i<obj->pointsCount; i++)
		{
			pmap[i] = 100.0;
			imap[i] = -1;
		}
		for(i=0; i<prevMatchedCount; i++)
		{
			if(obj->index == matches[i].second->refIndex)
			{
				IPointMatch *pfm = &matches[i];
				pmap[ pfm->second->localIndex ] = pfm->confidence;
				imap[ pfm->second->localIndex ] = obj->matchedPointsCount;
				memcpy(obj->matches + obj->matchedPointsCount, pfm, sizeof(IPointMatch));
				obj->matchedPointsCount++;
			}
		}
	
		for(i = 0; i < num1; i++)
		{
			if(screen_dist[i] < 3.0) continue;
			
			IPoint *curr = &currP[i];
			
			if(curr->matched) continue;
			
			vl_kdforest_query (refTree, &*neighbors, 2, curr->descriptor);
			
			const double ratio = neighbors[0].distance / neighbors[1].distance;
			const int match_idx = neighbors[0].index;
			const double curr_dist = neighbors[0].distance;
			
			if(ratio < point_match_factor && pmap[match_idx] > curr_dist)
			{
				if(t > 0 && screen_dist[i] > curr_dist && prev_obj_ind[i] > -1)
				{
					RefObject *obj_prv = &refObjectsMap[ prev_obj_ind[i] ];
					if(obj_prv->matchedPointsCount - 1 > prev_match_ind[i])
					{
						memcpy( obj_prv->matches+prev_match_ind[i], obj_prv->matches+prev_match_ind[i]+1, (obj_prv->matchedPointsCount-prev_match_ind[i]) * sizeof(IPointMatch) );
					}
					obj_prv->matchedPointsCount--;
				}
				
				if(imap[match_idx] > -1)
				{
					real_match_idx = imap[match_idx];
				}
				else
				{
					real_match_idx = obj->matchedPointsCount;					
					imap[match_idx] = real_match_idx;
					
					obj->matchedPointsCount++;
					matchedCount++;
				}
				
				IPointMatch *match = &obj->matches[ real_match_idx ];
				match->first = curr;
				match->second = &obj->points[match_idx];
				match->confidence = curr_dist;
				match->lk_tracked = 0;
				match->tracked = 0;
				match->prev = 0;
				match->wasGood = 0;
				
				if(screen_dist[i] > curr_dist) 
				{
					prev_match_ind[i] = real_match_idx;
					screen_dist[i] = curr_dist;
					prev_obj_ind[i] = t;
				}
				
				pmap[match_idx] = curr_dist;
			}
		}
		
		if(obj->matchedPointsCount > 4)
		{
			double max_dist = 0.0;
			for(i = 0; i < obj->matchedPointsCount; i++) 
			{
				if(max_dist < obj->matches[i].confidence) 
				{
					max_dist = obj->matches[i].confidence;
				}
			}
			for(i = 0; i < obj->matchedPointsCount; i++)
			{
				obj->matches[i].normConfidence = obj->matches[i].confidence / max_dist;
			}
		}
	}
	
	return matchedCount;
}
