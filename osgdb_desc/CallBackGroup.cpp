
#include "CallBackGroup.h"
using namespace VGE;
//
//template <class CallBackType>
//void CallBackGroup<CallBackType>::removeCallBack( CallBackType* callback )
//{
//	if(!callback)
//		return;
//	if(g_mCallBackMap.find(callback) != g_mCallBackMap.end())
//	{
//	   g_mCallBackMap.erase(callback);
//	   delete callback;
//	}
//}
//
//template <class CallBackType>
//void CallBackGroup<CallBackType>::addCallBack( CallBackType* callback )
//{
//	if(!callback)
//		return;
//	g_mCallBackMap[callback]=callback;
//}
//
//template <class CallBackType>
//void VGE::CallBackGroup<CallBackType>::clear()
//{
//	CallBackGroup::Map::iterator iter = g_mCallBackMap.begin();
//	while(iter != g_mCallBackMap.end())
//	{
//		if(iter->second)
//			delete iter->second;
//		iter++;
//	}
//	g_mCallBackMap.clear();
//}



void CallBackGroup::removeCallBack( CallBackObject* callback )
{
	if(!callback)
		return;
	if(g_mCallBackMap.find(callback) != g_mCallBackMap.end())
		g_mCallBackMap.erase(callback);
}

void CallBackGroup::addCallBack( CallBackObject* callback )
{
	if(!callback)
		return;
	g_mCallBackMap[callback]=callback;
}

void CallBackGroup::clear()
{
	CallBackGroup::Map::iterator iter = g_mCallBackMap.begin();
	while(iter != g_mCallBackMap.end())
	{
		if(iter->second)
			delete iter->second;
		iter++;
	}
	g_mCallBackMap.clear();
}
