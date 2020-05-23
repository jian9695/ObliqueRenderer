#pragma once
//#include  "vge_foundation_global.h"
#include <map>
namespace VGE
{
	//template <class CallBackType>
	//class VGE_FOUNDATION_EXPORT CallBackMangerTemplate
	//{
	//
	//public:
	//     void addCallBack(CallBackType* callback);
	//     void removeCallBack(CallBackType* callback);
	//	 typedef std::map<CallBackType*,CallBackType*> Map;
	//	 std::map<CallBackType*,CallBackType*> CallBackMap() { return g_mCallBackMap; }
	//	 void clear();
	//protected:
	//	std::map<CallBackType*,CallBackType*> g_mCallBackMap;
	//};

	class CallBackObject
	{
	public:
		CallBackObject(void){}
		virtual ~CallBackObject(void){}
	};
	class CallBackGroup
	{

	public:
		void addCallBack(CallBackObject* callback);
		void removeCallBack(CallBackObject* callback);
		typedef std::map<CallBackObject*,CallBackObject*> Map;
		std::map<CallBackObject*,CallBackObject*> CallBackMap() { return g_mCallBackMap; }
		void clear();
	protected:
		std::map<CallBackObject*,CallBackObject*> g_mCallBackMap;
	};
}

