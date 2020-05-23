
//#include "pch.h"
#include "ReadFileCallBackManager.h"
using namespace VGE;


ReaderWriter::ReadResult ReadFileCallBackManager::openArchive( const std::string& filename,ReaderWriter::ArchiveStatus status, unsigned int indexBlockSizeHint, const Options* useObjectCache )
{

	CallBackGroup::Map::iterator iter = g_mCallBackMap.begin();
	osgDB::ReaderWriter::ReadResult rr;
	while(iter != g_mCallBackMap.end())
	{
		if(dynamic_cast<OnReadFileCallBack*>(iter->second))
			(dynamic_cast<OnReadFileCallBack*>(iter->second))->beforeOpenArchive(rr,filename,status,indexBlockSizeHint,useObjectCache);
		iter++;
	}
	if(!rr.success())
	{
		rr = osgDB::ReadFileCallback::openArchive(filename,status,indexBlockSizeHint,useObjectCache);
	}
	iter = g_mCallBackMap.begin();
	while(iter != g_mCallBackMap.end())
	{
		if(dynamic_cast<OnReadFileCallBack*>(iter->second))
			(dynamic_cast<OnReadFileCallBack*>(iter->second))->afterOpenArchive(rr,filename,status,indexBlockSizeHint,useObjectCache);
		iter++;
	}
	return rr;
}

ReaderWriter::ReadResult ReadFileCallBackManager::readObject( const std::string& filename, const Options* options )
{
	CallBackGroup::Map::iterator iter = g_mCallBackMap.begin();
	osgDB::ReaderWriter::ReadResult rr;
	while(iter != g_mCallBackMap.end())
	{
		if(dynamic_cast<OnReadFileCallBack*>(iter->second))
			(dynamic_cast<OnReadFileCallBack*>(iter->second))->beforeReadObject(rr,filename,options);
		iter++;
	}
	if(!rr.success())
	{
		rr = osgDB::ReadFileCallback::readObject(filename,options);
	}
	iter = g_mCallBackMap.begin();
	while(iter != g_mCallBackMap.end())
	{
		if(dynamic_cast<OnReadFileCallBack*>(iter->second))
			(dynamic_cast<OnReadFileCallBack*>(iter->second))->afterReadObject(rr,filename,options);
		iter++;
	}
	return rr;
}

ReaderWriter::ReadResult ReadFileCallBackManager::readImage( const std::string& filename, const Options* options )
{
	CallBackGroup::Map::iterator iter = g_mCallBackMap.begin();
	osgDB::ReaderWriter::ReadResult rr;
	while(iter != g_mCallBackMap.end())
	{
		if(dynamic_cast<OnReadFileCallBack*>(iter->second))
			(dynamic_cast<OnReadFileCallBack*>(iter->second))->beforeReadImage(rr,filename,options);
		iter++;
	}
	if(!rr.success())
	{
		rr = osgDB::ReadFileCallback::readImage(filename,options);
	}
	iter = g_mCallBackMap.begin();
	while(iter != g_mCallBackMap.end())
	{
		if(dynamic_cast<OnReadFileCallBack*>(iter->second))
			(dynamic_cast<OnReadFileCallBack*>(iter->second))->afterReadImage(rr,filename,options);
		iter++;
	}
	return rr;
}

ReaderWriter::ReadResult ReadFileCallBackManager::readHeightField( const std::string& filename, const Options* options )
{
	CallBackGroup::Map::iterator iter = g_mCallBackMap.begin();
	osgDB::ReaderWriter::ReadResult rr;
	while(iter != g_mCallBackMap.end())
	{
		if(dynamic_cast<OnReadFileCallBack*>(iter->second))
			(dynamic_cast<OnReadFileCallBack*>(iter->second))->beforeReadHeightField(rr,filename,options);
		iter++;
	}
	if(!rr.success())
	{
		rr = osgDB::ReadFileCallback::readHeightField(filename,options);
	}
	iter = g_mCallBackMap.begin();
	while(iter != g_mCallBackMap.end())
	{
		if(dynamic_cast<OnReadFileCallBack*>(iter->second))
			(dynamic_cast<OnReadFileCallBack*>(iter->second))->afterReadHeightField(rr,filename,options);
		iter++;
	}
	return rr;
}

ReaderWriter::ReadResult ReadFileCallBackManager::readNode( const std::string& filename, const Options* options )
{
	CallBackGroup::Map::iterator iter = g_mCallBackMap.begin();
	osgDB::ReaderWriter::ReadResult rr;
	while(iter != g_mCallBackMap.end())
	{
		if(dynamic_cast<OnReadFileCallBack*>(iter->second))
			(dynamic_cast<OnReadFileCallBack*>(iter->second))->beforeReadNode(rr,filename,options);
		iter++;
	}
	if(!rr.success())
	{
		rr = osgDB::ReadFileCallback::readNode(filename,options);
	}
	iter = g_mCallBackMap.begin();
	while(iter != g_mCallBackMap.end())
	{
		if(dynamic_cast<OnReadFileCallBack*>(iter->second))
		   (dynamic_cast<OnReadFileCallBack*>(iter->second))->afterReadNode(rr,filename,options);
		iter++;
	}
	return rr;
}

ReaderWriter::ReadResult ReadFileCallBackManager::readShader( const std::string& filename, const Options* options )
{
	CallBackGroup::Map::iterator iter = g_mCallBackMap.begin();
	osgDB::ReaderWriter::ReadResult rr;
	while(iter != g_mCallBackMap.end())
	{
		if(dynamic_cast<OnReadFileCallBack*>(iter->second))
			(dynamic_cast<OnReadFileCallBack*>(iter->second))->beforeReadShader(rr,filename,options);
		iter++;
	}
	if(!rr.success())
	{
		rr = osgDB::ReadFileCallback::readShader(filename,options);
	}
	iter = g_mCallBackMap.begin();
	while(iter != g_mCallBackMap.end())
	{
		if(dynamic_cast<OnReadFileCallBack*>(iter->second))
			(dynamic_cast<OnReadFileCallBack*>(iter->second))->afterReadShader(rr,filename,options);
		iter++;
	}
	return rr;
}
ReadFileCallBackManager* ReadFileCallBackManager::g_pInstance;
ReadFileCallBackManager* ReadFileCallBackManager::instance()
{
	if(!g_pInstance)
		g_pInstance = new ReadFileCallBackManager;
	osgDB::Registry::instance()->setReadFileCallback(g_pInstance);
	return g_pInstance;
}
