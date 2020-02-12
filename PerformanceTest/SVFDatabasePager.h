#pragma once
#include <Windows.h>
#include <osgViewer/Renderer>
#include <osgGA/TrackballManipulator>
#include <osgDB/WriteFile>
#include <osgViewer/ViewerEventHandlers>
using namespace osgDB;
using namespace OpenThreads;

class VGEDatabasePager : public osgDB::DatabasePager
{
public:
	VGEDatabasePager()
		:osgDB::DatabasePager(){}
	void frame();
	void pause();
	void resume();
};
