#pragma once
#include <Windows.h>
#include <osg/Texture2D>
#include <osg/Geometry>
#include <osg/Camera>
#include <osg/Geode>
#include <osgViewer/Viewer>
class ScreenOverlay : public osg::Geode
{
public:
	ScreenOverlay(osgViewer::Viewer* viewer);
	~ScreenOverlay();
	void setTextureLayer(osg::Texture* tex);
};

