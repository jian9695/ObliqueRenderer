#pragma once
#include <Windows.h>
#include <osgGA/GUIEventHandler>
#include <osgViewer/Viewer>
#include "osg/Texture2D"
#include "osg/Image"

class MeasureTools3D : public osgGA::GUIEventHandler
{
public:
	MeasureTools3D(void);
	~MeasureTools3D(void);
	osg::ref_ptr<osg::Geode> PointGeode()  { return g_pPointGeode; }
	std::string m_filename;
	void loadFromShapeFile(std::string filename);
protected:
	bool handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter& aa);
	void updateGeometry();
protected:
	std::vector<osg::Vec3d> g_mPoints;
	std::vector<osg::Vec3d> g_mNormals;
	osg::ref_ptr<osg::Geode> g_pPointGeode;

};

