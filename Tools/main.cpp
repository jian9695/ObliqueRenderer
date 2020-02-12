
#include <QApplication>
#include <QDesktopWidget>
#include <QSplashScreen>
#include <QTextCodec>


#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

// Open files in binary mode
#include <fcntl.h> /*  _O_BINARY */
#include <windows.h>
#include <dbghelp.h>
#include <time.h>

#include <QFileInfo>
#include "osgDB/WriteFile"
#include "qdir.h"
#include "qfileinfo.h"
#include "osgDB/WriteFile"
#include "osgEarth/URI"

#include "osg/ComputeBoundsVisitor"
#include "osgGA/StateSetManipulator"
#include "osgViewer/ViewerEventHandlers"
#include "osgUtil/Optimizer"
#include "osg/MatrixTransform"
#include "osgGA/TrackballManipulator"
#include "osg/Texture2D"
#include "osg/Image"
#include "osgDB/writeFile"
#include "osgDB/readFile"
#include "osg/ClampColor"
#include "osg/Depth"
#include "osgUtil/SmoothingVisitor"
#include "ogrsf_frmts.h"
#include <QFileinfo>
#include <map>
#include <fstream>
#include <osg/ShapeDrawable>
using namespace osgDB;
using namespace OpenThreads;

int main(int argc, char** argv)
{
	osg::Geode* geode = new osg::Geode();
	osg::TessellationHints* hints = new osg::TessellationHints;
	hints->setDetailRatio(0.5f);
	geode->addDrawable(new osg::ShapeDrawable(new osg::Box(osg::Vec3(0.0f, 0.0f, 100.0f), 1000, 1000, 200), hints));
	osg::ref_ptr<osg::Program> program = new osg::Program;
	char vertexShaderSource[] =
		"varying vec3 normal;\n"
		"void main(void)\n"
		"{\n"
		"normal = gl_Normal;\n"
		"gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\n"
		"}\n";
	char fragmentShaderSource[] =

		"varying vec3 normal;\n"
		"void main(void) \n"
		"{\n"
		"gl_FragColor.a = 1;\n"
		"gl_FragColor.rgb = (normal + 1) * 0.5;\n"
		"}\n";
	
	program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));
	program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexShaderSource));
	geode->getOrCreateStateSet()->setAttribute(program.get(), osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
	//geode->addDrawable(new osg::ShapeDrawable(new osg::Sphere(osg::Vec3(0.0f, 0.0f, 0.0f), 500), hints));
	//geode->addDrawable(new osg::ShapeDrawable(new osg::Cone(osg::Vec3(0.0f, 0.0f, 100.0f), 500, 200), hints));
	//geode->addDrawable(new osg::ShapeDrawable(new osg::Cylinder(osg::Vec3(0.0f, 0.0f, 100.0f), 500, 200), hints));
	osgDB::writeNodeFile(*geode, "box.osg");
	osgViewer::Viewer viewer;
	//viewer.getScene()->setDatabasePager()
	// add the state manipulator
	viewer.addEventHandler(new osgGA::StateSetManipulator(viewer.getCamera()->getOrCreateStateSet()));

	// add the thread model handler
	viewer.addEventHandler(new osgViewer::ThreadingHandler);

	// add the window size toggle handler
	viewer.addEventHandler(new osgViewer::WindowSizeHandler);

	// add the stats handler
	viewer.addEventHandler(new osgViewer::StatsHandler);
	viewer.setSceneData(geode);
	return viewer.run();
}
