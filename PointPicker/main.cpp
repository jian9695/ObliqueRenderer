
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
#include "PointPicker.h"
#include <osg/ShapeDrawable>
using namespace osgDB;
using namespace OpenThreads;
class BoundManager
{
public:

	static bool readBound(osg::BoundingBox& bound,std::string filename)
	{
		std::ifstream ifs(filename.data());
		if (!ifs.good())
			return false;
		ifs >> bound.xMin() >> bound.yMin() >> bound.zMin() >> bound.xMax() >> bound.yMax() >> bound.zMax();
		ifs.close();
		return true;
	}
	static void writeBound(osg::BoundingBox bound, std::string outfilename)
	{
		std::ofstream ofs(outfilename.data());
		ofs << bound.xMin() << " " << bound.yMin() << " " << bound.zMin() << " " << bound.xMax() << " " << bound.yMax() << " " << bound.zMax();
		ofs.close();
	}
};


class GroundGeo : public osg::Geode
{
public:
	osg::BoundingBox _mybb;
	GroundGeo(osg::BoundingBox bb)
		:osg::Geode()
	{
		_mybb = bb;
		if (_mybb.zMin() == _mybb.zMax())
		{
			_mybb.zMin() -= 0.0001;
			_mybb.zMax() += 0.0001;
		}
	}
	/** Return the Geode's bounding box, which is the union of all the
	* bounding boxes of the geode's drawables.*/
	inline const osg::BoundingBox& getBoundingBox() const
	{
		return _mybb;
	}
	virtual osg::BoundingSphere computeBound() const
	{
		return osg::BoundingSphere(_mybb);
	}

};
//osg::Geode* createRect(osg::BoundingBox bb,double baseHeight)
//{
//	osg::Geode* geode = new GroundGeo(bb);
//	osg::ref_ptr<osg::Program> program = new osg::Program;
//
//	char vertexShaderSource[] =
//		"void main(void)\n"
//		"{\n"
//		"   gl_Position    = gl_ModelViewProjectionMatrix *  gl_Vertex;\n"
//		"}\n";
//	char fragmentShaderSource[] =
//		"void main(void) \n"
//		"{\n"
//		"     gl_FragColor = vec4(0.5,0.5,0.5,1);\n"
//		"}\n";
//
//	program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));
//	program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexShaderSource));
//	osg::StateSet* ss = geode->getOrCreateStateSet();
//
//	ss->setAttribute(program.get(), osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
//
//	geode->setCullingActive(false);
//
//	osg::ref_ptr<osg::Vec3dArray> vertices = new osg::Vec3dArray;
//	osg::ref_ptr<osg::Vec3dArray> normals = new osg::Vec3dArray;
//	osg::ref_ptr<osg::Geometry> polyGeom = new osg::Geometry();
//	
//	vertices->push_back(osg::Vec3d(bb.xMin(), bb.yMin(), baseHeight));
//	vertices->push_back(osg::Vec3d(bb.xMin(), bb.yMax(), baseHeight));
//	vertices->push_back(osg::Vec3d(bb.xMax(), bb.yMax(), baseHeight));
//
//
//	vertices->push_back(osg::Vec3d(bb.xMin(), bb.yMin(), baseHeight));
//	vertices->push_back(osg::Vec3d(bb.xMax(), bb.yMax(), baseHeight));
//	vertices->push_back(osg::Vec3d(bb.xMax(), bb.yMin(), baseHeight));
//	normals->push_back(osg::Vec3d(0, 0, 1));
//
//	polyGeom->setVertexArray(vertices.get());
//	polyGeom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES, 0, vertices->size()));
//
//	polyGeom->setNormalArray(normals.get());
//	polyGeom->setNormalBinding(osg::Geometry::BIND_OVERALL);
//
//	//osg::CullFace* cull = new osg::CullFace(osg::CullFace::BACK;
//	geode->addDrawable(polyGeom.get());
//	//osg::ComputeBoundsVisitor cbv;
//	//geode->accept(cbv);
//	//osg::BoundingBox geodebb = cbv.getBoundingBox();
//	osg::BoundingBox geodebb = geode->getBoundingBox();
//	//polyGeom->getOrCreateStateSet()->setMode(GL_CULL_FACE, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
//	return geode;
//}

osg::Node* createRect(osg::BoundingBox bb, double baseHeight)
{
	osg::ref_ptr<osg::Geode> geode = new osg::Geode;
	osg::ref_ptr<osg::Program> program = new osg::Program;

	char vertexShaderSource[] =
		"void main(void)\n"
		"{\n"
		"   gl_Position    = gl_ModelViewProjectionMatrix *  gl_Vertex;\n"
		"}\n";
	char fragmentShaderSource[] =
		"void main(void) \n"
		"{\n"
		"     gl_FragColor = vec4(0.5,0.5,0.5,1);\n"
		"}\n";

	program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));
	program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexShaderSource));
	osg::StateSet* ss = geode->getOrCreateStateSet();

	ss->setAttribute(program.get(), osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);

	geode->addDrawable(new osg::ShapeDrawable(new osg::Box(osg::Vec3(0, 0, 0), 1)));

	osg::MatrixTransform* mat = new osg::MatrixTransform;
	mat->setMatrix(
		osg::Matrix::scale(bb.xMax() - bb.xMin(), bb.yMax() - bb.yMin(), 0.00001)
		* osg::Matrix::translate(bb.center()));
	mat->addChild(geode.get());
	return mat;
}

int main(int argc, char** argv)
{

	osg::ArgumentParser arguments(&argc, argv);
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
	PointPicker* pointPicker = new PointPicker;
	viewer.addEventHandler(pointPicker);
	osg::ref_ptr<osg::Group> root = new osg::Group;

	std::string infile = argv[1];
	std::string inshpfile = argv[2];
	
	if(QFileInfo(inshpfile.data()).exists())
		pointPicker->loadFromShapeFile(inshpfile.data());
	pointPicker->m_filename = inshpfile;
	printf("%s\n", inshpfile.data());
	osg::ref_ptr<osg::Node> node;
	if (argc > 3)
	{
		double baseHeight = atof(argv[3]);
		osg::BoundingBox bb;
		std::string bbfile = infile + ".bb";
		if (!BoundManager::readBound(bb, bbfile))
		{
			node = osgDB::readNodeFile(infile);
			osg::ComputeBoundsVisitor cbv;
			node->accept(cbv);
			bb = cbv.getBoundingBox();
			BoundManager::writeBound(bb, bbfile);
		}
		bb.zMin() = bb.zMax() = baseHeight;
		osg::ref_ptr<osg::Node> ground = createRect(bb, baseHeight);
		root->addChild(ground.get());
	}
	if(!node.valid())
		node = osgDB::readNodeFile(infile);

	if (!node || !node.valid())
	{
		return 1;
	}


	root->addChild(node.get());


	root->addChild(pointPicker->PointGeode().get());
	viewer.setSceneData(root.get());
	//viewer.setUpViewInWindow(50, 50, 1024, 768);

    viewer.run();
	exit(0);
	return 0;
}
