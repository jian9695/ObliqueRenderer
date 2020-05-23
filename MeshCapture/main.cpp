
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
#include "GeometryVisitor.h"
using namespace osgDB;
using namespace OpenThreads;

int main(int argc, char** argv)
{

	QDir input_dir("F:/Oblique_Photogrammetry/weihai/Data/Tile_+000_-014/");
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	QFileInfoList list = input_dir.entryInfoList();

	osg::Geode* geode = new osg::Geode;
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		QString input_file = input_dir.absolutePath() + "/" + fileInfo.fileName();
		if (!input_file.endsWith(".osgb"))
			continue;
		std::string filename = fileInfo.baseName().toLocal8Bit().data();
		int len = std::string("Tile_+000_-014_").length();
		if (filename.length() < len + 3)
			continue;
		if (filename.substr(len, 3) != "L21")
			continue;
		printf("%s\n", input_file.toLocal8Bit().data());
		osg::Node* node = osgDB::readNodeFile(input_file.toLocal8Bit().data());

		GeometryVisitor visitor;
		node->accept(visitor);
		osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
		geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES, 0, visitor.Vertices->size()));
		geom->setVertexArray(visitor.Vertices.get());
		geom->setNormalArray(visitor.Normals.get());
		geom->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);
		geode->addDrawable(geom.get());
	}
	osgDB::writeNodeFile(*geode, "ljm.osg");

}
