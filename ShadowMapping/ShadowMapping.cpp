
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
#include "osg/LineSegment"
#include <osg/ShapeDrawable>
#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include "ogrsf_frmts.h"
#include "ModelLoader.h"
#include "ScreenOverlay.h"
#include "osgUtil/SmoothingVisitor"
#include "GrassSolar.h"
#include "osg/Material"
#include "GDAL_DS.h"
#include "osg/Texture1D"

using namespace osgDB;
using namespace OpenThreads;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  SortFileRequestFunctor
//

void printfVec3(osg::Vec3 p1, osg::Vec3 p2)
{
		printf("(%f,%f,%f),(%f,%f,%f)\n", p1.x(), p1.y(), p1.z(), p2.x(), p2.y(), p2.z());
}

class FaceVisitor
{
public:
		FaceVisitor(osg::PrimitiveSet* pset)
		{
				pSet = pset;
		}
		std::vector<unsigned int> getFaceIndices()
		{

				std::vector<unsigned int> indices;
				unsigned int idx = pSet->getNumIndices();
				unsigned int numofprims;


				if (pSet->getMode() == osg::PrimitiveSet::TRIANGLE_FAN)
				{
						numofprims = idx - 2;
						for (unsigned int i = 0; i < numofprims; i++)
						{
								indices.push_back(pSet->index(0));
								indices.push_back(pSet->index(i + 1));
								indices.push_back(pSet->index(i + 2));
						}
				}
				else if (pSet->getMode() == osg::PrimitiveSet::TRIANGLES)
				{
						for (unsigned int i = 0; i < idx; i++)
						{
								indices.push_back(pSet->index(i));
						}
				}
				else if (pSet->getMode() == osg::PrimitiveSet::TRIANGLE_STRIP)
				{
						numofprims = (idx - 2) / 2;
						for (unsigned int i = 0; i < numofprims; i++)
						{
								indices.push_back(pSet->index(i * 2));
								indices.push_back(pSet->index(i * 2 + 1));
								indices.push_back(pSet->index(i * 2 + 2));

								indices.push_back(pSet->index(i * 2 + 2));
								indices.push_back(pSet->index(i * 2 + 1));
								indices.push_back(pSet->index(i * 2 + 3));
						}

						if (numofprims + 2 % 2 != 0)
						{
								indices.push_back(pSet->index(idx - 3));
								indices.push_back(pSet->index(idx - 2));
								indices.push_back(pSet->index(idx - 1));
						}
				}
				return indices;
		}

private:
		osg::PrimitiveSet* pSet;
};

class GeometryVisitor : public osg::NodeVisitor
{
private:
		osg::ref_ptr<osg::Group> m_group;
public:
		GeometryVisitor() :osg::NodeVisitor(TRAVERSE_ALL_CHILDREN) { setNodeMaskOverride(0xffffffff); }

		void initialize()
		{
				m_group = new osg::Group();
		}

		void save(std::string filename)
		{
				osgDB::writeNodeFile(*m_group, filename);
		}

		virtual void apply(osg::Geode& node)
		{
				osg::MatrixList matlist = node.getWorldMatrices();
				osg::Matrix matWorld = osg::Matrix::identity();

				for (unsigned int i = 0; i < matlist.size(); i++)
				{
						matWorld = matlist[i] * matWorld;
				}
				bool isIdentify = (matWorld.isIdentity());
				osg::Texture* texture = dynamic_cast<osg::Texture*>(node.getStateSet()->getTextureAttribute(0, osg::StateAttribute::TEXTURE));
				osg::Texture* texture2d = dynamic_cast<osg::Texture2D*>(texture);
				m_group->addChild(&node);
				for (int i = 0; i < node.getNumDrawables(); i++)
				{
						osg::Drawable* drawable = node.getDrawable(i);
						osg::Geometry* geom = dynamic_cast<osg::Geometry*>(drawable);
						if (!geom) continue;
						//osg::ref_ptr<osg::Geometry> geomCpy = new osg::Geometry;

						osg::Vec3Array* vertices = (osg::Vec3Array*) geom->getVertexArray();
						osg::Vec2Array* uvs = (osg::Vec2Array*) geom->getTexCoordArray(0);
						osg::Vec3Array* normals = (osg::Vec3Array*) geom->getNormalArray();
						if (!normals)
						{
								osgUtil::SmoothingVisitor::smooth(*geom);
						}
						if (!texture2d)
						{
								texture = dynamic_cast<osg::Texture*>(geom->getStateSet()->getTextureAttribute(0, osg::StateAttribute::TEXTURE));
								texture2d = dynamic_cast<osg::Texture2D*>(texture);
						}
						for (int j = 0; j < geom->getNumPrimitiveSets(); ++j)
						{
								osg::PrimitiveSet* primitiveSet = geom->getPrimitiveSet(j);
								FaceVisitor visitor(primitiveSet);
								std::vector<unsigned int> indices = visitor.getFaceIndices();
								
								for (unsigned int k = 0; k < indices.size(); k++)
								{
										unsigned int index = indices[k];
										osg::Vec3 pos = (*vertices)[index];
										if (!isIdentify)
												pos = pos * matWorld;
										osg::Vec3 normal = (*normals)[index];
										osg::Vec2 uv = (*uvs)[index];
								}
						}
				}
		}
};

class ObliqueCamera : public osg::Camera
{
public:
		double _altAngle;
		double _azimuthAngle;
		osg::Vec3d _eye;
		osg::Vec3d _center;
		ObliqueCamera()
		{
				this->setReferenceFrame(osg::Transform::ABSOLUTE_RF);
		}

		static osg::Vec3 getLightPos(double altAngle = 0, double azimuthAngle = 0)
		{
				osg::Vec3 lightDir;
				lightDir.z() = cos(osg::DegreesToRadians(90.0 - altAngle));
				double projectedLenghOnXY = cos(osg::DegreesToRadians(altAngle));
				lightDir.y() = projectedLenghOnXY * cos(osg::DegreesToRadians(azimuthAngle));
				//lightDir.x() = sqrt((projectedLenghOnXY * projectedLenghOnXY) - lightDir.y() * lightDir.y());
				lightDir.x() = projectedLenghOnXY * cos(osg::DegreesToRadians(90 - azimuthAngle));
				lightDir.normalize();
				return lightDir;
		}

		void reset(osg::BoundingBoxd bb, double altAngle = 0, double azimuthAngle = 0)
		{
				//altAngle = 0;
				//bb = osg::BoundingBoxd(-100, -100, -100, 100, 100, 100);
				osg::Vec3d minBB(bb.xMin(), bb.yMin(), bb.zMin());
				osg::Vec3d maxBB(bb.xMax(), bb.yMax(), bb.zMax());
				altAngle = osg::DegreesToRadians(altAngle);
			 azimuthAngle = osg::DegreesToRadians(azimuthAngle);
				osg::Vec3d forward(0, -cos(altAngle), sin(altAngle));
				osg::Vec3d right(1, 0, 0);
				osg::Vec3d up = forward ^ right;
				minBB -= bb.center();
				maxBB -= bb.center();
				BBWrapper localBB(minBB, maxBB);
				BBWrapper cameraSpaceBB;
				cameraSpaceBB.init();
				osg::Vec3d eye = -forward * bb.radius();
				osg::Matrixd viewMatrix = osg::Matrixd::lookAt(eye, osg::Vec3d(0, 0, 0), up);
				for (size_t i = 0; i < 8; i++)
				{
						cameraSpaceBB.expandBy(localBB.corner(i) * viewMatrix);
				}

				osg::Matrixd inverseViewMatrix = osg::Matrixd::inverse(viewMatrix);
				eye = (cameraSpaceBB.center() - osg::Vec3d(0, 0, cameraSpaceBB.zhalfsize())) * inverseViewMatrix + bb.center();
				osg::Vec3d viewDirection = ((cameraSpaceBB.center() + osg::Vec3d(0, 0, cameraSpaceBB.zhalfsize())) * inverseViewMatrix + bb.center()) - eye;
				osg::Vec3d center = cameraSpaceBB.center() * inverseViewMatrix + bb.center();
				this->setViewMatrixAsLookAt(eye, center, up);
				double xsize = (((cameraSpaceBB.center() - osg::Vec3d(cameraSpaceBB.xhalfsize(), 0, 0)) * inverseViewMatrix) - ((cameraSpaceBB.center() + osg::Vec3d(cameraSpaceBB.xhalfsize(), 0, 0)) * inverseViewMatrix)).length();
				double ysize = (((cameraSpaceBB.center() - osg::Vec3d(0, cameraSpaceBB.yhalfsize(), 0)) * inverseViewMatrix) - ((cameraSpaceBB.center() + osg::Vec3d(0, cameraSpaceBB.yhalfsize(), 0)) * inverseViewMatrix)).length();
				double zsize = (((cameraSpaceBB.center() - osg::Vec3d(0, 0, cameraSpaceBB.zhalfsize())) * inverseViewMatrix) - ((cameraSpaceBB.center() + osg::Vec3d(0, 0, cameraSpaceBB.zhalfsize())) * inverseViewMatrix)).length();
				this->setProjectionMatrixAsOrtho(-xsize * 0.5, xsize*0.5, -ysize * 0.5, ysize*0.5, 0, zsize*3);
				_eye = eye;
				_center = center;
		}

		void reset2(osg::BoundingBoxd bb, double altAngle = 0, double azimuthAngle = 0)
		{
				osg::Vec3d lightDir;
				lightDir.z() = cos(osg::DegreesToRadians(90.0 - altAngle));
				double projectedLenghOnXY = cos(osg::DegreesToRadians(altAngle));
				lightDir.y() = projectedLenghOnXY * cos(osg::DegreesToRadians(azimuthAngle));
				//lightDir.x() = sqrt((projectedLenghOnXY * projectedLenghOnXY) - lightDir.y() * lightDir.y());
				lightDir.x() = projectedLenghOnXY * cos(osg::DegreesToRadians(90 - azimuthAngle));
				lightDir.normalize();

				osg::Vec3d lightXY(lightDir.x(), lightDir.y(), 0);
				lightXY.normalize();
				osg::Vec3d right = lightXY * osg::Matrix::rotate(osg::DegreesToRadians(90.0), osg::Vec3d(0, 0, 1));
				//printfVec3(lightXY, right);

				//altAngle = 0;
				//bb = osg::BoundingBoxd(-100, -100, -100, 100, 100, 100);
				osg::Vec3d minBB(bb.xMin(), bb.yMin(), bb.zMin());
				osg::Vec3d maxBB(bb.xMax(), bb.yMax(), bb.zMax());
				altAngle = osg::DegreesToRadians(altAngle);
				azimuthAngle = osg::DegreesToRadians(azimuthAngle);
				
				osg::Vec3d forward = lightDir;
				osg::Vec3d up = forward ^ right;
				minBB -= bb.center();
				maxBB -= bb.center();
				BBWrapper localBB(minBB, maxBB);
				BBWrapper cameraSpaceBB;
				cameraSpaceBB.init();
				osg::Vec3d eye = -forward * bb.radius();
				osg::Matrixd viewMatrix = osg::Matrixd::lookAt(eye, osg::Vec3d(0, 0, 0), up);
				for (size_t i = 0; i < 8; i++)
				{
						cameraSpaceBB.expandBy(localBB.corner(i) * viewMatrix);
				}

				osg::Matrixd inverseViewMatrix = osg::Matrixd::inverse(viewMatrix);
				eye = (cameraSpaceBB.center() - osg::Vec3d(0, 0, cameraSpaceBB.zhalfsize())) * inverseViewMatrix + bb.center();
				osg::Vec3d viewDirection = ((cameraSpaceBB.center() + osg::Vec3d(0, 0, cameraSpaceBB.zhalfsize())) * inverseViewMatrix + bb.center()) - eye;
				osg::Vec3d center = cameraSpaceBB.center() * inverseViewMatrix + bb.center();
				this->setViewMatrixAsLookAt(eye, center, up);
				double xsize = (((cameraSpaceBB.center() - osg::Vec3d(cameraSpaceBB.xhalfsize(), 0, 0)) * inverseViewMatrix) - ((cameraSpaceBB.center() + osg::Vec3d(cameraSpaceBB.xhalfsize(), 0, 0)) * inverseViewMatrix)).length();
				double ysize = (((cameraSpaceBB.center() - osg::Vec3d(0, cameraSpaceBB.yhalfsize(), 0)) * inverseViewMatrix) - ((cameraSpaceBB.center() + osg::Vec3d(0, cameraSpaceBB.yhalfsize(), 0)) * inverseViewMatrix)).length();
				double zsize = (((cameraSpaceBB.center() - osg::Vec3d(0, 0, cameraSpaceBB.zhalfsize())) * inverseViewMatrix) - ((cameraSpaceBB.center() + osg::Vec3d(0, 0, cameraSpaceBB.zhalfsize())) * inverseViewMatrix)).length();
				this->setProjectionMatrixAsOrtho(-xsize * 0.5, xsize*0.5, -ysize * 0.5, ysize*0.5, 0, zsize * 3);
				_eye = eye;
				_center = center;
		}
};


osg::Vec3d getLightDir(double altAngle = 0, double azimuthAngle = 0)
{
		osg::Vec3d lightDir;
		lightDir.z() = cos(osg::DegreesToRadians(90.0 - altAngle));
		double projectedLenghOnXY = cos(osg::DegreesToRadians(altAngle));
		lightDir.y() = projectedLenghOnXY * cos(osg::DegreesToRadians(azimuthAngle));
		//lightDir.x() = sqrt((projectedLenghOnXY * projectedLenghOnXY) - lightDir.y() * lightDir.y());
		lightDir.x() = projectedLenghOnXY * cos(osg::DegreesToRadians(90 - azimuthAngle));
		lightDir.normalize();
		return lightDir;
}

osg::Texture1D* createRandomColors(int num)
{
		srand(time(0));
		osg::ref_ptr<osg::Image> colorImage = new osg::Image;
		colorImage->allocateImage(num, 1, 1, GL_RGB, GL_UNSIGNED_BYTE);
		unsigned char* colorData = colorImage->data();
		for (long i = 0; i < num; i++)
		{
				*colorData++ = (char)(rand() % 256);
				*colorData++ = (char)(rand() % 256);
				*colorData++ = (char)(rand() % 256);
		}
		osg::Texture1D* texColor = new osg::Texture1D(colorImage.get());
		texColor->setResizeNonPowerOfTwoHint(false);
		return texColor;
}

class HeightVisitor : public osg::NodeVisitor
{
public:
		BBWrapper BB;
		std::vector<double> HeightField;
		int Rows;
		int Columns;
		double CellSize;
		double NoData;
		HeightVisitor() :osg::NodeVisitor(TRAVERSE_ALL_CHILDREN)
		{
				setNodeMaskOverride(0xffffffff);
				GDALAllRegister();
		}

		void Create(osg::BoundingBoxd bb, double cellsize)
		{
				NoData = -9999;
				CellSize = cellsize;
				BB = BBWrapper(bb);
				Rows = (int)(BB.ysize() / cellsize);
				while (Rows * cellsize < BB.ysize())
				{
						Rows++;
				}

				Columns = (int)(BB.xsize() / cellsize);
				while (Columns * cellsize < BB.xsize())
				{
						Columns++;
				}

				HeightField.resize(Rows * Columns);
				for (size_t i = 0; i < Rows * Columns; i++)
				{
						HeightField[i] = NoData;
				}
		}

		bool GetIndex(double x, double y, double& col, double& row)
		{
				if (x < BB.xMin() || x > BB.xMax() || y < BB.yMin() || y > BB.yMax())
						return false;
				col = (int)((x - BB.xMin()) / CellSize);
				row = (int)((BB.yMax() - y) / CellSize);
				if (col > Columns - 1)
						col = Columns - 1;
				if (row > Rows - 1)
						row = Rows - 1;
				return true;
		}

		void GetHeightRange(const osg::BoundingBoxd& bb, double& min, double& max) 
		{
				GetHeightRange(bb.xMin(), bb.xMax(), bb.yMin(), bb.yMax(), min, max);
		}

		void GetHeightRange(const BBWrapper& bb, double& min, double& max)
		{
				GetHeightRange(bb.xMin(), bb.xMax(), bb.yMin(), bb.yMax(), min, max);
		}

		void GetHeightRange(double xmin, double xmax, double ymin, double ymax, double& min, double& max)
		{
				min = DBL_MAX;
				max = -DBL_MAX;
				double x = xmin;
				while (x <= xmax)
				{
						double y = ymin;
						while (y <= ymax)
						{
								double col, row;
								if (!GetIndex(x, y, col, row))
								{
										y += CellSize;
										continue;
								}

								int index = col + row * Columns;
								double height = HeightField[index];
								if (height <= NoData)
								{
										y += CellSize;
										continue;
								}
								if (min > height)
										min = height;
								if (max < height)
										max = height;
								y += CellSize;
						}
						x += CellSize;
				}
		}

		virtual void apply(osg::Geode& node)
		{
				osg::MatrixList matlist = node.getWorldMatrices();
				osg::Matrix matWorld = osg::Matrix::identity();

				/*计算goede的坐标变换 ？*/
				for (unsigned int i = 0; i < matlist.size(); i++)
				{
						matWorld = matlist[i] * matWorld;
				}

				for (int i = 0; i < node.getNumDrawables(); i++)
				{
						osg::Drawable* drawable = node.getDrawable(i);
						osg::Geometry* geom = dynamic_cast<osg::Geometry*>(drawable);
						if (!geom) continue;

						if (dynamic_cast<osg::Vec3Array*>(geom->getVertexArray()))
						{
								osg::Vec3Array* vertices = (osg::Vec3Array *)geom->getVertexArray();
								for (unsigned int k = 0; k < vertices->size(); k++)
								{
										osg::Vec3 pos = (*vertices)[k];
										pos = pos * matWorld;
										if (pos.z() <= NoData)
												continue;
										double col, row;
										if (!GetIndex(pos.x(), pos.y(), col, row))
												continue;
										int index = col + row * Columns;
										if (HeightField[index] < pos.z())
												HeightField[index] = pos.z();
								}
						}
						else if (dynamic_cast<osg::Vec3dArray*>(geom->getVertexArray()))
						{
								osg::Vec3dArray* vertices = (osg::Vec3dArray*)geom->getVertexArray();
								for (unsigned int k = 0; k < vertices->size(); k++)
								{
										osg::Vec3d pos = (*vertices)[k];
										pos = pos * matWorld;
										if (pos.z() <= NoData)
												continue;
										double col, row;
										if (!GetIndex(pos.x(), pos.y(), col, row))
												continue;
										int index = col + row * Columns;
										if (HeightField[index] < pos.z())
												HeightField[index] = pos.z();
								}
						}

				}
		}

		void Write(std::string outfile)
		{
				GDAL_DS<double>* ds = new GDAL_DS<double>();
				GDAL_DSInfo info;
			 memset(info.adfGeoTransform, 0, sizeof(double) * 6);
				info.adfGeoTransform[0] =	BB.xMin();
				info.adfGeoTransform[3] = BB.yMax();
				info.adfGeoTransform[1] = CellSize;
				info.adfGeoTransform[5] = -CellSize;
				info.filename = outfile;
				info.ncols = Columns;
				info.nrows = Rows;
				info.numbands = 1;
				ds->setDSInfo(&info);
				ds->create(outfile);
				ds->writeData(1, &HeightField[0], NoData);
				delete ds;
		}

private:
};


bool _canUpdateLightDir = false;

class MyKeyboardHandler : public osgGA::GUIEventHandler
{
public:

		MyKeyboardHandler() : osgGA::GUIEventHandler()
		{

		}

		bool handle(const osgGA::GUIEventAdapter& ea, osgGA::GUIActionAdapter& aa)
		{
				osgViewer::Viewer* viewer = dynamic_cast<osgViewer::Viewer*>(&aa);
				if (!viewer)
						return false;

				switch (ea.getEventType())
				{
						case(osgGA::GUIEventAdapter::KEYDOWN):
						{
								if (ea.getKey() == 'm')
								{
										_canUpdateLightDir = true;
								}
						}
				}

				return false;
		}
};

#include "osg/Point"
#include "osg/LineWidth"

osg::Geometry* createLine(std::vector<osg::Vec3> points)
{
		// create Geometry object to store all the vertices and lines primitive.
		osg::Geometry* linesGeom = new osg::Geometry();

		// this time we'll preallocate the vertex array to the size
		// and then use an iterator to fill in the values, a bit perverse
		// but does demonstrate that we have just a standard std::vector underneath.
		osg::Vec3Array* vertices = new osg::Vec3Array();
		for (size_t i = 0; i < points.size(); i++)
		{
				vertices->push_back(points[i]);
		}

		// pass the created vertex array to the points geometry object.
		linesGeom->setVertexArray(vertices);

		// set the colors as before, plus using the above
		osg::Vec4Array* colors = new osg::Vec4Array;
		colors->push_back(osg::Vec4(1.0f, 1.0f, 0.0f, 1.0f));
		linesGeom->setColorArray(colors, osg::Array::BIND_OVERALL);

		// Set the normal in the same way as the color (see note at POINTS, above).
		//osg::Vec3Array* normals = new osg::Vec3Array;
		//normals->push_back(osg::Vec3(0.0f, -1.0f, 0.0f));
		//linesGeom->setNormalArray(normals, osg::Array::BIND_OVERALL);


		// This time we simply use primitive, and hardwire the number of coords to use
		// since we know up front,
		linesGeom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_STRIP, 0, points.size()));

		osg::ref_ptr<osg::LineWidth> lineWid = new osg::LineWidth(3.0);
		linesGeom->getOrCreateStateSet()->setAttribute(lineWid);
		return linesGeom;
}

osg::Geometry* createLineSegments(std::vector<osg::Vec3> points)
{
		// create Geometry object to store all the vertices and lines primitive.
		osg::Geometry* linesGeom = new osg::Geometry();

		// this time we'll preallocate the vertex array to the size
		// and then use an iterator to fill in the values, a bit perverse
		// but does demonstrate that we have just a standard std::vector underneath.
		osg::Vec3Array* vertices = new osg::Vec3Array();
		for (size_t i = 0; i < points.size(); i++)
		{
				vertices->push_back(points[i]);
		}

		// pass the created vertex array to the points geometry object.
		linesGeom->setVertexArray(vertices);

		// set the colors as before, plus using the above
		osg::Vec4Array* colors = new osg::Vec4Array;
		colors->push_back(osg::Vec4(1.0f, 1.0f, 0.0f, 1.0f));
		linesGeom->setColorArray(colors, osg::Array::BIND_OVERALL);

		// Set the normal in the same way as the color (see note at POINTS, above).
		//osg::Vec3Array* normals = new osg::Vec3Array;
		//normals->push_back(osg::Vec3(0.0f, -1.0f, 0.0f));
		//linesGeom->setNormalArray(normals, osg::Array::BIND_OVERALL);


		// This time we simply use primitive, and hardwire the number of coords to use
		// since we know up front,
		linesGeom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES, 0, points.size()));

		osg::ref_ptr<osg::LineWidth> lineWid = new osg::LineWidth(3.0);
		linesGeom->getOrCreateStateSet()->setAttribute(lineWid);
		return linesGeom;
}

osg::Geometry* createPoints(std::vector<osg::Vec3> points)
{
		// create Geometry object to store all the vertices and lines primitive.
		osg::Geometry* linesGeom = new osg::Geometry();

		// this time we'll preallocate the vertex array to the size
		// and then use an iterator to fill in the values, a bit perverse
		// but does demonstrate that we have just a standard std::vector underneath.
		osg::Vec3Array* vertices = new osg::Vec3Array();
		for (size_t i = 0; i < points.size(); i++)
		{
				vertices->push_back(points[i]);
		}

		// pass the created vertex array to the points geometry object.
		linesGeom->setVertexArray(vertices);

		// set the colors as before, plus using the above
		osg::Vec4Array* colors = new osg::Vec4Array;
		colors->push_back(osg::Vec4(0.0f, 1.0f, 1.0f, 1.0f));
		linesGeom->setColorArray(colors, osg::Array::BIND_OVERALL);

		// Set the normal in the same way as the color (see note at POINTS, above).
		osg::Vec3Array* normals = new osg::Vec3Array;
	//	//normals->push_back(osg::Vec3(0.0f, -1.0f, 0.0f));
		//linesGeom->setNormalArray(normals, osg::Array::BIND_OVERALL);


		// This time we simply use primitive, and hardwire the number of coords to use
		// since we know up front,
		linesGeom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS, 0, points.size()));
		osg::Point* pointSize = new osg::Point;
		pointSize->setSize(10.0);
		linesGeom->getOrCreateStateSet()->setAttribute(pointSize);


		return linesGeom;
}

void bindColorShader(osg::StateSet* stateSet)
{
		char fragmentShaderSource[] =
				"uniform vec3 color;\n"
				"void main(void) \n"
				"{\n"
				"  gl_FragColor = vec4(color, 1);\n"
				"}\n";

		stateSet->addUniform(new osg::Uniform("color", osg::Vec3(0,0,0)));
		osg::ref_ptr<osg::Program>  program = new osg::Program;
		program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));
		stateSet->setAttribute(program.get(), osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
}

//void ComputeViewshed(double px, double py)
//{
//		int viewshedSize = 32;
//		osg::ref_ptr<osg::Image> colorImage = new osg::Image;
//		colorImage->allocateImage(viewshedSize, viewshedSize, 1, GL_RGB, GL_UNSIGNED_BYTE);
//		double resol = 1.0 / viewshedSize;
//		double totalarea = 0;
//		double skyarea = 0;
//		int npixel = 0;
//		double diagonal = sqrt(0.5 * 0.5 + 0.5 * 0.5);
//		for (unsigned int row = 0; row < viewshedSize; row++)
//		{
//				double y = 0.5 - row * resol - 0.5 * resol;
//				for (unsigned int col = 0; col < viewshedSize; col++)
//				{
//						double x = -0.5 + col * resol + 0.5 * resol;
//						double dist2ori = sqrt(x * x + y * y);
//						if (dist2ori >= 0.5)
//						{
//								continue;
//						}
//						double zenithD = sqrt(x * x + y * y) * 2.0 * 90.0;//in degrees
//						if (zenithD <= 0.000000001)
//								zenithD = 0.000000001;
//						double zenithR = zenithD * 3.1415926 / 180.0;
//						double x2 = 0.0;
//						double y2 = 1.0;
//						double	cosa = (x * x2 + y * y2) / sqrt((x * x + y * y) * (x2 * x2 + y2 * y2));
//						double lon = acos(cosa) * 180.0 / 3.1415926;
//						lon = 360.0 - lon;
//						lon = 1.0 - (lon / 360.0);
//						double lat = dist2ori;
//				}
//		}
//}


osg::BoundingBox MaskBB;

//void PushVec(std::vector<double>& vecArr, const osg::Vec3d& vec)
//{
	//	vecArr.push_back(vec.)
//}

enum ValType
{
		d = 0,
		vec = 1,
		s = 2,
		i = 3,
};


struct OuputVariable
{
public:
		double _d;
		osg::Vec3d _vec;
		std::string _s;
		long _i;
		ValType _type;
		OuputVariable(double val)
		{
				_d = val;
				_type = ValType::d;
		}

		OuputVariable(const osg::Vec3d& val)
		{
				_vec = val;
				_type = ValType::vec;
		}

		OuputVariable(const std::string& val)
		{
				_s = val;
				_type = ValType::s;
		}

		OuputVariable(int val)
		{
				_i = val;
				_type = ValType::i;
		}

		OuputVariable(long val)
		{
				_i = val;
				_type = ValType::i;
		}

		void Out(std::ofstream& ofs)
		{
				if (_type == ValType::d)
				{
						ofs << _d;
				}
				else if (_type == ValType::s)
				{
						ofs << "\"" << _s << "\"";
				}
				else if (_type == ValType::i)
				{
						ofs << _i;
				}
				else if (_type == ValType::vec)
				{
						ofs << _vec.x() << "," << _vec.y() << "," << _vec.z();
				}
		}
};


void Test()
{
		GDAL_DS<float>* kechuangDSM = new GDAL_DS<float>();
		kechuangDSM->open("E:/Data/models2/kechuang_models/kechuang/kechuang_dsm_25cm.tif", GA_Update);
		float* data = kechuangDSM->readData(1);
		float min = 10000;
		float max = 0;
		float nodata = -9999;
		for (size_t i = 0; i < kechuangDSM->slice; i++)
		{
				float val = data[i];
				if (val < 0 || val > 1000)
				{
						val = nodata;
				}
				else
				{
						if (min > val)
								min = val;
						if (max < val)
								max = val;
				}
				data[i] = val;
		}
		kechuangDSM->writeData(1, data, nodata);
		kechuangDSM->close();
		delete kechuangDSM;
		delete[] data;
}


void GenerateDSM()
{
		osg::ref_ptr<osg::Node> node = osgDB::readNodeFile("E:/Data/models2/kechuang_models/kechuang/kechuang.ive");
		osg::ComputeBoundsVisitor visitor;
		node->accept(visitor);
		HeightVisitor heightField;
  heightField.Create(visitor.getBoundingBox(), 0.25);
  node->accept(heightField);
		heightField.Write("E:/Data/models2/kechuang_models/kechuang/kechuang.tif");
}

std::vector<float> LoadRaster(std::string filename)
{
		GDAL_DS<float>* ds = new GDAL_DS<float>();
		ds->open(filename);
		float* data = ds->readData(1);
		std::vector<float> vec;
		vec.resize(ds->slice);
		memcpy(&vec[0], data, ds->slice * sizeof(float));
		delete[] data;
		delete ds;
		return vec;
}

std::vector<SolarRadiation> LoadReference(std::string dir)
{
		std::vector<float> global = LoadRaster(dir + "Global.tif");
		std::vector<float> beam = LoadRaster(dir + "Beam.tif");
		std::vector<float> diffuse = LoadRaster(dir + "Diffuse.tif");
		std::vector<float> reflected = LoadRaster(dir + "Reflected.tif");
		std::vector<SolarRadiation> vec;
		for (size_t i = 0; i < global.size(); i++)
		{
				SolarRadiation rad;
				rad.global = global[i];
				rad.beam = beam[i];
				rad.diffuse = diffuse[i];
				rad.reflected = reflected[i];
				vec.push_back(rad);
		}
		return vec;
}


//"E:/Code/WeihaiStudyArea/WeihaiDSM.tif"
void ShiftRaster(std::string filename)
{
		GDAL_DS<float>* dt = new GDAL_DS<float>();
		dt->open(filename, GA_Update);
		//UTM
		//dt->adfGeoTransform[0] = 422264.88;
		//dt->adfGeoTransform[3] = 4152158.06;
		//Local
		dt->adfGeoTransform[0] = -3574.17;
		dt->adfGeoTransform[3] = -8555.26;
		dt->m_dataset->SetGeoTransform(dt->adfGeoTransform);
		std::string proj = "";
		dt->m_dataset->SetProjection(proj.data());
		dt->close();
}

void GenerateFlatTerrain(std::string templateFile)
{
		GDAL_DS<float>* dt = new GDAL_DS<float>();
		dt->open(templateFile);
		OGREnvelope bb = dt->bound;
		bb.MaxX = bb.MinX + dt->adfGeoTransform[1] * 9.5;
		bb.MinY = bb.MaxY - dt->adfGeoTransform[1] * 9.5;
		dt->crop(bb, "E:/Code/WeihaiStudyArea/Sub.tif");
		delete dt;
		dt = new GDAL_DS<float>();
		dt->open("E:/Code/WeihaiStudyArea/Sub.tif", GA_Update);
		float* data = dt->readData(1);
		for (size_t i = 0; i < dt->slice; i++)
		{
				data[i] = 0;
		}
		dt->writeData(1, data, -9999);
		dt->close();
		delete dt;
}

void Write(std::string templateFile, std::string outfile, float* data)
{
		GDAL_DS<float>* dt = new GDAL_DS<float>();
		dt->open(templateFile);
		dt->close();
		dt->create(outfile);
		dt->writeData(1, data, 0);
		dt->close();
		delete dt;
}

osg::BoundingBoxd LoadBound(std::string filename)
{
		GDAL_DS<float>* dt = new GDAL_DS<float>();
		dt->open(filename);
		double nodata = dt->getNoData(1);
		float* data = dt->readData(1);
		float minz = 1000000;
		float maxz = -1000000;
		for (size_t i = 0; i < dt->slice; i++)
		{
				float val = data[i];
				if (val == nodata)
						continue;
				if (minz > val)
						minz = val;
				if (maxz < val)
						maxz = val;
		}
		osg::BoundingBoxd bb = osg::BoundingBoxd(dt->bound.MinX, dt->bound.MinY, minz, dt->bound.MaxX, dt->bound.MaxY, maxz);
		delete dt;
		return bb;
}

int main(int argc, char** argv)
{
		GDALAllRegister();
		//GenerateFlatTerrain("E:/Code/WeihaiStudyArea/Elevation.tif");
		//ShiftRaster("E:/Code/WeihaiStudyArea/Aspect.tif");
		//ShiftRaster("E:/Code/WeihaiStudyArea/Slope.tif");
		//ShiftRaster("E:/Code/WeihaiStudyArea/Elevation.tif");
		//ShiftRaster("E:/Code/WeihaiStudyArea/Global.tif");
		//ShiftRaster("E:/Code/WeihaiStudyArea/Diffuse.tif");
		//ShiftRaster("E:/Code/WeihaiStudyArea/Reflected.tif");
		//ShiftRaster("E:/Code/WeihaiStudyArea/Beam.tif");
		//Test();
		//GenerateDSM();
		//ModelLoader::ConvertAspectQGIS2Grass("E:/Code/Weihai_DSM_025_StudyArea_Aspect.tif", "E:/Code/Weihai_DSM_025_StudyArea_Aspect_GrassGIS.tif");
		//ComputeViewshed(0, 0);
		std::vector<std::string> tilenames;
		tilenames.push_back("Tile_-002_-016");
		tilenames.push_back("Tile_-002_-017");
		tilenames.push_back("Tile_-003_-016");
		tilenames.push_back("Tile_-003_-017");
		BBWrapper bb = ModelLoader::CalBound("E:/Data/weihai/Data/", tilenames);
		bb = LoadBound("E:/Code/Weihai_DSM_025_StudyArea.tif");
	//	GDAL_DS<float>* dsm025 = new GDAL_DS<float>();
	//	dsm025->open("E:/Code/Weihai_DSM_025.tif");
		//dsm025->crop(bb, "E:/Code/Weihai_DSM_025_StudyArea.tif");
		//delete dsm025;

		//ModelLoader::CopyLeafTiles("E:/Data/weihai/Data/", tilenames, "E:/Code/StudyArea.osgb");

	//	ModelLoader::Test();
	//	ModelLoader::TileBoundary2Shapefile("E:/Data/weihai/Data/","TileMap.shp");
		osg::Vec3 maskBBCenter = osg::Vec3(-3945.2300, -7068.9270, 0);
		osg::Vec2 maskBBSize = osg::Vec2(1000, 1000);
		MaskBB = osg::BoundingBox(maskBBCenter.x() - maskBBSize.x() * 0.5, maskBBCenter.y() - maskBBSize.y() * 0.5, -10000,
				maskBBCenter.x() + maskBBSize.x() * 0.5, maskBBCenter.y() + maskBBSize.y() * 0.5, 10000);

		SolarParam param;
		param.aspect = 270;
		param.slope = 0;
		param.lon = 122.1204;
		param.lat = 37.5131;
		param.day = 183;
		param.time_step = 0.5;
		GrassSolar grassSolar;
		osg::ref_ptr<osg::Image> colorRampImage = osgDB::readImageFile("E:/Code/ColoEsriTemp.png");
	 //std::vector<SunVector> sunVectors =	grassSolar.getSunVectors(param);
		//sunVectors = grassSolar.getSunVectors(param);
		SolarRadiation rad = grassSolar.calculateSolarRadiation(param);
		printf("%f\n", rad.global);
		//std::vector<SolarRadiation> reference = LoadReference("E:/Code/WeihaiStudyArea/");
		//printfVec3(getLightDir(45, 0), getLightDir(45, 0));
		//printfVec3(getLightDir(45, 90), getLightDir(45, 90));
		//printfVec3(getLightDir(45, 180), getLightDir(45, 180));
		//printfVec3(getLightDir(45, 270), getLightDir(45, 270));
		//printf("%f, %f\n", 90.0, GrassSolar::calculateAspect(GrassSolar::solarAngle2Vector(45, 90)));

		//printf("%f, %f\n", 0.0, GrassSolar::calculateAspect(GrassSolar::solarAngle2Vector(45, 0)));
		//printf("%f, %f\n", 90.0, GrassSolar::calculateAspect(GrassSolar::solarAngle2Vector(45, 90)));
		//printf("%f, %f\n", 180.0, GrassSolar::calculateAspect(GrassSolar::solarAngle2Vector(45, 180)));
		//printf("%f, %f\n", 270.0, GrassSolar::calculateAspect(GrassSolar::solarAngle2Vector(45, 270)));
		//printf("%f, %f\n", 300.0, GrassSolar::calculateAspect(GrassSolar::solarAngle2Vector(45, 300.0)));
		//printf("%f, %f\n", 359.0, GrassSolar::calculateAspect(GrassSolar::solarAngle2Vector(45, 359)));

		//param.slope = 0;
		//while (param.slope <= 90)
		//{
		//		printf("%f, %f\n", param.slope, GrassSolar::calculateSlope(GrassSolar::solarAngle2Vector(param.slope, 359)));
		//		param.slope += 10;
		//}

		//printf("%f, %f\n", 85.0, GrassSolar::calculateSlope(GrassSolar::solarAngle2Vector(85.0, 359)));
		//param.slope = 0;

		//while (param.slope <= 90)
		//{
		//		rad = grassSolar.calculateSolarRadiation(param);
		//		printf("%f,%f\n", param.slope, rad.global);
		//		param.slope += 10;
		//}

		//while (param.aspect <= 360)
		//{
		//		rad = grassSolar.calculateSolarRadiation(param);
		//		printf("%f,%f\n", param.aspect, rad.global);
		//		param.aspect += 10;
		//}

	 //osg::Vec3 grassLight = GrassSolar::solarAngle2Vector(45, 0);
		//osg::Vec3 myLight = getLightDir(45, 0);
		//grassLight =	GrassSolar::solarAngle2Vector(45, 270);
		//grassLight = GrassSolar::solarAngle2Vector(45, 90);
		//grassLight = GrassSolar::solarAngle2Vector(45, 90);
		osg::DisplaySettings::instance()->setNumOfDatabaseThreadsHint(1);
		osg::DisplaySettings::instance()->setNumOfHttpDatabaseThreadsHint(0);

		std::string infile = argv[1];
		QDir qoutdir = QDir(argv[2]);
		if (!qoutdir.exists())
				qoutdir.mkdir(".");
		//DeleteDirectory(qoutdir.absolutePath());
		std::string outdir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();
		double resol = atof(argv[3]);
		double altAngle = atof(argv[4]);
		double azimuthAngle = atof(argv[5]);
		int texSizeX = 2048;
		int texSizeY = 2048;


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
		//osg::ref_ptr<osg::Node> node = ModelLoader::Load3DTiles(infile.data());
		osg::ref_ptr<osg::Node> node = ModelLoader::Load3DTiles(infile.data(), tilenames, false);
		osg::ref_ptr<osg::Node> studyArea = osgDB::readNodeFile("E:/Code/StudyArea.osgb");
		((osg::Group*)studyArea.get())->addChild(studyArea.get());
		node->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE | osg::StateAttribute::PROTECTED);
		studyArea->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE | osg::StateAttribute::PROTECTED);
		viewer.setUpViewInWindow(100, 100, 1024, 768);
		viewer.realize();
		viewer.setCameraManipulator(new osgGA::TrackballManipulator);
		//viewer.setSceneData(node.get());


		//node->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
		//node->getOrCreateStateSet()->setMode( GL_BLEND, osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);
		osg::ref_ptr<osg::MatrixTransform> sceneNode = new osg::MatrixTransform;
		sceneNode->addChild(node.get());

		//ObliqueCameraWrapper* orthoCamera = new ObliqueCameraWrapper(viewer.getCamera(),bb, 45);
		//viewer.getCamera()->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
		viewer.getCamera()->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		//viewer.getCamera()->setClearColor(osg::Vec4(0, 0, 0, 0));

		osg::ref_ptr<ObliqueCamera> orthoCamera = new ObliqueCamera();
		orthoCamera->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
		orthoCamera->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		orthoCamera->setClearColor(osg::Vec4(0, 0, 0, 0));

		orthoCamera->setReferenceFrame(osg::Transform::ABSOLUTE_RF_INHERIT_VIEWPOINT);
		orthoCamera->setViewport(0, 0, texSizeX, texSizeY);
		orthoCamera->setRenderOrder(osg::Camera::PRE_RENDER);
		orthoCamera->setRenderTargetImplementation(osg::Camera::FRAME_BUFFER_OBJECT);
		orthoCamera->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		//viewer.setUpViewInWindow(0, 0, 512, 512);
		//orthoCamera->attach(osg::Camera::COLOR_BUFFER0, tex.get());
		//orthoCamera->attach(osg::Camera::COLOR_BUFFER0, img.get());
		orthoCamera->addChild(sceneNode.get());

		//osg::Matrix mat = osg::Matrix::rotate(osg::DegreesToRadians(azimuthAngle), osg::Vec3(0, 0, 1));
		//sceneNode->setMatrix(mat);
		//osg::ComputeBoundsVisitor visitor;
		//sceneNode->accept(visitor);
		//osg::BoundingBox bb = visitor.getBoundingBox();
		osg::Vec3 center = bb.center();
		osg::Vec2 bbSize(bb.xMax() - bb.xMin(), bb.yMax() - bb.yMin());

		//osg::ComputeBoundsVisitor rawvisitor;
		//node->accept(rawvisitor);

#pragma region 

		char shadowVertexShaderSource[] =
				"#extension GL_EXT_gpu_shader4 : enable\n"
				"uniform mat4 osg_ViewMatrixInverse;\n"
				"uniform mat4 osg_ViewMatrix;\n"
				"varying vec4 pos;\n"
				"varying vec4 normal;\n"
				"void main(void)\n"
				"{\n"
				"gl_TexCoord[0] = gl_MultiTexCoord0;\n"
				"pos =  gl_Vertex;\n"
				"normal = vec4(gl_Normal, 1);\n"
				"gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\n"
				"}\n";

	char shadowFragmentShaderSource[] =
			"#extension GL_EXT_gpu_shader4 : enable\n"
			"uniform sampler2D tex;\n"
			"uniform vec3 minbound;\n"
			"uniform vec3 maxbound;\n"
			"varying vec4 pos;\n"
			"varying vec4 normal;\n"
			"void main(void) \n"
			"{\n"
			"  vec4 color = texture2D(tex, gl_TexCoord[0].xy); \n"
			"  gl_FragData[0] = pos;\n"
			"  //gl_FragData[0] = normal;\n"
			"}\n";

	char vertexShaderSource[] =
			"#extension GL_EXT_gpu_shader4 : enable\n"
			"varying vec4 pos;\n"
			"void main(void)\n"
			"{\n"
			"gl_TexCoord[0] = gl_MultiTexCoord0;\n"
			"pos =  gl_Vertex;\n"
			"gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\n"
			"}\n";

	char fragmentShaderSource[] =
			"#extension GL_EXT_gpu_shader4 : enable\n"
			"uniform sampler2D tex;\n"
			"uniform sampler2D shadowTex;\n"
			"uniform vec3 minbound;\n"
			"uniform vec3 maxbound;\n"
			"uniform vec3 lightPos;\n"
			"uniform mat4 lightMatrix;\n"
			"varying vec4 pos;\n"

			"vec2 poissonDisk[4] = vec2[](\n"
			"vec2(-0.94201624, -0.39906216), \n"
			"vec2(0.94558609, -0.76890725), \n"
			"vec2(-0.094184101, -0.92938870), \n"
			"vec2(0.34495938, 0.29387760)\n"
			"); \n"

			"void main(void) \n"
			"{\n"
			"  vec3 color = texture2D(tex, gl_TexCoord[0].xy); \n"
			"  vec4 shadowUV = lightMatrix * pos; \n"
			"  if(pos.x > minbound.x && pos.x < maxbound.x && pos.y > minbound.y && pos.y < maxbound.y)\n"
			"  {\n"
			"    color = color + vec3(0.3,0,0);\n"
			"    float visibility = 1.0;\n"
			"    //for (int i = 0; i < 4; i++) \n"
			"    //{\n"
			"      //vec4 shadowPos = texture2D(shadowTex, shadowUV.xy + poissonDisk[i] / 700.0); \n"
			"      //float shadowToCamera = length(shadowPos.xyz-lightPos); \n"
			"      //float distToCamera = length(pos.xyz-lightPos); \n"
			"    	 //if(distToCamera - shadowToCamera > 0.2) \n"
			"    	 //{\n"
			"    	 //	 visibility -= 0.2;\n"
			"     	//}\n"
			"   //}\n"
			"      vec4 shadowPos = texture2D(shadowTex, shadowUV.xy); \n"
			"      float shadowToCamera = length(shadowPos.xyz-lightPos); \n"
			"      float distToCamera = length(pos.xyz-lightPos); \n"
			"    	 if(distToCamera - shadowToCamera > 0.2) \n"
			"    	 {\n"
			"         color = color * 0.6;\n"
			"     	}\n"
			"  }\n"
			"  gl_FragColor = vec4(color, 1);\n"
			"}\n";

	osg::Vec3 minbound = osg::Vec3(bb.xMin(), bb.yMin(), bb.zMin());
	osg::Vec3 maxbound = osg::Vec3(bb.xMax(), bb.yMax(), bb.zMax());
	sceneNode->getOrCreateStateSet()->addUniform(new osg::Uniform("tex", 0));
	sceneNode->getOrCreateStateSet()->addUniform(new osg::Uniform("minbound", minbound));
	sceneNode->getOrCreateStateSet()->addUniform(new osg::Uniform("maxbound", maxbound));
	osg::ref_ptr<osg::Program> program = new osg::Program;
	program->addShader(new osg::Shader(osg::Shader::FRAGMENT, shadowFragmentShaderSource));
	program->addShader(new osg::Shader(osg::Shader::VERTEX, shadowVertexShaderSource));
	sceneNode->getOrCreateStateSet()->setAttribute(program.get(), osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
	osg::ClampColor* clamp = new osg::ClampColor();
	clamp->setClampVertexColor(GL_FALSE);
	clamp->setClampFragmentColor(GL_FALSE);
	clamp->setClampReadColor(GL_FALSE);
	sceneNode->getOrCreateStateSet()->setAttribute(clamp, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
#pragma endregion

		osg::ref_ptr<osg::Texture2D> tex = new osg::Texture2D;
		tex->setTextureSize(texSizeX, texSizeY);
		tex->setResizeNonPowerOfTwoHint(false);
		tex->setFilter(osg::Texture2D::MIN_FILTER, osg::Texture2D::LINEAR);
		tex->setFilter(osg::Texture2D::MAG_FILTER, osg::Texture2D::LINEAR);
		tex->setWrap(osg::Texture2D::WRAP_S, osg::Texture2D::REPEAT);
		tex->setWrap(osg::Texture2D::WRAP_T, osg::Texture2D::REPEAT);
		//rtTexture->setDataVariance(osg::Object::DYNAMIC);
		//tex->setInternalFormat(GL_RGBA);
		tex->setSourceFormat(GL_RGBA);
		tex->setInternalFormat(GL_RGBA32F_ARB);
		tex->setSourceType(GL_FLOAT);
		//tex->setSourceType(GL_UNSIGNED_BYTE);
		osg::ref_ptr<osg::Image> img = new osg::Image;
		img->allocateImage(texSizeX, texSizeY,1,GL_RGBA,GL_FLOAT);
		//img->allocateImage(texSizeX, texSizeY, 1, GL_RGBA, GL_UNSIGNED_BYTE);
		tex->setImage(img.get());
		orthoCamera->attach(osg::Camera::COLOR_BUFFER0, tex.get());
		//orthoCamera->attach(osg::Camera::COLOR_BUFFER0, img.get());

		//viewer.setUpViewInWindow(100, 100, 1024, 768);
		viewer.realize();
		viewer.setCameraManipulator(new osgGA::TrackballManipulator);
		osg::ref_ptr<osg::Group> nodeWrapper = new osg::Group;
		nodeWrapper->addChild(node.get());
		program = new osg::Program;
		program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));
		program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexShaderSource));
		//nodeWrapper->getOrCreateStateSet()->setAttribute(program.get(), osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
		sceneNode->getOrCreateStateSet()->addUniform(new osg::Uniform("minbound", minbound));
		sceneNode->getOrCreateStateSet()->addUniform(new osg::Uniform("maxbound", maxbound));
		sceneNode->getOrCreateStateSet()->addUniform(new osg::Uniform("tex", 0));
		//nodeWrapper->getOrCreateStateSet()->addUniform(new osg::Uniform("lightMatrix", osg::Matrix::identity()));
		nodeWrapper->getOrCreateStateSet()->addUniform(new osg::Uniform(osg::Uniform::FLOAT_MAT4, "lightMatrix"));
		nodeWrapper->getOrCreateStateSet()->addUniform(new osg::Uniform("lightPos", osg::Vec3(0,0,0)));
		nodeWrapper->getOrCreateStateSet()->addUniform(new osg::Uniform("minbound", minbound));
		nodeWrapper->getOrCreateStateSet()->addUniform(new osg::Uniform("maxbound", maxbound));
		nodeWrapper->getOrCreateStateSet()->addUniform(new osg::Uniform("tex", 0));
		nodeWrapper->getOrCreateStateSet()->addUniform(new osg::Uniform("shadowTex", 1));
		nodeWrapper->getOrCreateStateSet()->setTextureAttribute(1, tex.get());
		//osg::ref_ptr<osg::Texture2D> tex2 = new osg::Texture2D;
		//osg::ref_ptr<osg::Image> img2 = osgDB::readImageFile("ortho.png");
		//tex2->setImage(img2.get());
		//sceneNode->getOrCreateStateSet()->addUniform(new osg::Uniform("tex", 2));

		//HeightVisitor heightField;
		//heightField.Create(rawvisitor.getBoundingBox(), 10);
		//node->accept(heightField);

		osgDB::DatabasePager* pager = viewer.getScene()->getDatabasePager();
		osg::ref_ptr<osg::Group> root = new osg::Group;
		root->addChild(nodeWrapper.get());
		root->addChild(orthoCamera.get());
		viewer.addEventHandler(new MyKeyboardHandler());
		viewer.setSceneData(root.get());
		/*center = osg::Vec3(-3945.2300, -7068.9270, 0);
		bbSize = osg::Vec2(1000, 1000);
		double minz, maxz;
		BBWrapper localbb = BBWrapper(center.x() - bbSize.x() * 0.5, center.y() - bbSize.y() * 0.5, bb.zMin(),
				center.x() + bbSize.x() * 0.5, center.y() + bbSize.y() * 0.5, bb.zMax());
		heightField.GetHeightRange(localbb, minz, maxz);
		localbb = BBWrapper(center.x() - bbSize.x() * 0.5, center.y() - bbSize.y() * 0.5, minz,
				center.x() + bbSize.x() * 0.5, center.y() + bbSize.y() * 0.5, maxz);
		double scale = 1.0;
		minbound = osg::Vec3(center.x() - localbb.xhalfsize()*scale, center.y() - localbb.yhalfsize()*scale, center.z() - localbb.zhalfsize());
		maxbound = osg::Vec3(center.x() + localbb.xhalfsize()*scale, center.y() + localbb.yhalfsize()*scale, center.z() + localbb.zhalfsize());
		nodeWrapper->getOrCreateStateSet()->getUniform("minbound")->set(minbound);
		nodeWrapper->getOrCreateStateSet()->getUniform("maxbound")->set(maxbound);
		sceneNode->getOrCreateStateSet()->getUniform("minbound")->set(minbound);
		sceneNode->getOrCreateStateSet()->getUniform("maxbound")->set(maxbound);*/
		//return viewer.run();
	//	ModelLoader::CopyLeafTiles(infile, "E:/Code/WeihaiLeafTiles/", localbb);
		//ModelLoader::Test();
		//BBWrapper localbb = bb;
		osg::Vec3 startPoint = bb._min;
		osg::Vec3 endPoint = bb._max;
		startPoint.z() = endPoint.z() = center.z();
		osg::Vec3 startEnd = endPoint - startPoint;
		startEnd.normalize();
		float dist = (endPoint - startPoint).length();
		float curDist = 0;

		orthoCamera->reset2(bb, 90, 0);
		viewer.frame();
		while (pager->getRequestsInProgress())
		{
				viewer.frame();
		}

		//heightField.Create(bb, 2);
		//sceneNode->accept(heightField);
		//heightField.Write("E:/Code/ljm.tif");
		osg::ref_ptr<osg::Geode> shapes = new osg::Geode;
		std::vector<osg::Vec3> points;
		std::vector<osgUtil::LineSegmentIntersector::Intersection> intersections;
		std::vector<SolarRadiation> rads;
		//osg::Vec3d lightdir = getLightDir(10, 90);
		//while (curDist <= dist)
		//{
		//		osg::Vec3 curPoint = startPoint + startEnd * curDist;
		//		osg::Vec3 start = curPoint + lightdir * 1000;
		//		osg::Vec3 end = curPoint - lightdir * 1000;
		//		//osg::ref_ptr<osgUtil::LineSegmentIntersector> intersector = new osgUtil::LineSegmentIntersector(curPoint + osg::Vec3(0, 0, 1000), curPoint - osg::Vec3(0, 0, 1000));
		//		osg::ref_ptr<osgUtil::LineSegmentIntersector> intersector = new osgUtil::LineSegmentIntersector(start, end);
		//		osgUtil::IntersectionVisitor intersectVisitor(intersector.get());
		//		sceneNode->accept(intersectVisitor);

		//		if (intersector->containsIntersections())
		//		{
		//				//for (osgUtil::LineSegmentIntersector::Intersections::iterator hitr = intersector->getIntersections().begin();
		//				//		hitr != intersector->getIntersections().end();
		//				//		++hitr)
		//				//{
		//				//		osg::TessellationHints* hints = new osg::TessellationHints;
		//				//		hints->setDetailRatio(0.5f);
		//				//		//	shapes->addDrawable(new osg::ShapeDrawable(new osg::Sphere(intersector->getFirstIntersection().getWorldIntersectPoint(), 1), hints));
		//				//		printfVec3(hitr->getWorldIntersectPoint(), hitr->getWorldIntersectNormal());
		//				//		points.push_back(hitr->getWorldIntersectPoint());
		//				//}
		//				osg::TessellationHints* hints = new osg::TessellationHints;
		//				hints->setDetailRatio(0.5f);
		//				//shapes->addDrawable(new osg::ShapeDrawable(new osg::Sphere(intersector->getFirstIntersection().getWorldIntersectPoint(), 1), hints));
		//				printfVec3(intersector->getFirstIntersection().getWorldIntersectPoint(), intersector->getFirstIntersection().getWorldIntersectNormal());
		//				points.push_back(intersector->getFirstIntersection().getWorldIntersectPoint());
		//				intersections.push_back(intersector->getFirstIntersection());
		//				rads.push_back(grassSolar.calculateSolarRadiation(param, sceneNode, intersector->getFirstIntersection()));
		//		}
		//		curDist += 20;
		//}

		ShadowCaster shadowCaster;
		std::vector<std::string> rasterFiles;
		std::vector<osg::BoundingBoxd> bounds;
		rasterFiles.push_back("E:/Code/Weihai_DSM_025_StudyArea.tif");
		//rasterFiles.push_back("E:/Code/Weihai_DSM_1m.tif");
		bounds.push_back(bb);
		//bounds.push_back(osg::BoundingBoxd());
		shadowCaster.Initialize(rasterFiles, bounds, "E:/Code/WeihaiStudyArea/Slope.tif", "E:/Code/WeihaiStudyArea/Aspect.tif");



		std::vector <osg::Vec3> linePoints;
	// std::ofstream ofs("e:/compare_outliers.csv");
		bool fixedslope = false;
		std::string outname = "e:/compare_output_025";
		if (fixedslope)
		{
				outname += "_fixedslope";
		}
		outname += ".csv";
		std::ofstream ofs(outname);
		if (!fixedslope)
		{
				ofs << "day,x2d,y2d,z2d,x3d,y3d,z3d,row,col,slope2d,slope3d,aspect2d,aspect3d,global2d,global3d,beam2d,beam3d,dif2d,dif3d,shadow2d,shadow3d\n";
		}
		else
		{
				ofs << "day,x2d,y2d,z2d,x3d,y3d,z3d,row,col,global2d,global3d,beam2d,beam3d,dif2d,dif3d,shadow2d,shadow3d\n";
		}
		srand(time(NULL));
		long n = 0;
		GDAL_DS<float>* ds = (GDAL_DS<float>*)shadowCaster.m_rasters[0];
		double zmax = shadowCaster.m_bounds[0].zMax();
		param.linke = 3.0;
		param.time_step = 0.5;
		param.lon = 122.1204;
		param.lat = 37.5131;
		param.day = 183;
		param.slope = 45;
		param.aspect = 270;
		double elevatedHeight = 0.1;

		std::vector<float> globalRads;
		std::vector<float> beamRads;
		std::vector<float> diffuseRads;
		globalRads.resize(ds->slice);
		beamRads.resize(ds->slice);
		diffuseRads.resize(ds->slice);

		param.time_step= 1;
		float lat = ds->bound.MaxY - ds->adfGeoTransform[1] * 0.5;
		int indx = 0;

		int startRow = 3000;
		int startCol = 1000;
		int rows = 200;
		int cols = 200;
		int count = 0;
		for (int i = startRow; i < startRow + rows; i++)
		{
				lat = ds->bound.MaxY - ds->adfGeoTransform[1] * i - ds->adfGeoTransform[1] * 0.5;
				float lon = ds->bound.MinX + ds->adfGeoTransform[1] * 0.5;
				for (int j = startCol; j < startCol + cols; j++)
				{
						lon = ds->bound.MinX + ds->adfGeoTransform[1] * j + ds->adfGeoTransform[1] * 0.5;
						indx = j + i * ds->ncols;
						float elev = ds->m_cache[indx];
						param.slope = shadowCaster.m_slopeData[indx];
						param.aspect = shadowCaster.m_aspectData[indx];
						SolarRadiation rad2d;
						if (!isnan(param.slope) && !isnan(param.aspect) && !isnan(elev))
						{
								osg::Vec3d intersectPos2d(lon, lat, elev + 0.1);
								param.elev = intersectPos2d.z();
							 rad2d = grassSolar.calculateSolarRadiation(param, &shadowCaster, intersectPos2d);
						}
						else
						{
								rad2d.global = 0;
								rad2d.beam = 0;
								rad2d.diffuse = 0;
						}
						globalRads[indx] = rad2d.global;
						beamRads[indx] = rad2d.beam;
						diffuseRads[indx] = rad2d.diffuse;
						count++;
						if (rad2d.beam == 0) 
						{
								printf("%d/%d,%f\n", count, rows*cols, rad2d.beam);
						}
						printf("%d/%d,%f\n", count, rows*cols,rad2d.beam);
				}
		}
		Write("E:/Code/Weihai_DSM_025_StudyArea.tif", "E:/Code/Weihai_Global.tif", &globalRads[0]);
		Write("E:/Code/Weihai_DSM_025_StudyArea.tif", "E:/Code/Weihai_Beam.tif", &beamRads[0]);
		Write("E:/Code/Weihai_DSM_025_StudyArea.tif", "E:/Code/Weihai_Diffuse.tif", &diffuseRads[0]);

		while (n < 1000)
		{
				int col, row, index;
				double xrand = rand() / (double)RAND_MAX;
				double yrand = rand() / (double)RAND_MAX;
				double x = bb.xMin() + bb.xsize() * xrand;
				double y = bb.yMin() + bb.ysize() * xrand;
				std::tie(col, row, index) =	ds->getCellIndex(x, y);
				if (col == 0 || row == 0 || col == ds->ncols - 1 || row == ds->nrows - 1)
						continue;

				//if (index >= reference.size() || isnan(reference[index].global))
					//	continue;

				//param.day = rand() / 365 + 1;
				param.day = 183;
				//BBWrapper bb;
				//double cellSize;
				//std::tie(col, row, index, dsIndex, cellSize, bb) = shadowCaster.GetCellIndex(x, y);
				double z = zmax;
				//osg::Vec3d lightDir = getLightDir(45, 0);
				osg::Vec3d end = osg::Vec3d(x, y, -1000);
				osg::Vec3d	start = osg::Vec3d(x, y, 1000);
				Ray ray(start, end - start);
				double intersectDist, slope2d, aspect2d, slope3d, aspect3d;
				bool intersects = shadowCaster.Intersects(ray, intersectDist, slope2d, aspect2d);
				osg::Vec3d intersectPos2d = ray.orig + ray.dir * intersectDist;
				osg::ref_ptr<osgUtil::LineSegmentIntersector> intersector = new osgUtil::LineSegmentIntersector(start, end);
				osgUtil::IntersectionVisitor intersectVisitor(intersector.get());
				sceneNode->accept(intersectVisitor);
				if (!intersector->containsIntersections())
						continue;

				std::vector<OuputVariable> outputVariables;
				std::vector<osg::Vec3d> results;

				if (!fixedslope)
				{
						param.slope = slope2d;
						param.aspect = aspect2d;
				}

				intersectPos2d = intersectPos2d + osg::Vec3d(0, 0, elevatedHeight);
				param.elev = intersectPos2d.z();

				SolarRadiation rad2d = grassSolar.calculateSolarRadiation(param, &shadowCaster, intersectPos2d);
				osg::Vec3d intersectPos3d = intersector->getFirstIntersection().getWorldIntersectPoint();
				osg::Vec3d intersectNormal = intersector->getFirstIntersection().getWorldIntersectNormal();
				intersectNormal.normalize();
				slope3d = GrassSolar::calculateSlope(intersectNormal);
				aspect3d = GrassSolar::calculateAspect(intersectNormal);
				if (!fixedslope)
				{
						param.slope = slope3d;
						param.aspect = aspect3d;
				}
				intersectPos3d = intersectPos3d + osg::Vec3d(0, 0, elevatedHeight);
				param.elev = intersectPos3d.z();
				//grassSolar.calculateSolarRadiation(param, sceneNode, &shadowCaster, intersectPos2 + osg::Vec3d(0, 0, 2), linePoints);
				SolarRadiation rad3d = grassSolar.calculateSolarRadiation(param, sceneNode, intersectPos3d);
				double radDiff = abs(rad3d.global - rad2d.global) / rad3d.global * 100;
				printf("%d:%d\n", n, (int)radDiff);
				//SolarRadiation radRef = reference[index];
				//ofs << "day,x2d,y2d,z2d,x3d,y3d,z3d,row,col,slope2d,slope3d,aspect2d,aspect3d,global2d,global3d,beam2d,beam3d,dif2d,dif3d,shadow2d,shadow3d\n";
				outputVariables.push_back(OuputVariable(param.day));
				outputVariables.push_back(OuputVariable(intersectPos2d));
				outputVariables.push_back(OuputVariable(intersectPos3d));
				outputVariables.push_back(OuputVariable(row));
				outputVariables.push_back(OuputVariable(col));
				if (!ofs)
				{
						outputVariables.push_back(OuputVariable(slope2d));
						outputVariables.push_back(OuputVariable(slope3d));
						outputVariables.push_back(OuputVariable(aspect2d));
						outputVariables.push_back(OuputVariable(aspect3d));
				}
				//outputVariables.push_back(OuputVariable(radRef.global));
				outputVariables.push_back(OuputVariable(rad2d.global));
				outputVariables.push_back(OuputVariable(rad3d.global));
				//outputVariables.push_back(OuputVariable(radRef.beam));
				outputVariables.push_back(OuputVariable(rad2d.beam));
				outputVariables.push_back(OuputVariable(rad3d.beam));
				//outputVariables.push_back(OuputVariable(radRef.diffuse));
				outputVariables.push_back(OuputVariable(rad2d.diffuse));
				outputVariables.push_back(OuputVariable(rad3d.diffuse));
				outputVariables.push_back(OuputVariable(rad2d.shadowMasks));
				outputVariables.push_back(OuputVariable(rad3d.shadowMasks));

				for (int v = 0; v < outputVariables.size(); v++)
				{
						outputVariables[v].Out(ofs);
						if (v != outputVariables.size())
						{
								ofs << ",";
						}
				}
				ofs << "\n";
				//double diff = abs(intersectPos.z() - intersectPos2.z()) / intersectPos2.z() * 100;
				//if (diff > 50)
				//{
				//		linePoints.push_back(start);
				//		linePoints.push_back(intersectPos2);
				//		linePoints.push_back(start);
				//		linePoints.push_back(intersectPos);
				//		n++;
				//	 	ofs << start.x() << "," << start.y() << "," << start.z() << "," << intersectPos.x() << "," << intersectPos.y() << "," << intersectPos.z() << "," << intersectPos2.x() << "," << intersectPos2.y() << "," << intersectPos2.z() << "\n";
				//}
				//osg::TessellationHints* hints = new osg::TessellationHints;
				//hints->setDetailRatio(0.5f);
				//shapes->addDrawable(new osg::ShapeDrawable(new osg::Sphere(intersector->getFirstIntersection().getWorldIntersectPoint(), 1), hints));
				//printfVec3(intersector->getFirstIntersection().getWorldIntersectPoint(), intersector->getFirstIntersection().getWorldIntersectNormal());
				//printfVec3(intersector->getFirstIntersection().getWorldIntersectPoint(), intersectPos);
				//points.push_back(intersector->getFirstIntersection().getWorldIntersectPoint());
				//intersections.push_back(intersector->getFirstIntersection());
				//rads.push_back(grassSolar.calculateSolarRadiation(param, sceneNode, intersector->getFirstIntersection()));
				n++;
				if(n % 1000 == 0)
						ofs.flush();
		}
		ofs.flush();
		ofs.close();
		bindColorShader(shapes->getOrCreateStateSet());
		//osg::ref_ptr<osg::Geometry> line = createLine(points);
		//osg::ref_ptr<osg::Geometry> pointgeo = createPoints(points);
		//shapes->addDrawable(line.get());
		shapes->addDrawable(createLineSegments(linePoints));
		root->addChild(shapes.get());
		//for (size_t n = 0; n < 200; n++)
//{
//		double xrand = rand() / (double)RAND_MAX;
//		double yrand = rand() / (double)RAND_MAX;
//		double x = bb.xMin() + bb.xsize() * xrand;
//		double y = bb.yMin() + bb.ysize() * xrand;
//		int col, row, index, dsIndex;
//		BBWrapper bb;
//		double cellSize;
//		std::tie(col, row, index, dsIndex, cellSize, bb) = shadowCaster.GetCellIndex(x, y);
//		double z = bb.zMax() + 2;
//		osg::Vec3d start(x, y, z);
//		osg::Vec3d	end = start + getLightDir(45, 0) * 1000;
//		osg::Vec3d lightDir = getLightDir(45, 0);
//		Ray ray(start, end - start);
//		double intersectDist;
//		bool intersects = shadowCaster.Intersects(ray, intersectDist);
//		osg::Vec3d intersectPos = ray.orig + ray.dir * intersectDist;
//		osg::ref_ptr<osgUtil::LineSegmentIntersector> intersector = new osgUtil::LineSegmentIntersector(start, end);
//		osgUtil::IntersectionVisitor intersectVisitor(intersector.get());
//		sceneNode->accept(intersectVisitor);
//		if (intersector->containsIntersections())
//		{
//				osg::Vec3d intersectPos2 = intersector->getFirstIntersection().getWorldIntersectPoint();
//				//osg::TessellationHints* hints = new osg::TessellationHints;
//			 //hints->setDetailRatio(0.5f);
//				//shapes->addDrawable(new osg::ShapeDrawable(new osg::Sphere(intersector->getFirstIntersection().getWorldIntersectPoint(), 1), hints));
//				//printfVec3(intersector->getFirstIntersection().getWorldIntersectPoint(), intersector->getFirstIntersection().getWorldIntersectNormal());
//				//printfVec3(intersector->getFirstIntersection().getWorldIntersectPoint(), intersectPos);
//			 //points.push_back(intersector->getFirstIntersection().getWorldIntersectPoint());
//				//intersections.push_back(intersector->getFirstIntersection());
//				//rads.push_back(grassSolar.calculateSolarRadiation(param, sceneNode, intersector->getFirstIntersection()));
//				ofs << intersectPos.x() << "," << intersectPos.y() << "," << intersectPos.z() << "," << intersectPos2.x() << "," << intersectPos2.y() << "," << intersectPos2.z() << "\n";
//		}
//		else
//		{
//				printf("%s,%s\n", intersects ? "true" : "false", "false");
//		}
//}
		//double y = minbound.y();
		//while (y <= maxbound.y())
		//{
		//		double x = minbound.x();
		//		while (x <= maxbound.x())
		//		{
		//				osg::Vec3d end(x, y, minbound.z() - 1000);
		//				osg::Vec3d start(x, y, maxbound.z() + 1000);
		//				start = osg::Vec3d(-3221.2, -8797.5, 30);
		//				end = start + getLightDir(45, 0) * 1000;
		//					
		//				Ray ray(start, end - start);
		//				double intersectDist;
		//			 bool intersects =	shadowCaster.Intersects(ray, intersectDist);
		//				osg::Vec3d intersectPos = ray.orig + ray.dir * intersectDist;
		//				printfVec3(ray.orig, intersectPos);
		//				//osg::ref_ptr<osgUtil::LineSegmentIntersector> intersector = new osgUtil::LineSegmentIntersector(curPoint + osg::Vec3(0, 0, 1000), curPoint - osg::Vec3(0, 0, 1000));
		//				osg::ref_ptr<osgUtil::LineSegmentIntersector> intersector = new osgUtil::LineSegmentIntersector(start, end);
		//				osgUtil::IntersectionVisitor intersectVisitor(intersector.get());
		//				sceneNode->accept(intersectVisitor);
		//				if (intersector->containsIntersections())
		//				{
		//						osg::TessellationHints* hints = new osg::TessellationHints;
		//						hints->setDetailRatio(0.5f);
		//						//shapes->addDrawable(new osg::ShapeDrawable(new osg::Sphere(intersector->getFirstIntersection().getWorldIntersectPoint(), 1), hints));
		//						//printfVec3(intersector->getFirstIntersection().getWorldIntersectPoint(), intersector->getFirstIntersection().getWorldIntersectNormal());
		//						printfVec3(intersector->getFirstIntersection().getWorldIntersectPoint(), intersectPos);
		//						points.push_back(intersector->getFirstIntersection().getWorldIntersectPoint());
		//						intersections.push_back(intersector->getFirstIntersection());
		//						rads.push_back(grassSolar.calculateSolarRadiation(param, sceneNode, intersector->getFirstIntersection()));
		//				}
		//				x += 50;
		//		}
		//		y += 50;
		//}

		return viewer.run();

		double maxRad = -9999999;
		double minRad = 9999999;
		for (long i = 0; i < rads.size(); i++)
		{
				if (maxRad < rads[i].global)
						maxRad = rads[i].global;
				if (minRad > rads[i].global)
						minRad = rads[i].global;
		}

		osg::ref_ptr<osg::TessellationHints> hints = new osg::TessellationHints;
		hints->setDetailRatio(0.5f);
		osg::Vec3 startColor = osg::Vec3(0, 1, 0);
		osg::Vec3 endColor = osg::Vec3(1, 0, 0);
		for (long i = 0; i < rads.size(); i++)
		{
				float linearPos = (rads[i].global - minRad) / maxRad;
				//osg::Vec4 color = colorRampImage->getColor((unsigned int)(linearPos * colorRampImage->s()), colorRampImage->t() / 2);

				/*osg::ref_ptr<osg::Material> mat = new osg::Material;
				mat->setColorMode(osg::Material::DIFFUSE);
				mat->setAmbient(osg::Material::FRONT_AND_BACK, osg::Vec4(0.2, 0.2, 0.2, 1.0));
				mat->setDiffuse(osg::Material::FRONT_AND_BACK, color);
				mat->setSpecular(osg::Material::FRONT_AND_BACK, osg::Vec4(1.0, 1.0, 1.0, 1.0));
				mat->setShininess(osg::Material::FRONT_AND_BACK, 64);*/
				osg::Vec3 color = startColor + (endColor - startColor) * linearPos;
				osg::ref_ptr< osg::ShapeDrawable> shapeDrawable = new osg::ShapeDrawable(new osg::Sphere(intersections[i].getWorldIntersectPoint(), 3), hints.get());
				//shapeDrawable->getOrCreateStateSet()->setAttributeAndModes(mat.get(), osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
				//shapeDrawable->getOrCreateStateSet()->addUniform(new osg::Uniform("color", osg::Vec3(color.r(), color.g(), color.b())));
				shapeDrawable->getOrCreateStateSet()->addUniform(new osg::Uniform("color", color));
				shapes->addDrawable(shapeDrawable.get());
		}
		bindColorShader(shapes->getOrCreateStateSet());
		//osg::ref_ptr<osg::Geometry> line = createLine(points);
		//osg::ref_ptr<osg::Geometry> pointgeo = createPoints(points);
		//shapes->addDrawable(line.get());
		//shapes->addDrawable(pointgeo.get());
		root->addChild(shapes.get());
		altAngle = 45;
		azimuthAngle = 0;

		int frameCount = 0;
		//while (!viewer.done())
		//{
		//		float sum = 0;
		//	//	frameCount++;
		//				//if (frameCount % 25 == 0)
		//		if (_canUpdateLightDir)
		//		{
		//				_canUpdateLightDir = false;
		//				//osgDB::writeImageFile(*img,"ortho.png");
		//				/*
		//				osg::Vec4* pShadowPos = (osg::Vec4*) img->data();
		//				double sum = 0;
		//				for (size_t i = 0; i < img->s() * img->t(); i++)
		//				{
		//						if (isnan(pShadowPos->z()))
		//								continue;
		//						sum += (pShadowPos->z());
		//				}
		//				printf("%f\n", sum/(img->s() * img->t()));*/

		//				//azimuthAngle = 0;
		//				azimuthAngle += 1;
		//				if (azimuthAngle >= 360)
		//				{
		//						altAngle += 1;
		//						azimuthAngle = 0;
		//				}

		//				if (altAngle >= 90)
		//				{
		//						altAngle = 45;
		//				}

		//				orthoCamera->reset2(bb, altAngle, azimuthAngle);

		//				osg::Matrix lightMatrix = mat * orthoCamera->getViewMatrix() * orthoCamera->getProjectionMatrix()
		//						* osg::Matrix::scale(0.5, 0.5, 0.5)* osg::Matrix::translate(0.5, 0.5, 0.5);

		//				nodeWrapper->getOrCreateStateSet()->getUniform("lightMatrix")->set(lightMatrix);
		//				osg::Vec3 centerUV = orthoCamera->_center * lightMatrix;
		//				centerUV = minbound * lightMatrix;
		//				centerUV = maxbound * lightMatrix;
		//				nodeWrapper->getOrCreateStateSet()->getUniform("lightPos")->set(osg::Vec3(orthoCamera->_eye * mat));
		//		}
		//		//	if (altAngle > 85)
		//		//			altAngle = 15;
		//		//	mat = osg::Matrix::rotate(osg::DegreesToRadians(azimuthAngle), osg::Vec3(0, 0, 1));
		//		//	sceneNode->setMatrix(mat);

		//		/*	osg::BoundingBoxd transformedBB;
		//			transformedBB.init();
		//			for (size_t i = 0; i < 8; i++)
		//			{
		//					transformedBB.expandBy(localbb.corner(i) * mat);
		//			}
		//			localbb = BBWrapper(transformedBB);*/

		//		viewer.frame();
		//}
		//printf("%f\n", 0.0);

		//altAngle = 90;
		//bool isFirst = true;

		//altAngle = 90;
		//while (altAngle > 5 && !viewer.done())
		//{
		//		azimuthAngle = 0;
		//		while (azimuthAngle < 360 && !viewer.done())
		//		{


		//				viewer.frame();
		//				//if (isFirst)
		//				//{
		//				//		while (pager->getRequestsInProgress())
		//				//		{
		//				//				viewer.frame();
		//				//		}
		//				//		isFirst = false;
		//				//}
		//				//viewer.frame();
		//				//orthoCamera->dirtyAttachmentMap();
		//				int count = 0;
		//				/*while (pager->getRequestsInProgress())
		//				{
		//						viewer.frame();
		//						count++;
		//						if (count > 5)
		//								break;
		//				}*/
		//				//std::stringstream ssoutname;
		//				//ssoutname << outdir << altAngle << "_" << azimuthAngle << ".png";
		//				//printf("%s\n", ssoutname.str().data());
		//				//osgDB::writeImageFile(*img, ssoutname.str().data());
		//				azimuthAngle += 12;
		//		}
		//			altAngle -= 3;
		//}

		//for (size_t i = 0; i < imageArr.size(); i++)
		//{
		//		ImageOutput* imageOutput = imageArr[i];
		//		imageOutput->Write(outdir);
		//		delete imageOutput;
		//}


		//GDALDataset *poSrcDS = (GDALDataset *)GDALOpen((outdir + "oblique.vrt").data(), GA_ReadOnly);

		//char **papszOptions = NULL;

		//papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");
		//papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "LZW");
		//GDALDataset * poDstDS = poDriver->CreateCopy((outdir + "oblique.tif").data(), poSrcDS, FALSE,
		//	papszOptions, GDALTermProgress, NULL);
		///* Once we're done, close properly the dataset */
		//if (poDstDS != NULL)
		//	GDALClose((GDALDatasetH)poDstDS);
		//CSLDestroy(papszOptions);
		//GDALClose((GDALDatasetH)poSrcDS);

		//poDriver->Delete((outdir + "oblique.vrt").data());
		//for (size_t i = 0; i < tilefiles.size(); i++)
		//{
		//	if (QFileInfo(tilefiles[i].data()).exists())
		//		QFile::remove(tilefiles[i].data());
		//}

		return viewer.run();
		//exit(0);
}

