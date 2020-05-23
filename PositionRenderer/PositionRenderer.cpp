
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
#include "osg/LineSegment"
#include <osg/ShapeDrawable>
#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include "ogrsf_frmts.h"
#include "ModelLoader.h"
#include "ScreenOverlay.h"
#include "osgUtil/SmoothingVisitor"

using namespace osgDB;
using namespace OpenThreads;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  SortFileRequestFunctor
//

class BBWrapper : public osg::BoundingBoxd
{
public:
		BBWrapper() :
				osg::BoundingBoxd()
		{

		}

		BBWrapper(const osg::BoundingBoxd& bb) :
				osg::BoundingBoxd(bb)
		{

		}

		BBWrapper(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax) :
				osg::BoundingBoxd(xmin, ymin, zmin, xmax, ymax, zmax)
		{

		}

		BBWrapper(osg::Vec3d min, osg::Vec3d max) :
				osg::BoundingBoxd(min, max)
		{

		}

		double xsize() { return xMax() - xMin(); }
		double ysize() { return yMax() - yMin(); }
		double zsize() { return zMax() - zMin(); }
		double xhalfsize() { return xsize() * 0.5; }
		double yhalfsize() { return ysize() * 0.5; }
		double zhalfsize() { return zsize() * 0.5; }
		osg::Vec3d size() { return _max - _min; }
		osg::Vec3d halfsize() { return size() * 0.5; }
};

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
		int m_imageID;
		std::string m_outdir;
public:
		BBWrapper m_BB;

		GeometryVisitor() :osg::NodeVisitor(TRAVERSE_ALL_CHILDREN) { setNodeMaskOverride(0xffffffff); }

		void initialize(std::string outdir)
		{
				m_group = new osg::Group();
				m_imageID = 0;
				m_outdir = outdir;
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
				osg::Texture2D* texture2d = nullptr;
				if (node.getStateSet())
				{
						osg::Texture2D* texture2d = dynamic_cast<osg::Texture2D*>(node.getStateSet()->getTextureAttribute(0, osg::StateAttribute::TEXTURE));
				}

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
								normals = (osg::Vec3Array*) geom->getNormalArray();
						}
						if (!texture2d)
						{
								texture2d = dynamic_cast<osg::Texture2D*>(geom->getStateSet()->getTextureAttribute(0, osg::StateAttribute::TEXTURE));
								//std::string filename = texture->getName();
								if (texture2d)
								{
										osg::Image* img = texture2d->getImage();
										if (img)
										{
												std::string filename = img->getFileName();
												filename = img->getFileName();
												//std::stringstream ss;
												//ss << m_imageID++ << ".jpg";
												osgDB::writeImageFile(*img, m_outdir + filename);
												//img->setFileName(ss.str());
										}
								}
						}
						for (int j = 0; j < geom->getNumPrimitiveSets(); ++j)
						{
								osg::PrimitiveSet* primitiveSet = geom->getPrimitiveSet(j);
								FaceVisitor visitor(primitiveSet);
								//std::vector<unsigned int> indices = visitor.getFaceIndices();

								//for (unsigned int k = 0; k < indices.size(); k++)
								//{
								//		unsigned int index = indices[k];
								//		osg::Vec3 pos = (*vertices)[index];
								//		if (!isIdentify)
								//				pos = pos * matWorld;
								//		osg::Vec3 normal = (*normals)[index];
								//		osg::Vec2 uv = (*uvs)[index];
								//}
						}
				}
		}

		virtual void apply(osg::Node& node)
		{
	
		}
		
		virtual void traverseGeode(osg::Geode* node)
		{
				osg::MatrixList matlist = node->getWorldMatrices();
				osg::Matrix matWorld = osg::Matrix::identity();

				for (unsigned int i = 0; i < matlist.size(); i++)
				{
						matWorld = matlist[i] * matWorld;
				}
				bool isIdentify = (matWorld.isIdentity());
				osg::Texture2D* texture2d = nullptr;
				if (node->getStateSet())
				{
						osg::Texture2D* texture2d = dynamic_cast<osg::Texture2D*>(node->getStateSet()->getTextureAttribute(0, osg::StateAttribute::TEXTURE));
				}
				for (int i = 0; i < node->getNumDrawables(); i++)
				{
						osg::Drawable* drawable = node->getDrawable(i);
						osg::Geometry* geom = dynamic_cast<osg::Geometry*>(drawable);
						if (!geom) continue;
						//osg::ref_ptr<osg::Geometry> geomCpy = new osg::Geometry;

						osg::Vec3Array* vertices = (osg::Vec3Array*) geom->getVertexArray();
						osg::Vec2Array* uvs = (osg::Vec2Array*) geom->getTexCoordArray(0);
						//osg::Vec3Array* normals = (osg::Vec3Array*) geom->getNormalArray();
						//if (!normals)
						//{
						//		osgUtil::SmoothingVisitor::smooth(*geom);
						//		//normals = (osg::Vec3Array*) geom->getNormalArray();
						//}
						if (!texture2d)
						{
								texture2d = dynamic_cast<osg::Texture2D*>(geom->getStateSet()->getTextureAttribute(0, osg::StateAttribute::TEXTURE));
								//std::string filename = texture->getName();
								if (texture2d)
								{
										osg::Image* img = texture2d->getImage();
										if (img)
										{
												std::string filename = img->getFileName();
												filename = img->getFileName();
												//std::stringstream ss;
												//ss << m_imageID++ << ".jpg";
												osgDB::writeImageFile(*img, m_outdir + filename);
												//img->setFileName(ss.str());
										}
								}
						}
						for (int j = 0; j < geom->getNumPrimitiveSets(); ++j)
						{
								osg::PrimitiveSet* primitiveSet = geom->getPrimitiveSet(j);
								FaceVisitor visitor(primitiveSet);
								//std::vector<unsigned int> indices = visitor.getFaceIndices();

								//for (unsigned int k = 0; k < indices.size(); k++)
								//{
								//		unsigned int index = indices[k];
								//		osg::Vec3 pos = (*vertices)[index];
								//		if (!isIdentify)
								//				pos = pos * matWorld;
								//		osg::Vec3 normal = (*normals)[index];
								//		osg::Vec2 uv = (*uvs)[index];
								//}
						}
				}

				std::stringstream ss;
				ss << m_imageID++ << ".osg";
				osgDB::writeNodeFile(*node, m_outdir + ss.str());
				m_group->addChild(node);
		}

		//virtual void traverseGroup(osg::Node* node)
		//{
		//		osg::Group* group = dynamic_cast<osg::Group*>(node);
		//		osg::PagedLOD* pageLOD = dynamic_cast<osg::PagedLOD*>(node);
		//		osg::Geode* geode = dynamic_cast<osg::Geode*>(node);
		//		osg::Geometry* geom = dynamic_cast<osg::Geometry*>(node);
		//		for (size_t i = 0; i < group->getNumChildren(); i++)
		//		{
		//				osg::Node* child = group->getChild(i);
		//				if (m_BB.contains(child->getBound().center()))
		//				{
		//						printf("");
		//				}
		//				std::string name = child->getName();
		//				int startPos = name.size();
		//				while (startPos >= 0)
		//				{
		//						if (name[startPos] == '\\' || name[startPos] == '/' || name[startPos] == '//')
		//						{
		//								break;
		//						}
		//						startPos--;
		//				}
		//				name = name.substr(startPos + 1, name.size() - startPos);
		//				//osgDB::writeNodeFile(*child, name + ".osg");
		//		}
		//}

		virtual void traverseGroup(osg::Node* node)
		{
				osg::Group* group = dynamic_cast<osg::Group*>(node);
				osg::PagedLOD* pageLOD = dynamic_cast<osg::PagedLOD*>(node);
				osg::Geode* geode = dynamic_cast<osg::Geode*>(node);
				osg::Geometry* geom = dynamic_cast<osg::Geometry*>(node);
				if (group && !geode)
				{
						if (group->getStateSet())
						{
								osg::Texture2D* texture2d = dynamic_cast<osg::Texture2D*>(group->getStateSet()->getTextureAttribute(0, osg::StateAttribute::TEXTURE));
								//std::string filename = texture->getName();
								if (texture2d)
								{
										osg::Image* img = texture2d->getImage();
										if (img)
										{
												std::string filename = img->getFileName();
												filename = img->getFileName();
												//std::stringstream ss;
												//ss << m_imageID++ << ".jpg";
												osgDB::writeImageFile(*img, m_outdir + filename);
												//img->setFileName(ss.str());
										}
								}
						}

						for (size_t i = 0; i < group->getNumChildren(); i++)
						{
								traverseGroup(group->getChild(i));
						}
						return;
				}

				if (geode)
				{
						if (m_BB.contains(geode->getBound().center()))
						{
								traverseGeode(geode);
						}
				}
		}
};

void printfVec3(osg::Vec3 p1, osg::Vec3 p2)
{
		printf("(%f,%f,%f),(%f,%f,%f)\n", p1.x(), p1.y(), p1.z(), p2.x(), p2.y(), p2.z());
}

class ObliqueCamera : public osg::Camera
{
public:
		double _altAngle;
		double _azimuthAngle;
		double _viewDistance;
		osg::Vec3d _viewDirection;
		osg::Vec3d _centerEyePoint;
		osg::Vec3d _center;
		osg::Vec3d _upDirection;
		double _znear;
		double _zfar;
		int _ncols;
		int _nrows;
		int imgWidth;
		int imgHeight;
		double _srcTileSizeY;
		double _srcTileSizeX;
		osg::BoundingBox _bb;
		GDALDriver *_poDriver;
		int _curRow;
		int _curCol;
		double _curAdfGeoTransform[6];
		ObliqueCamera()
		{
				this->setReferenceFrame(osg::Transform::ABSOLUTE_RF);
				GDALAllRegister();
				const char *pszFormat = "GTiff";
				_poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
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
		}

		//void reset(const osg::BoundingBox& bb, double altAngle = 0, double azimuthAngle = 0)
		//{
		//		_altAngle = osg::DegreesToRadians(altAngle);
		//		_azimuthAngle = osg::DegreesToRadians(azimuthAngle);
		//		osg::Vec3d bbmin(bb.xMin(), bb.yMin(), bb.zMin());
		//		osg::Vec3d bbmax(bb.xMax(), bb.yMax(), bb.zMin());
		//		osg::Vec3d bbsize = bbmax - bbmin;

		//		float maxsize = bbsize.x() > bbsize.y() ? bbsize.x() : bbsize.y();
		//		_viewDistance = maxsize * 10;

		//		osg::Vec3d translation(0, -_viewDistance * cos(_altAngle), _viewDistance * sin(_altAngle));
		//		if (_altAngle == 0 || osg::RadiansToDegrees(_altAngle) == 90)
		//		{
		//				translation = osg::Vec3d(0, 0, _viewDistance);
		//		}
		//		osg::Vec3d bbminNear = bbmin + translation;
		//		osg::Vec3d bbmaxNear = bbmax + translation;
		//		osg::Vec3d bbminFar = bbmin - translation;
		//		osg::Vec3d bbmaxFar = bbmax - translation;

		//		_znear = _viewDistance - bb.radius() * 2;
		//		_zfar = _viewDistance + bb.radius() * 2;
		//		float  right = (bb.xMax() - bb.xMin())*0.5;



		//		//_upDirection = osg::Vec3d(0.0, 1.0, 0.0);
		//		//_upDirection.normalize();
		//		osg::Vec3d vforward(0, -cos(_altAngle), sin(_altAngle));
		//		osg::Vec3d vright(1, 0, 0);
		//		_upDirection = vforward ^ vright;
		//		//_upDirection = vright ^ vforward;

		//		_viewDirection = -translation;
		//		_viewDirection.normalize();
		//		_center = (bbmin + bbmax) * 0.5;

		//		osg::Vec3d eyePoint = _center - _viewDirection * _viewDistance;
		//		//this->setViewMatrixAsLookAt(eyePoint, _center, _upDirection);

		//		_bb = bb;
		//}

		void setupGrid(double destTileSizeX, double destTileSizeY, int imgW, int imgH)
		{
				imgWidth = imgW;
				imgHeight = imgH;
				double zspan = _bb.zMax() - _bb.zMin();
				double tanAngle = tan(_altAngle);
				double projectedLen = zspan / tanAngle;
				if (_altAngle == 0 || osg::RadiansToDegrees(_altAngle) == 90)
				{
						projectedLen = 0;
						tanAngle = 1;
				}

				_srcTileSizeY = destTileSizeY * tanAngle;
				_srcTileSizeX = destTileSizeX;

				double xspan = _bb.xMax() - _bb.xMin();
				double yspan = _bb.yMax() - _bb.yMin();

				double yspanProjected = yspan + projectedLen;


				_nrows = (int)(yspanProjected / _srcTileSizeY);
				while (_nrows * _srcTileSizeY < yspanProjected)
						_nrows = _nrows + 1;

				_ncols = (int)(xspan / _srcTileSizeX);
				while (_ncols * _srcTileSizeX < xspan)
						_ncols = _ncols + 1;
				_center.y() = _center.y() + projectedLen * 0.5;
				xspan = _ncols * _srcTileSizeX;

				_bb.xMin() = _center.x() - _ncols * _srcTileSizeX * 0.5;
				_bb.xMax() = _center.x() + _ncols * _srcTileSizeX * 0.5;
				_bb.yMin() = _center.y() - _nrows * _srcTileSizeY * 0.5;
				_bb.yMax() = _center.y() + _nrows * _srcTileSizeY * 0.5;
		}
		void setupGrid(double destTileSize, int imgSize)
		{
				return setupGrid(destTileSize, destTileSize, imgSize, imgSize);
			/*	double zspan = _bb.zMax() - _bb.zMin();
				double tanAngle = tan(_altAngle);
				double projectedLen = zspan / tanAngle;
				if (_altAngle == 0 || osg::RadiansToDegrees(_altAngle) == 90)
						tanAngle = 1;
				_srcTileSizeY = destTileSize * tanAngle;
				_srcTileSizeX = destTileSize;
				double xspan = _bb.xMax() - _bb.xMin();
				double yspan = _bb.yMax() - _bb.yMin();

				double yspanProjected = yspan + projectedLen;


				_nrows = (int)(yspanProjected / _srcTileSizeY);
				while (_nrows * _srcTileSizeY < yspanProjected)
						_nrows = _nrows + 1;

				_ncols = (int)(xspan / _srcTileSizeX);
				while (_ncols * _srcTileSizeX < xspan)
						_ncols = _ncols + 1;
				_center.y() = _center.y() + projectedLen * 0.5;
				xspan = _ncols * _srcTileSizeX;
				osg::BoundingBox bb = _bb;
				bb.xMin() = _center.x() - _ncols * _srcTileSizeX * 0.5;
				bb.xMax() = _center.x() + _ncols * _srcTileSizeX * 0.5;
				bb.yMin() = _center.y() - _nrows * _srcTileSizeY * 0.5;
				bb.yMax() = _center.y() + _nrows * _srcTileSizeY * 0.5;

				_bb = bb;*/
		}

		void getRowCol(osg::Vec3d pos, int& nrow, int& ncol)
		{
				osg::Vec3d l1 = pos + _viewDirection * _viewDistance;
				osg::Vec3d l2 = pos - _viewDirection * _viewDistance;
				osg::ref_ptr<osg::LineSegment> line = new osg::LineSegment;
				line->set(l1, l2);
				osg::Vec3d center = _center;
				center.z() = _bb.zMin();

				osg::Vec3d p1 = center + (osg::Vec3d(_bb.xMin(), _bb.yMin(), _bb.zMin()) - _center) * 10;
				osg::Vec3d p2 = center + (osg::Vec3d(_center.x(), _bb.yMax(), _bb.zMin()) - _center) * 10;
				osg::Vec3d p3 = center + (osg::Vec3d(_bb.xMax(), _bb.yMin(), _bb.zMin()) - _center) * 10;
				double ratio;
				line->intersect(p1, p2, p3, ratio);
				pos = line->start() + (line->end() - line->start()) * ratio;

				nrow = (int)((_bb.yMax() - pos.y()) / _srcTileSizeY);
				ncol = (int)((pos.x() - _bb.xMin()) / _srcTileSizeY);
		}
		void setupCameraAtCell(int nrow, int ncol)
		{
				double sinAngle = sin(_altAngle);
				if (_altAngle == 0 || osg::RadiansToDegrees(_altAngle) == 90)
						sinAngle = 1;
				double tanAngle = tan(_altAngle);
				if (_altAngle == 0 || osg::RadiansToDegrees(_altAngle) == 90)
						tanAngle = 1;
				osg::Vec3d center = _center;
				center.x() = _bb.xMin() + ncol * _srcTileSizeX + _srcTileSizeX * 0.5;
				center.y() = _bb.yMax() - nrow * _srcTileSizeY - _srcTileSizeY * 0.5;
				//center.z() = (_bb.zMin() + _bb.zMax()) * 0.5;
				osg::Vec3d eyePoint = center - _viewDirection * _viewDistance;
				float left = -_srcTileSizeX * 0.5;
				float right = _srcTileSizeX * 0.5;

				float bottom = -_srcTileSizeY * sinAngle * 0.5;
				float top = _srcTileSizeY * sinAngle * 0.5;

				//##########################################################
						//center = osg::Vec3d(0, 50, -50);
						//osg::Vec3d forward(0, 0.707106781, -0.707106781);
						//eyePoint = osg::Vec3d(0, -294.9489706, 294.9489706);

						//bottom = -67.67766953;
						//top = 67.67766953;
						//left = -50;
						//right = 50;
						//_znear = 173.2050781;
						//_zfar = 519.6152344;

						//float length = 346.4101563;
						//float dist = 346.4101563;
				//#############################################################
				this->setViewMatrixAsLookAt(eyePoint, center, _upDirection);
				//float shiftdist = (_bb.zMax() - _bb.zMin()) * 0.5;
				this->setProjectionMatrixAsOrtho(left, right, bottom, top, _znear, _zfar);
				_curRow = nrow;
				_curCol = ncol;

				_curAdfGeoTransform[0] = _bb.xMin() + ncol * _srcTileSizeX;///* top left x */
				_curAdfGeoTransform[1] = _srcTileSizeX / imgWidth;///* w-e pixel resolution */
				_curAdfGeoTransform[2] = 0;///* 0 */
				_curAdfGeoTransform[3] = _bb.yMax() - nrow * _srcTileSizeY / tanAngle;// /* top left y */
				_curAdfGeoTransform[4] = 0;///* 0 */
				_curAdfGeoTransform[5] = -_srcTileSizeX / imgHeight;///* n-s pixel resolution (negative value) */
		}
		bool imagExists(int nrow, int ncol, osg::Image* img, std::string dir = "./img/")
		{
				std::stringstream ssDest;
				ssDest << dir << nrow << "_" << ncol << ".tif";
				if (QFileInfo(ssDest.str().data()).exists())
						return true;
				return false;
		}

		std::string produceImage(int nrow, int ncol, osg::Image* img, std::string dir = "./img/")
		{
				QDir qdir(dir.data());
				if (!qdir.exists())
						qdir.mkpath(".");
				std::stringstream ssName;
				ssName << nrow << "_" << ncol;
				return produceImage(nrow, ncol, img, dir, ssName.str());
		}

		std::string produceImage(int nrow, int ncol, osg::Image* img, std::string dir, std::string name)
		{
				QDir qdir(dir.data());
				if (!qdir.exists())
						qdir.mkpath(".");
				std::stringstream ssSrc;
				ssSrc << dir << name << "_2.tif";
				std::stringstream ssDest;
				ssDest << dir << name << ".tif";

				osgDB::writeImageFile(*img, ssSrc.str().data());
				//printf("%s\n", ssSrc.str().data());
				GDALDataset *poSrcDS = (GDALDataset *)GDALOpen(ssSrc.str().data(), GA_ReadOnly);

				double tanAngle = tan(_altAngle);
				if (_altAngle == 0 || osg::RadiansToDegrees(_altAngle) == 90)
						tanAngle = 1;
				//GDALDataset *poDstDS = (GDALDataset *)GDALOpen(ssDest.str().data(), GA_Update);

				//_curAdfGeoTransform[0] = _bb.xMin() + ncol * _srcTileSizeX;///* top left x */
				//_curAdfGeoTransform[1] = _srcTileSizeX / img->s();///* w-e pixel resolution */
				//_curAdfGeoTransform[2] = 0;///* 0 */
				//_curAdfGeoTransform[3] = _bb.yMax() - nrow * _srcTileSizeY / tanAngle;// /* top left y */
				//_curAdfGeoTransform[4] = 0;///* 0 */
				//_curAdfGeoTransform[5] = -_srcTileSizeX / img->t();///* n-s pixel resolution (negative value) */



				char **papszOptions = NULL;

				papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");
				papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "LZW");
				GDALDataset * poDstDS = _poDriver->CreateCopy(ssDest.str().data(), poSrcDS, FALSE,
						papszOptions, GDALTermProgress, NULL);

				poDstDS->SetGeoTransform(_curAdfGeoTransform);
				//GDALClose((GDALDatasetH)poDstDS);

				/* Once we're done, close properly the dataset */
				if (poDstDS != NULL)
						GDALClose((GDALDatasetH)poDstDS);
				CSLDestroy(papszOptions);
				GDALClose((GDALDatasetH)poSrcDS);
				QFile::remove(ssSrc.str().data());
				return ssDest.str();
		}

		void produceImage(int id, osg::Image* img)
		{

				std::stringstream ssSrc;
				ssSrc << "img/" << id << "_2.tif";
				std::stringstream ssDest;
				ssDest << "img/" << id << ".tif";

				osgDB::writeImageFile(*img, ssSrc.str().data());
				GDALDataset *poSrcDS = (GDALDataset *)GDALOpen(ssSrc.str().data(), GA_ReadOnly);
				double tanAngle = tan(_altAngle);
				if (_altAngle == 0 || osg::RadiansToDegrees(_altAngle) == 90)
						tanAngle = 1;
				//GDALDataset *poDstDS = (GDALDataset *)GDALOpen(ssDest.str().data(), GA_Update);
				//double adfGeoTransform[6];
				//adfGeoTransform[0] = _bb.xMin() + _curCol * _srcTileSizeX;///* top left x */
				//adfGeoTransform[1] = _srcTileSizeX / img->s();///* w-e pixel resolution */
				//adfGeoTransform[2] = 0;///* 0 */
				//adfGeoTransform[3] = _bb.yMax() - _curRow * _srcTileSizeY / tanAngle;// /* top left y */
				//adfGeoTransform[4] = 0;///* 0 */
				//adfGeoTransform[5] = -_srcTileSizeX / img->t();///* n-s pixel resolution (negative value) */


				char **papszOptions = NULL;

				papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");
				papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "LZW");
				GDALDataset * poDstDS = _poDriver->CreateCopy(ssDest.str().data(), poSrcDS, FALSE,
						papszOptions, GDALTermProgress, NULL);

				poDstDS->SetGeoTransform(_curAdfGeoTransform);
				//GDALClose((GDALDatasetH)poDstDS);

				/* Once we're done, close properly the dataset */
				if (poDstDS != NULL)
						GDALClose((GDALDatasetH)poDstDS);
				CSLDestroy(papszOptions);
				GDALClose((GDALDatasetH)poSrcDS);
				QFile::remove(ssSrc.str().data());

		}

		bool findPixelPosition(osg::Image* img, osg::Vec3d& destpos)
		{
				//double adfGeoTransform[6];
				//double halfsize = _srcTileSizeY * 0.5;
				//double tanAngle = tan(_altAngle);
				//if (_altAngle == 0 || osg::RadiansToDegrees(_altAngle) == 90)
				//	tanAngle = 1;
				//adfGeoTransform[0] = _bb.xMin() + _curCol * _srcTileSizeY;///* top left x */
				//adfGeoTransform[1] = _srcTileSizeY / img->s();///* w-e pixel resolution */
				//adfGeoTransform[2] = 0;///* 0 */
				//adfGeoTransform[3] = _bb.yMax() - _curRow * _srcTileSizeY / tanAngle;// /* top left y */
				//adfGeoTransform[4] = 0;///* 0 */
				//adfGeoTransform[5] = -_srcTileSizeY / img->t();///* n-s pixel resolution (negative value) */
				std::vector<osg::Vec3d> pixels;
				unsigned char* data = img->data();
				for (size_t i = 0; i < img->t(); i++)
				{
						double y = _curAdfGeoTransform[3] + _curAdfGeoTransform[5] * (img->t() - i - 1) + _curAdfGeoTransform[5] * 0.5;
						for (size_t j = 0; j < img->s(); j++)
						{
								double x = _curAdfGeoTransform[0] + _curAdfGeoTransform[1] * j + _curAdfGeoTransform[1] * 0.5;
								if (*data > 0.5)
										pixels.push_back(osg::Vec3d(x, y, 0));
								data += 4;
						}
				}
				if (pixels.size() == 0)
						return false;
				destpos = osg::Vec3d(0, 0, 0);
				for (size_t i = 0; i < pixels.size(); i++)
				{
						destpos += pixels[i];
				}

				destpos = destpos / pixels.size();

				return true;
		}

};

#include "osg/Texture1D"
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
		HeightVisitor() :osg::NodeVisitor(TRAVERSE_ALL_CHILDREN) { setNodeMaskOverride(0xffffffff); }

		void Create(osg::BoundingBoxd bb, double cellsize)
		{
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
						HeightField[i] = DBL_MIN;
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
								if (height <= DBL_MIN)
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
										osg::Vec3 pos = (*vertices)[k] * matWorld;
										double col, row;
										if (!GetIndex(pos.x(), pos.y(), col, row))
												continue;
										int index = col + row * Columns;
										if (HeightField[index] == DBL_MIN || HeightField[index] < pos.z())
												HeightField[index] = pos.z();
								}
						}
						else if (dynamic_cast<osg::Vec3dArray*>(geom->getVertexArray()))
						{
								osg::Vec3dArray* vertices = (osg::Vec3dArray*)geom->getVertexArray();
								for (unsigned int k = 0; k < vertices->size(); k++)
								{
										osg::Vec3d pos = (*vertices)[k] * matWorld;
										double col, row;
										if (!GetIndex(pos.x(), pos.y(), col, row))
												continue;
										int index = col + row * Columns;
										if (HeightField[index] == DBL_MIN || HeightField[index] < pos.z())
												HeightField[index] = pos.z();
								}
						}

				}
		}


private:
};

int main(int argc, char** argv)
{
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
		int texSizeX = 4096;
		int texSizeY = 4096;


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


		osg::ref_ptr<osg::Node> node = ModelLoader::Load3DTiles(infile.data());

		node->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
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
		viewer.setUpViewInWindow(0, 0, texSizeX, texSizeY);
		//orthoCamera->attach(osg::Camera::COLOR_BUFFER0, tex.get());
		//orthoCamera->attach(osg::Camera::COLOR_BUFFER0, img.get());
		orthoCamera->addChild(sceneNode.get());

		osg::Matrix mat = osg::Matrix::rotate(osg::DegreesToRadians(azimuthAngle), osg::Vec3(0, 0, 1));
		sceneNode->setMatrix(mat);
		osg::ComputeBoundsVisitor visitor;
		sceneNode->accept(visitor);
		osg::BoundingBox bb = visitor.getBoundingBox();
		osg::Vec3 center = bb.center();
		osg::Vec2 bbSize(bb.xMax() - bb.xMin(), bb.yMax() - bb.yMin());

		osg::ComputeBoundsVisitor rawvisitor;
		node->accept(rawvisitor);
		HeightVisitor heightField;
		heightField.Create(rawvisitor.getBoundingBox(), 10);
		node->accept(heightField);
		//if (argc > 9)
		//{
		//		double xmin, ymin, xmax, ymax;
		//		xmin = atof(argv[6]); xmax = atof(argv[7]);
		//		ymin = atof(argv[8]); ymax = atof(argv[9]);
		//		std::vector<osg::Vec3d> corners;
		//		corners.push_back(osg::Vec3d(xmin, ymin, bb.zMax()));
		//		corners.push_back(osg::Vec3d(xmin, ymax, bb.zMin()));
		//		corners.push_back(osg::Vec3d(xmax, ymin, bb.zMin()));
		//		corners.push_back(osg::Vec3d(xmax, ymax, bb.zMax()));
		//		osg::BoundingBox newbb;
		//		for (size_t i = 0; i < corners.size(); i++)
		//		{
		//				osg::Vec3d newCorner = corners[i] * mat;
		//				newbb.expandBy(newCorner);
		//		}

		//		bb = newbb;
		//}

#pragma region Voxel Shaders
		//char vertexShaderSource[] =
	//		"#extension GL_EXT_gpu_shader4 : enable\n"
	//		"uniform mat4 osg_ViewMatrixInverse;\n"
	//		"varying vec4 pos;\n"
	//		"void main(void)\n"
	//		"{\n"
	//		"gl_TexCoord[0] = gl_MultiTexCoord0;\n"
	//		"pos = gl_Vertex;\n"
	//		"gl_Position = gl_ModelViewProjectionMatrix * pos;\n"
	//		"}\n";

	//char fragmentShaderSource[] =
	//		"#extension GL_EXT_gpu_shader4 : enable\n"
	//		"vec4 unpackColor(float f)\n"
	//		"{\n"
	//		"		vec3 color;\n"

	//		"		color.r = floor(f / 256.0 / 256.0);\n"
	//		"		color.g = floor((f - color.r * 256.0 * 256.0) / 256.0);\n"
	//		"		color.b = floor(f - color.r * 256.0 * 256.0 - color.g * 256.0);\n"

	//		// now we have a vec3 with the 3 components in range [0..256]. Let's normalize it!
	//		"		return vec4(color / 256.0, 0);\n"
	//		"	}\n"

	//		"uniform sampler2D tex;\n"
	//		"uniform sampler1D colorRamp; \n"
	//		"uniform vec3 voxelSize;\n"
	//		"uniform vec3 voxelDims;\n"
	//		"uniform vec3 minbound;\n"
	//		"uniform vec3 maxbound;\n"
	//		"varying vec4 pos;\n"
	//		"void main(void) \n"
	//		"{\n"
	//		"  vec3 voxelIndexXYZ = floor((pos.xyz - minbound) / voxelSize);\n"
	//		"  int voxelIndex = int(voxelDims.x * voxelDims.y * voxelIndexXYZ.z +  voxelDims.x * voxelIndexXYZ.y + voxelIndexXYZ.x);\n"
	//		"  vec3 color = texelFetch(colorRamp, voxelIndex, 0).xyz; \n"
	//		"  //color = normalize(voxelIndex);\n"
	//		"  gl_FragColor = texture2D(tex, gl_TexCoord[0].xy) * vec4(color, 1);\n"
	//		"}\n";


	//osg::Vec3 minbound = osg::Vec3(bb.xMin(), bb.yMin(), bb.zMin());
	//osg::Vec3 maxbound = osg::Vec3(bb.xMax(), bb.yMax(), bb.zMax());
	//osg::Vec3 voxelSize = osg::Vec3((bb.xMax() - bb.xMin()) / 20, (bb.yMax() - bb.yMin()) / 20, (bb.zMax() - bb.zMin()) / 50);
	//osg::Vec3 voxelDims = (maxbound - minbound);
	//voxelDims = osg::Vec3(ceil(voxelDims.x() / voxelSize.x()), ceil(voxelDims.y() / voxelSize.y()), ceil(voxelDims.z() / voxelSize.z()));

	//sceneNode->getOrCreateStateSet()->addUniform(new osg::Uniform("tex", 0));
	//sceneNode->getOrCreateStateSet()->addUniform(new osg::Uniform("colorRamp", 1));
	//sceneNode->getOrCreateStateSet()->addUniform(new osg::Uniform("minbound", minbound));
	//sceneNode->getOrCreateStateSet()->addUniform(new osg::Uniform("maxbound", maxbound));
	//sceneNode->getOrCreateStateSet()->addUniform(new osg::Uniform("voxelSize", voxelSize));
	//sceneNode->getOrCreateStateSet()->addUniform(new osg::Uniform("voxelDims", voxelDims));
	//osg::ref_ptr <osg::Texture1D> colorRampTex = createRandomColors(voxelDims.x()*voxelDims.y()*voxelDims.z());
	//colorRampTex->setResizeNonPowerOfTwoHint(false);
	//sceneNode->getOrCreateStateSet()->setTextureAttribute(1, colorRampTex.get());
	//printfVec3(voxelSize, voxelDims);
	//osg::ref_ptr<osg::Program> program = new osg::Program;
	//program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));
	//program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexShaderSource));
	/*osg::ref_ptr<osg::Image> maskImg = osgEarth::URI(baseDir+"mask.png").getImage();
	if(!maskImg || !maskImg.valid())
	{
		maskImg = new osg::Image;
		maskImg->allocateImage(1,1,1,GL_RGBA,GL_UNSIGNED_BYTE);
		unsigned char* data = maskImg->data();
		data[0]=255;data[1]=255;data[2]=255;data[3]=255;

	}*/

	//sceneNode->getOrCreateStateSet()->setAttribute(program.get(), osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
#pragma endregion

#pragma region 
		char vertexShaderSource[] =
			"#extension GL_EXT_gpu_shader4 : enable\n"
			"uniform mat4 osg_ViewMatrixInverse;\n"
			"varying vec4 pos;\n"
			"void main(void)\n"
			"{\n"
			"gl_TexCoord[0] = gl_MultiTexCoord0;\n"
			"pos = gl_Vertex;\n"
			"gl_Position = gl_ModelViewProjectionMatrix * pos;\n"
			"}\n";

	char fragmentShaderSource[] =
			"#extension GL_EXT_gpu_shader4 : enable\n"
			"uniform sampler2D tex;\n"
			"uniform vec3 minbound;\n"
			"uniform vec3 maxbound;\n"
			"varying vec4 pos;\n"
			"void main(void) \n"
			"{\n"
			"  vec3 color = texture2D(tex, gl_TexCoord[0].xy); \n"
			"  if(pos.x > minbound.x && pos.x < maxbound.x && pos.y > minbound.y && pos.y < maxbound.y)\n"
			"  {\n"
			"     color = color + vec3(0.3,0,0);\n"
			"  }\n"
			"  gl_FragColor = vec4(color, 1);\n"
			"}\n";


	osg::Vec3 minbound = osg::Vec3(bb.xMin(), bb.yMin(), bb.zMin());
	osg::Vec3 maxbound = osg::Vec3(bb.xMax(), bb.yMax(), bb.zMax());
	sceneNode->getOrCreateStateSet()->addUniform(new osg::Uniform("tex", 0));
	sceneNode->getOrCreateStateSet()->addUniform(new osg::Uniform("minbound", minbound));
	sceneNode->getOrCreateStateSet()->addUniform(new osg::Uniform("maxbound", maxbound));
	osg::ref_ptr<osg::Program> program = new osg::Program;
	program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));
	program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexShaderSource));
	sceneNode->getOrCreateStateSet()->setAttribute(program.get(), osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
#pragma endregion


		
		//azimuthAngle = 0;
		//while (azimuthAngle < 360)
		//{
		//		mat = osg::Matrix::rotate(osg::DegreesToRadians(azimuthAngle), osg::Vec3(0, 0, 1));
		//		sceneNode->setMatrix(mat);
		//		center = osg::Vec3(-3945.2300, -7068.9270, 0);
		//		bbSize = osg::Vec2(1000, 1000);
		//		double minz, maxz;
		//		BBWrapper localbb = BBWrapper(center.x() - bbSize.x() * 0.5, center.y() - bbSize.y() * 0.5, bb.zMin(),
		//				center.x() + bbSize.x() * 0.5, center.y() + bbSize.y() * 0.5, bb.zMax());
		//		heightField.GetHeightRange(localbb, minz, maxz);
		//		localbb = BBWrapper(center.x() - bbSize.x() * 0.5, center.y() - bbSize.y() * 0.5, minz,
		//				center.x() + bbSize.x() * 0.5, center.y() + bbSize.y() * 0.5, maxz);

		//		double scale = 1.0;
		//		minbound = osg::Vec3(center.x() - localbb.xhalfsize()*scale, center.y() - localbb.yhalfsize()*scale, center.z() - localbb.zhalfsize());
		//		maxbound = osg::Vec3(center.x() + localbb.xhalfsize()*scale, center.y() + localbb.yhalfsize()*scale, center.z() + localbb.zhalfsize());
		//		sceneNode->getOrCreateStateSet()->getUniform("minbound")->set(minbound);
		//		sceneNode->getOrCreateStateSet()->getUniform("maxbound")->set(maxbound);
		//		osg::BoundingBoxd transformedBB;
		//		transformedBB.init();
		//		for (size_t i = 0; i < 8; i++)
		//		{
		//				transformedBB.expandBy(localbb.corner(i) * mat);
		//		}
		//		localbb = BBWrapper(transformedBB);


		//		//visitor.reset();
		//		//sceneNode->accept(visitor);
		//		//bb = visitor.getBoundingBox();
		//		//osg::BoundingBox newBB;
		//		//newBB.init();
		//		//newBB.expandBy(minbound * mat);
		//		//newBB.expandBy(maxbound * mat);
		//		//bbSize = osg::Vec2(newBB.xMax() - newBB.xMin(), newBB.yMax() - newBB.yMin());
		//		//osg::Vec3 transformedCenter = center * mat;
		//		//bb = osg::BoundingBox(transformedCenter.x() - bbSize.x() * 0.5, transformedCenter.y() - bbSize.y() * 0.5, bb.zMin(),
		//		//		transformedCenter.x() + bbSize.x() * 0.5, transformedCenter.y() + bbSize.y() * 0.5, bb.zMax());

		//		orthoCamera->reset(localbb, altAngle, azimuthAngle);
		//		viewer.frame();
		//		while (pager->getRequestsInProgress())
		//		{
		//				viewer.frame();
		//		}
		//		std::stringstream ssoutname;
		//		ssoutname << outdir << altAngle << "_" << azimuthAngle << ".png";
		//		printf("%s\n", ssoutname.str().data());
		//		osgDB::writeImageFile(*img, ssoutname.str().data());
		//		azimuthAngle += 45;
		//}



		struct ImageOutput
		{
				osg::ref_ptr<osg::Texture2D> Tex;
				osg::ref_ptr<osg::Image> Image;
				int Azimuth;
				int Alt;

				void Write(std::string outdir)
				{
						std::stringstream ssoutname;
						ssoutname << outdir << Alt << "_" << Azimuth << ".png";
						printf("%s\n", ssoutname.str().data());
						osgDB::writeImageFile(*Image, ssoutname.str().data());
				}

				void Create(int texSizeX, int texSizeY, int azimuth, int alt)
				{
						Tex = new osg::Texture2D;
						Tex->setTextureSize(texSizeX, texSizeY);
						Tex->setResizeNonPowerOfTwoHint(false);
						Tex->setFilter(osg::Texture2D::MIN_FILTER, osg::Texture2D::LINEAR);
						Tex->setFilter(osg::Texture2D::MAG_FILTER, osg::Texture2D::LINEAR);
						Tex->setWrap(osg::Texture2D::WRAP_S, osg::Texture2D::REPEAT);
						Tex->setWrap(osg::Texture2D::WRAP_T, osg::Texture2D::REPEAT);
						//rtTexture->setDataVariance(osg::Object::DYNAMIC);
						Tex->setInternalFormat(GL_RGBA);
						Tex->setSourceFormat(GL_RGBA);
						Tex->setSourceType(GL_UNSIGNED_BYTE);

						Image = new osg::Image;
						Image->allocateImage(texSizeX, texSizeY, 1, GL_RGBA, GL_UNSIGNED_BYTE);

						Tex->setImage(Image.get());
						Azimuth = azimuth;
						Alt = alt;
				}
		};

		std::vector<ImageOutput*> imageArr;
		//osg::ref_ptr<osg::Texture2D> tex = new osg::Texture2D;

		//tex->setTextureSize(texSizeX, texSizeY);
		//tex->setResizeNonPowerOfTwoHint(false);
		//tex->setFilter(osg::Texture2D::MIN_FILTER, osg::Texture2D::LINEAR);
		//tex->setFilter(osg::Texture2D::MAG_FILTER, osg::Texture2D::LINEAR);
		//tex->setWrap(osg::Texture2D::WRAP_S, osg::Texture2D::REPEAT);
		//tex->setWrap(osg::Texture2D::WRAP_T, osg::Texture2D::REPEAT);
		////rtTexture->setDataVariance(osg::Object::DYNAMIC);
		//tex->setInternalFormat(GL_RGBA);
		//tex->setSourceFormat(GL_RGBA);
		//tex->setSourceType(GL_UNSIGNED_BYTE);
		ImageOutput* imageOutput = new ImageOutput;
		imageOutput->Create(texSizeX, texSizeY, azimuthAngle, altAngle);
		//osg::ref_ptr<osg::Image> img = new osg::Image;
		//img->allocateImage(texSize, texSize,1,GL_ALPHA,GL_FLOAT);
		//img->allocateImage(texSizeX, texSizeY, 1, GL_RGBA, GL_UNSIGNED_BYTE);
		//img->allocateImage(texSize, texSize,1,GL_RGBA,GL_FLOAT); 

		//tex->setImage(imageOutput->Image.get());
		orthoCamera->attach(osg::Camera::COLOR_BUFFER0, imageOutput->Tex.get());
		imageArr.push_back(imageOutput);
		orthoCamera->attach(osg::Camera::COLOR_BUFFER0, imageOutput->Image.get());
		osg::ref_ptr<ScreenOverlay> pOverlay = new ScreenOverlay(&viewer);
		pOverlay->setTextureLayer(imageOutput->Tex.get());
		osg::ref_ptr<osg::Group> root = new osg::Group;
		//root->addChild(node.get());
		root->addChild(orthoCamera.get());
		root->addChild(pOverlay.get());

		viewer.setSceneData(root.get());
		//viewer.setUpViewInWindow(0, 0, 1024, 768);
		viewer.realize();
		viewer.setCameraManipulator(new osgGA::TrackballManipulator);
		//viewer.getCamera()->getOrCreateStateSet()->setMode(GL_BLEND, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
		//orthoCamera->getOrCreateStateSet()->setAttribute(clamp, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
		//STPager* databasePager = new STPager;
		//databasePager->cancel();
		//viewer.getScene()->setDatabasePager(databasePager);
		osgDB::DatabasePager* pager = viewer.getScene()->getDatabasePager();
		//return 0;
		GDALAllRegister();

		azimuthAngle = 0;
		altAngle = 90;
		mat = osg::Matrix::rotate(osg::DegreesToRadians(azimuthAngle), osg::Vec3(0, 0, 1));
		sceneNode->setMatrix(mat);
		center = osg::Vec3(-3945.2300, -7068.9270, 0);
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
		sceneNode->getOrCreateStateSet()->getUniform("minbound")->set(minbound);
		sceneNode->getOrCreateStateSet()->getUniform("maxbound")->set(maxbound);
		osg::BoundingBoxd transformedBB;
		transformedBB.init();
		for (size_t i = 0; i < 8; i++)
		{
				transformedBB.expandBy(localbb.corner(i) * mat);
		}
		localbb = BBWrapper(transformedBB);
		orthoCamera->reset(localbb, altAngle, azimuthAngle);
		viewer.frame();
		while (pager->getRequestsInProgress())
		{
				viewer.frame();
		}
	
		GeometryVisitor geomVisitor;
		geomVisitor.initialize("E:/Code/TestMesh/");
		geomVisitor.m_BB = localbb;
		//node->accept(geomVisitor);
		geomVisitor.traverseGroup(node);
		geomVisitor.save("E:/Code/TestMesh/TestMesh.osg");

		altAngle = 90;
		bool isFirst = true;

		altAngle = 90;
		while (altAngle > 5 && !viewer.done())
		{
				azimuthAngle = 0;
				while (azimuthAngle < 360 && !viewer.done())
				{
						mat = osg::Matrix::rotate(osg::DegreesToRadians(azimuthAngle), osg::Vec3(0, 0, 1));
						sceneNode->setMatrix(mat);
						center = osg::Vec3(-3945.2300, -7068.9270, 0);
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
						sceneNode->getOrCreateStateSet()->getUniform("minbound")->set(minbound);
						sceneNode->getOrCreateStateSet()->getUniform("maxbound")->set(maxbound);
						osg::BoundingBoxd transformedBB;
						transformedBB.init();
						for (size_t i = 0; i < 8; i++)
						{
								transformedBB.expandBy(localbb.corner(i) * mat);
						}
						localbb = BBWrapper(transformedBB);
						orthoCamera->reset(localbb, altAngle, azimuthAngle);

						viewer.frame();
						//if (isFirst)
						//{
						//		while (pager->getRequestsInProgress())
						//		{
						//				viewer.frame();
						//		}
						//		isFirst = false;
						//}
						//viewer.frame();
						//orthoCamera->dirtyAttachmentMap();
						int count = 0;
						/*while (pager->getRequestsInProgress())
						{
								viewer.frame();
								count++;
								if (count > 5)
										break;
						}*/
						//std::stringstream ssoutname;
						//ssoutname << outdir << altAngle << "_" << azimuthAngle << ".png";
						//printf("%s\n", ssoutname.str().data());
						//osgDB::writeImageFile(*img, ssoutname.str().data());
						azimuthAngle += 12;
				}
					altAngle -= 3;
		}

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

