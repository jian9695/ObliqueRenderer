
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
using namespace osgDB;
using namespace OpenThreads;

class ShapeFile
{
public:

	ShapeFile()
	{
		poDS = NULL;
	}
	ShapeFile(std::string filename, int update = 0) {
		g_mFileName = filename;
		poDS = NULL;
		poDS = (GDALDataset*)GDALOpenEx(filename.data(), GDAL_OF_VECTOR | update, NULL, NULL, NULL);
		poLayer = poDS->GetLayer(0);
	}

	void close()
	{
		if (poDS)
			GDALClose(poDS);
		poDS = NULL;
	}
	void create(const char* filename, OGRSpatialReference* spatialRef = NULL, OGRFeatureDefn *poFDefn = NULL, OGRwkbGeometryType geotype = wkbPoint)
	{


		g_mFileName = filename;

		const char *pszDriverName = "ESRI Shapefile";
		GDALDriver *poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(
			pszDriverName);
		//poDriver->Delete(filename.data());

		if (poDS)
			close();
		//const char *pszDriverName = "ESRI Shapefile";
		//OGRSFDriver *poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(
		//	pszDriverName);

		if (QFileInfo(filename).exists())
		{
			poDriver->Delete(filename);
			//poDS = OGRSFDriverRegistrar::Open(filename, TRUE);
			//poLayer = poDS->GetLayer(0);
		}
		poDS = poDriver->Create(filename, 0, 0, 0, GDT_Unknown, NULL);
		poLayer = poDS->CreateLayer(QFileInfo(filename).baseName().toLocal8Bit().data(), spatialRef, geotype, NULL);

		if (poFDefn)
		{
			for (int iField = 0; iField < poFDefn->GetFieldCount(); iField++)
			{
				OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);
				poLayer->CreateField(poFieldDefn);
			}
		}
		//}
	}

	~ShapeFile() {
		close();
	}
public:
	OGRLayer       *poLayer;
	GDALDataset    *poDS;
private:

	std::string g_mFileName;
};
class OptionsReadFileCallback : public osgDB::ReadFileCallback
{
	osgDB::ReaderWriter::ReadResult
		readNode(const std::string& name, const osgDB::Options* options)
	{

		if(name.length() > 5 && name.substr(0,4) == "http")
		{
		std::string str = name;
		for ( std::string::iterator it=str.begin(); it!=str.end(); ++it)
		{
		if( *it == '+')
		{
		*it = '-';
		}
		}
		return osgDB::ReadFileCallback::readNode(str, options);
		}
		osgDB::ReaderWriter::ReadResult result = osgDB::ReadFileCallback::readNode(name, options);
	/*	osgUtil::SmoothingVisitor visitor;
		result.getNode()->accept(visitor);*/
		return result;
	}
	osgDB::ReaderWriter::ReadResult readImage(const std::string& filename, const osgDB::Options* options)
	{
		 osgDB::ReaderWriter::ReadResult result = osgDB::ReadFileCallback::readImage(filename, options);
		 return result;
	}
	virtual osgDB::ReaderWriter::ReadResult readObject(const std::string& filename, const osgDB::Options* options)
	{
		osgDB::ReaderWriter::ReadResult result = osgDB::ReadFileCallback::readObject(filename, options);
		return result;
	}
   
};
struct DatabasePager::DatabasePagerCompileCompletedCallback : public osgUtil::IncrementalCompileOperation::CompileCompletedCallback
{
	DatabasePagerCompileCompletedCallback(osgDB::DatabasePager* pager, osgDB::DatabasePager::DatabaseRequest* databaseRequest):
_pager(pager),
	_databaseRequest(databaseRequest) {}

virtual bool compileCompleted(osgUtil::IncrementalCompileOperation::CompileSet* /*compileSet*/)
{
	_pager->compileCompleted(_databaseRequest.get());
	return true;
}

osgDB::DatabasePager*                               _pager;
osg::ref_ptr<osgDB::DatabasePager::DatabaseRequest> _databaseRequest;
};

class DatabasePager::FindCompileableGLObjectsVisitor : public osgUtil::StateToCompile
{
public:
	FindCompileableGLObjectsVisitor(const DatabasePager* pager, osg::Object* markerObject):
	  osgUtil::StateToCompile(osgUtil::GLObjectsVisitor::COMPILE_DISPLAY_LISTS|osgUtil::GLObjectsVisitor::COMPILE_STATE_ATTRIBUTES, markerObject),
		  _pager(pager),
		  _changeAutoUnRef(false), _valueAutoUnRef(false),
		  _changeAnisotropy(false), _valueAnisotropy(1.0)
	  {
		  _assignPBOToImages = _pager->_assignPBOToImages;

		  _changeAutoUnRef = _pager->_changeAutoUnRef;
		  _valueAutoUnRef = _pager->_valueAutoUnRef;
		  _changeAnisotropy = _pager->_changeAnisotropy;
		  _valueAnisotropy = _pager->_valueAnisotropy;

		  switch(_pager->_drawablePolicy)
		  {
		  case DatabasePager::DO_NOT_MODIFY_DRAWABLE_SETTINGS:
			  // do nothing, leave settings as they came in from loaded database.
			  // OSG_NOTICE<<"DO_NOT_MODIFY_DRAWABLE_SETTINGS"<<std::endl;
			  break;
		  case DatabasePager::USE_DISPLAY_LISTS:
			  _mode = _mode | osgUtil::GLObjectsVisitor::SWITCH_ON_DISPLAY_LISTS;
			  _mode = _mode | osgUtil::GLObjectsVisitor::SWITCH_OFF_VERTEX_BUFFER_OBJECTS;
			  _mode = _mode & ~osgUtil::GLObjectsVisitor::SWITCH_ON_VERTEX_BUFFER_OBJECTS;
			  break;
		  case DatabasePager::USE_VERTEX_BUFFER_OBJECTS:
			  _mode = _mode | osgUtil::GLObjectsVisitor::SWITCH_ON_VERTEX_BUFFER_OBJECTS;
			  break;
		  case DatabasePager::USE_VERTEX_ARRAYS:
			  _mode = _mode & ~osgUtil::GLObjectsVisitor::SWITCH_ON_DISPLAY_LISTS;
			  _mode = _mode & ~osgUtil::GLObjectsVisitor::SWITCH_ON_VERTEX_BUFFER_OBJECTS;
			  _mode = _mode | osgUtil::GLObjectsVisitor::SWITCH_OFF_DISPLAY_LISTS;
			  _mode = _mode | osgUtil::GLObjectsVisitor::SWITCH_OFF_VERTEX_BUFFER_OBJECTS;
			  break;
		  }

		  if (osgDB::Registry::instance()->getBuildKdTreesHint()==osgDB::Options::BUILD_KDTREES &&
			  osgDB::Registry::instance()->getKdTreeBuilder())
		  {
			  _kdTreeBuilder = osgDB::Registry::instance()->getKdTreeBuilder()->clone();
		  }
	  }

	  META_NodeVisitor("osgDB","FindCompileableGLObjectsVisitor")

		  bool requiresCompilation() const { return !empty(); }

	  virtual void apply(osg::Geode& geode)
	  {
		  StateToCompile::apply(geode);

		  if (_kdTreeBuilder.valid())
		  {
			  geode.accept(*_kdTreeBuilder);
		  }
	  }

	  void apply(osg::Texture& texture)
	  {
		  StateToCompile::apply(texture);

		  if (_changeAutoUnRef)
		  {
			  texture.setUnRefImageDataAfterApply(_valueAutoUnRef);
		  }

		  if ((_changeAnisotropy && texture.getMaxAnisotropy() != _valueAnisotropy))
		  {
			  texture.setMaxAnisotropy(_valueAnisotropy);
		  }
	  }

	  const DatabasePager*                    _pager;
	  bool                                    _changeAutoUnRef;
	  bool                                    _valueAutoUnRef;
	  bool                                    _changeAnisotropy;
	  float                                   _valueAnisotropy;
	  osg::ref_ptr<osg::KdTreeBuilder>        _kdTreeBuilder;

protected:

	FindCompileableGLObjectsVisitor& operator = (const FindCompileableGLObjectsVisitor&) { return *this; }
};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  SortFileRequestFunctor
//
struct DatabasePager::SortFileRequestFunctor
{
	bool operator() (const osg::ref_ptr<DatabasePager::DatabaseRequest>& lhs,const osg::ref_ptr<DatabasePager::DatabaseRequest>& rhs) const
	{
		if (lhs->_timestampLastRequest>rhs->_timestampLastRequest) return true;
		else if (lhs->_timestampLastRequest<rhs->_timestampLastRequest) return false;
		else return (lhs->_priorityLastRequest>rhs->_priorityLastRequest);
	}
};

class STPager : public osgDB::DatabasePager
{

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//  FindCompileableGLObjectsVisitor
	//
	


public:
	STPager()
		:osgDB::DatabasePager()
	{

	}

	void frame()
	{


		//bool firstTime = true;

		osg::ref_ptr<DatabasePager::ReadQueue> read_queue;
		osg::ref_ptr<DatabasePager::ReadQueue> out_queue;


		read_queue = _fileRequestQueue;



		//read_queue->block();


		//
		// delete any children if required.
		//
		if (_deleteRemovedSubgraphsInDatabaseThread/* && !(read_queue->_childrenToDeleteList.empty())*/)
		{
			ObjectList deleteList;
			{
				// Don't hold lock during destruction of deleteList
				OpenThreads::ScopedLock<OpenThreads::Mutex> lock(read_queue->_requestMutex);
				if (!read_queue->_childrenToDeleteList.empty())
				{
					deleteList.swap(read_queue->_childrenToDeleteList);
					read_queue->updateBlock();
				}
			}
		}

		//
		// load any subgraphs that are required.
		//
		osg::ref_ptr<DatabaseRequest> databaseRequest;
		read_queue->takeFirst(databaseRequest);

		bool readFromFileCache = false;

		osg::ref_ptr<FileCache> fileCache = osgDB::Registry::instance()->getFileCache();
		osg::ref_ptr<FileLocationCallback> fileLocationCallback = osgDB::Registry::instance()->getFileLocationCallback();
		osg::ref_ptr<Options> dr_loadOptions;
		std::string fileName;
		int frameNumberLastRequest = 0;
		if (databaseRequest.valid())
		{
			{
				OpenThreads::ScopedLock<OpenThreads::Mutex> drLock(_dr_mutex);
				dr_loadOptions = databaseRequest->_loadOptions;
				fileName = databaseRequest->_fileName;
				frameNumberLastRequest = databaseRequest->_frameNumberLastRequest;
			}
			if (dr_loadOptions.valid())
			{
				if (dr_loadOptions->getFileCache()) fileCache = dr_loadOptions->getFileCache();
				if (dr_loadOptions->getFileLocationCallback()) fileLocationCallback = dr_loadOptions->getFileLocationCallback();

				dr_loadOptions = dr_loadOptions->cloneOptions();
			}
			else
			{
				dr_loadOptions = new osgDB::Options;
			}

			dr_loadOptions->setTerrain(databaseRequest->_terrain);

			// disable the FileCache if the fileLocationCallback tells us that it isn't required for this request.
			if (fileLocationCallback.valid() && !fileLocationCallback->useFileCache()) fileCache = 0;


			// check if databaseRequest is still relevant
			if ((_frameNumber-frameNumberLastRequest)<=1)
			{


				// do nothing as this thread can handle the load
				if (fileCache.valid() && fileCache->isFileAppropriateForFileCache(fileName))
				{
					if (fileCache->existsInCache(fileName))
					{
						readFromFileCache = true;
					}
				}


			}

		}
		if (databaseRequest.valid())
		{

			// load the data, note safe to write to the databaseRequest since once
			// it is created this thread is the only one to write to the _loadedModel pointer.
			//OSG_NOTICE<<"In DatabasePager thread readNodeFile("<<databaseRequest->_fileName<<")"<<std::endl;
			//osg::Timer_t before = osg::Timer::instance()->tick();


			// assume that readNode is thread safe...
			ReaderWriter::ReadResult rr = readFromFileCache ?
				fileCache->readNode(fileName, dr_loadOptions.get(), false) :
			Registry::instance()->readNode(fileName, dr_loadOptions.get(), false);

			osg::ref_ptr<osg::Node> loadedModel;
			if (rr.validNode()) loadedModel = rr.getNode();
			if (rr.error()) OSG_WARN<<"Error in reading file "<<fileName<<" : "<<rr.message() << std::endl;
			if (rr.notEnoughMemory()) OSG_INFO<<"Not enought memory to load file "<<fileName << std::endl;

			if (loadedModel.valid() &&
				fileCache.valid() &&
				fileCache->isFileAppropriateForFileCache(fileName) &&
				!readFromFileCache)
			{
				fileCache->writeNode(*(loadedModel), fileName, dr_loadOptions.get());
			}

			{
				OpenThreads::ScopedLock<OpenThreads::Mutex> drLock(_dr_mutex);
				if ((_frameNumber-databaseRequest->_frameNumberLastRequest)>1)
				{

					loadedModel = 0;
				}
			}

			//OSG_NOTICE<<"     node read in "<<osg::Timer::instance()->delta_m(before,osg::Timer::instance()->tick())<<" ms"<<std::endl;

			if (loadedModel.valid())
			{
				loadedModel->getBound();

				// find all the compileable rendering objects
				DatabasePager::FindCompileableGLObjectsVisitor stateToCompile(this, this->getMarkerObject());
				loadedModel->accept(stateToCompile);

				bool loadedObjectsNeedToBeCompiled = _doPreCompile &&
					_incrementalCompileOperation.valid() &&
					_incrementalCompileOperation->requiresCompile(stateToCompile);

				// move the databaseRequest from the front of the fileRequest to the end of
				// dataToCompile or dataToMerge lists.
				osg::ref_ptr<osgUtil::IncrementalCompileOperation::CompileSet> compileSet = 0;
				if (loadedObjectsNeedToBeCompiled)
				{
					// OSG_NOTICE<<"Using IncrementalCompileOperation"<<std::endl;

					compileSet = new osgUtil::IncrementalCompileOperation::CompileSet(loadedModel.get());
					compileSet->buildCompileMap(_incrementalCompileOperation->getContextSet(), stateToCompile);
					compileSet->_compileCompletedCallback = new DatabasePagerCompileCompletedCallback(this, databaseRequest.get());
					_incrementalCompileOperation->add(compileSet.get(), false);
				}
				{
					OpenThreads::ScopedLock<OpenThreads::Mutex> drLock(_dr_mutex);
					databaseRequest->_loadedModel = loadedModel;
					databaseRequest->_compileSet = compileSet;
				}
				// Dereference the databaseRequest while the queue is
				// locked. This prevents the request from being
				// deleted at an unpredictable time within
				// addLoadedDataToSceneGraph.
				if (loadedObjectsNeedToBeCompiled)
				{
					OpenThreads::ScopedLock<OpenThreads::Mutex> listLock(
						_dataToCompileList->_requestMutex);
					_dataToCompileList->addNoLock(databaseRequest.get());
					databaseRequest = 0;
				}
				else
				{
					OpenThreads::ScopedLock<OpenThreads::Mutex> listLock(
						_dataToMergeList->_requestMutex);
					_dataToMergeList->addNoLock(databaseRequest.get());
					databaseRequest = 0;
				}

			}

			// _pager->_dataToCompileList->pruneOldRequestsAndCheckIfEmpty();
		}

	
	}
	//virtual void setUpThreads(unsigned int totalNumThreads=2, unsigned int numHttpThreads=1)
	//{
	//	return;
	//}
	virtual int cancel()
	{
		int result = 0;

		for(DatabaseThreadList::iterator dt_itr = _databaseThreads.begin();
			dt_itr != _databaseThreads.end();
			++dt_itr)
		{
			(*dt_itr)->setDone(true);
		}

		return result;
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
	ObliqueCamera(const osg::BoundingBox& bb, double altAngle = 0, double azimuthAngle = 0)
	{
		_altAngle = osg::DegreesToRadians(altAngle);

		osg::Vec3d bbmin(bb.xMin(), bb.yMin(), bb.zMin());
		osg::Vec3d bbmax(bb.xMax(), bb.yMax(), bb.zMin());
		osg::Vec3d bbsize = bbmax - bbmin;

		float maxsize = bbsize.x() > bbsize.y() ? bbsize.x() : bbsize.y();
		_viewDistance = maxsize * 10;

		osg::Vec3d translation(0, -_viewDistance * cos(_altAngle), _viewDistance * sin(_altAngle));
		if (_altAngle == 0 || osg::RadiansToDegrees(_altAngle) == 90)
		{
			translation = osg::Vec3d(0, 0, _viewDistance);
		}
		osg::Vec3d bbminNear = bbmin + translation;
		osg::Vec3d bbmaxNear = bbmax + translation;
		osg::Vec3d bbminFar = bbmin - translation;
		osg::Vec3d bbmaxFar = bbmax - translation;

		_znear = _viewDistance - bb.radius() * 2;
		_zfar = _viewDistance + bb.radius() * 2;
		float  right = (bb.xMax() - bb.xMin())*0.5;

		this->setReferenceFrame(osg::Transform::ABSOLUTE_RF);

		//_upDirection = osg::Vec3d(0.0, 1.0, 0.0);
		//_upDirection.normalize();
		osg::Vec3d vforward(0, -cos(_altAngle), sin(_altAngle));
		osg::Vec3d vright(1, 0, 0);
		_upDirection = vforward ^ vright;
		//_upDirection = vright ^ vforward;

		_viewDirection = -translation;
		_viewDirection.normalize();
		_center = (bbmin + bbmax) * 0.5;

		osg::Vec3d eyePoint = _center - _viewDirection * _viewDistance;
		//this->setViewMatrixAsLookAt(eyePoint, _center, _upDirection);

		_bb = bb;

		GDALAllRegister();
		const char *pszFormat = "GTiff";

		char **papszMetadata;

		_poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
		char **papszOptions = NULL;
		//std::stringstream ssSrc;
		//ssSrc << nrow << "_" << ncol << ".tif";
	}
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
		double zspan = _bb.zMax() - _bb.zMin();
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

		_bb = bb;
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
		osg::Vec3d eyePoint = center - _viewDirection *_viewDistance;
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
		std::stringstream ssSrc;
		ssSrc << dir << nrow << "_" << ncol << "_2.tif";
		std::stringstream ssDest;
		ssDest << dir << nrow << "_" << ncol << ".tif";

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

//int main(int argc, char** argv)
//{
//	osg::DisplaySettings::instance()->setNumOfDatabaseThreadsHint(1);
//	osg::DisplaySettings::instance()->setNumOfHttpDatabaseThreadsHint(0);
//
//
//	osgViewer::Viewer viewer;
//	//viewer.getScene()->setDatabasePager()
//	// add the state manipulator
//	viewer.addEventHandler( new osgGA::StateSetManipulator(viewer.getCamera()->getOrCreateStateSet()) );
//
//	// add the thread model handler
//	viewer.addEventHandler(new osgViewer::ThreadingHandler);
//
//	// add the window size toggle handler
//	viewer.addEventHandler(new osgViewer::WindowSizeHandler);
//
//	// add the stats handler
//	viewer.addEventHandler(new osgViewer::StatsHandler);
//
//
//	// add the record camera path handler
//	//viewer.addEventHandler(new osgViewer::RecordCameraPathHandler);
//
//	// add the LOD Scale handler
//	//viewer.addEventHandler(new osgViewer::LODScaleHandler);
//
//	// add the screen capture handler
//	//viewer.addEventHandler(new osgViewer::ScreenCaptureHandler);
//	//osg::Node* node = osgDB::readNodeFile("F:/Oblique_Photogrammetry/weihai/Data/all.desc");
//	osg::Node* node = osgDB::readNodeFile("E:/efficient-sparse-voxel-octrees-1.4/scenes/kechuang.ive");
//	//osg::Node* node = osgDB::readNodeFile("F:/Oblique_Photogrammetry/xinji/Data/all.desc");
//	//osg::Node* node = osgDB::readNodeFile("F:/Oblique_Photogrammetry/weihai/Data/all.desc");
//	//osg::Node* node = osgDB::readNodeFile("F:/Oblique_Photogrammetry/hongpeng/Data/all.desc");
//	//osg::Node* node = osgDB::readNodeFile("F:/Oblique_Photogrammetry/suzhou_gaoxin/Data/all.desc");
//	/*osg::Geode* node = new osg::Geode();
//	node->addDrawable(new osg::ShapeDrawable(new osg::Box(osg::Vec3(0.0f, 0.0f, 0.0f), 2 * 1)));*/
//	//node->addDrawable(new osg::ShapeDrawable(new osg::Sphere(osg::Vec3(0.0f, 0.0f, 0.0f), 2)));
//	std::stringstream ss_outdir;
//	std::string outdir = "F:/Oblique_Photogrammetry/weihai/DOM025/";
//	//transformGCPs("../jiashan/georeference.txt",outdir + "georeference.txt",obliqueMat);
//	//osg::Node* sceneNode = 
//	//osg::Node* node = osgDB::readNodeFile("D:/Projects/Oblique_Photogrammetry/xinji/Data/all.desc");
//	node->getOrCreateStateSet()->setMode(GL_LIGHTING,osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
//	//node->getOrCreateStateSet()->setMode( GL_BLEND, osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);
//	osg::MatrixTransform* sceneNode = new osg::MatrixTransform;
//	sceneNode->addChild(node);
//	//sceneNode->setMatrix(osg::Matrix::scale(1, 5, 1));
//	osg::ComputeBoundsVisitor visitor;
//	sceneNode->accept(visitor);
//	osg::BoundingBox bb = visitor.getBoundingBox();
//	osg::BoundingBox tileBB;
//
//	double resolution = 0.5;
//	float bbW = bb.xMax() - bb.xMin();
//	float bbH = bb.yMax() - bb.yMin();
//	//int texSizeX = 4096;
//	//int texSizeY = texSizeX * (bbH / bbW);
//
//	int texSizeY = 1280;
//	int texSizeX = texSizeY * (bbW / bbH);
//
//	float tileW = texSizeX * 0.25;
//	float tileH = texSizeY * 0.25;
//
//	float curX = bb.xMin();
//	float curY = bb.yMax();
//	int nrows = (int)(bbH / tileH)+1;
//	int ncols = (int)(bbW / tileW)+1;
//	//osg::Camera* camera2 = createScreenQuadCamera(node,bb,rtTexture);
//	//camera2->addChild(screenQuad);
//
//	//orthoCamera->getOrCreateStateSet()->setRenderBinDetails(2,"RenderBin"); 
//	osg::ref_ptr<osg::Texture2D> tex = new osg::Texture2D;
//
//	tex->setTextureSize(texSizeX, texSizeY);
//	tex->setResizeNonPowerOfTwoHint(false);
//	tex->setFilter(osg::Texture2D::MIN_FILTER,osg::Texture2D::LINEAR);
//	tex->setFilter(osg::Texture2D::MAG_FILTER,osg::Texture2D::LINEAR);
//	tex->setWrap(osg::Texture2D::WRAP_S,osg::Texture2D::REPEAT);
//	tex->setWrap(osg::Texture2D::WRAP_T,osg::Texture2D::REPEAT);
//	//rtTexture->setDataVariance(osg::Object::DYNAMIC);
//	tex->setInternalFormat(GL_RGBA);
//	tex->setSourceFormat(GL_RGBA);
//	tex->setSourceType(GL_UNSIGNED_BYTE);
//	/*tex->setInternalFormat(GL_ALPHA32F_ARB);
//	tex->setSourceFormat(GL_ALPHA);
//	tex->setSourceType(GL_FLOAT);*/
//	//tex->setInternalFormat(GL_RGBA32F_ARB);
//	//tex->setSourceFormat(GL_RGBA);
//	//tex->setSourceType(GL_FLOAT);
//	osg::ref_ptr<osg::Image> img = new osg::Image;
//	//img->allocateImage(texSize, texSize,1,GL_ALPHA,GL_FLOAT);
//    img->allocateImage(texSizeX, texSizeY,1,GL_RGBA,GL_UNSIGNED_BYTE); 
//    //img->allocateImage(texSize, texSize,1,GL_RGBA,GL_FLOAT); 
//	tex->setImage(img.get());
//	osg::StateSet* stateset = sceneNode->getOrCreateStateSet();
//	osg::ref_ptr<osg::Program> program = new osg::Program;
//	char vertexShaderSource[] = 
//		"varying vec4 pos;\n"
//		"uniform vec4 bound;\n"
//		"void main(void)\n"
//		"{\n"
//		"gl_TexCoord[0] = gl_MultiTexCoord0;\n"
//		"pos = gl_Vertex;\n"
//		"vec4 min = vec4(bound.x,bound.z,0,1);\n"
//		"vec4 max = vec4(bound.y,bound.w,0,1);\n"
//		"gl_TexCoord[1] = vec4((pos.x-min.x)/(max.x-min.x),(pos.y-min.y)/(max.y-min.y),0,1);\n"
//		"gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\n"
//		"}\n";
//	char fragmentShaderSource[] = 
//
//		"uniform vec4 color;\n"
//		"uniform sampler2D tex;\n"
//		"uniform sampler2D texMask;\n"
//		"varying vec4 pos;\n"
//		"void main(void) \n"
//		"{\n"
//		"//gl_FragColor = texture2D(tex, gl_TexCoord[0].xy);\n"
//		"gl_FragColor = gl_TexCoord[1];\n"
//		"if(gl_TexCoord[1].y < 0.5)\n"
//		"   gl_FragColor = vec4(0,0,0,1); \n"
//
//
//		"}\n";
//
//	program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));
//	program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexShaderSource));
//	/*osg::ref_ptr<osg::Image> maskImg = osgEarth::URI(baseDir+"mask.png").getImage();
//	if(!maskImg || !maskImg.valid())
//	{
//		maskImg = new osg::Image;
//		maskImg->allocateImage(1,1,1,GL_RGBA,GL_UNSIGNED_BYTE);
//		unsigned char* data = maskImg->data();
//		data[0]=255;data[1]=255;data[2]=255;data[3]=255;
//
//	}*/
//	stateset->setAttribute(program.get(),osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE); 
//	visitor.reset();
//	node->accept(visitor);
//	osg::BoundingBox bound = visitor.getBoundingBox();
//
//	stateset->addUniform(new osg::Uniform("bound", osg::Vec4(bound.xMin(), bound.xMax(), bound.yMin(), bound.yMax())));
//	osg::ref_ptr<osg::Texture2D> maskTex = new osg::Texture2D;  
//	maskTex->setName( "maskTex" );
//	maskTex->setFilter(osg::Texture::MIN_FILTER, osg::Texture::LINEAR_MIPMAP_LINEAR);//_MIPMAP_LINEAR
//	maskTex->setFilter(osg::Texture::MAG_FILTER, osg::Texture::LINEAR);
//	//maskTex->setImage(maskImg.get());
//	//stateset->addUniform( new osg::Uniform("osg_ViewMatrixInverse", osg::Matrix::identity()) );
//	stateset->addUniform( new osg::Uniform("tex", 0) );
//	//stateset->addUniform( new osg::Uniform("texPolygon", 1) );
//	//stateset->addUniform( new osg::Uniform("texMask", 2) );
//	//stateset->addUniform( new osg::Uniform("bound", osg::Vec4(bound.xMin(),bound.xMax(),bound.yMin(),bound.yMax())) );
//	//stateset->setTextureAttributeAndModes( 2, maskTex.get(), osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);
//	//stateset->setTextureAttributeAndModes( 1, maskTex.get(), osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE );
//	//stateset->setTextureAttributeAndModes( 0, maskTex.get(), osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE );
//
//	//mat->setCullingActive(false);
//	//stateset->setAttribute(program.get(),osg::StateAttribute::ON |osg::StateAttribute::OVERRIDE);  
//	//osg::ClampColor* clamp = new osg::ClampColor();
//	//clamp->setClampVertexColor(GL_FALSE);
//	//clamp->setClampFragmentColor(GL_FALSE);
//	//clamp->setClampReadColor(GL_FALSE);
//	//stateset->setAttribute(clamp, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
//	osg::Camera* orthoCamera = new osg::Camera;
//	orthoCamera->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
//	orthoCamera->setClearMask( GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT );
//	orthoCamera->setClearColor( osg::Vec4(0,0,0,1));
//
//	//camera->setClearColor(osg::Vec4(0.53f, 0.85f, 1.0f, 0.9f));				// Background
//	orthoCamera->setReferenceFrame( osg::Transform::ABSOLUTE_RF_INHERIT_VIEWPOINT );
//	orthoCamera->setViewport( 0,0, texSizeX, texSizeY);
//	orthoCamera->setRenderOrder(osg::Camera::PRE_RENDER);
//	orthoCamera->setRenderTargetImplementation(osg::Camera::FRAME_BUFFER_OBJECT); 
//	orthoCamera->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//	viewer.setUpViewInWindow(0,0,texSizeX,texSizeY);
//	orthoCamera->attach(osg::Camera::COLOR_BUFFER0,tex.get());
//	orthoCamera->attach(osg::Camera::COLOR_BUFFER0,img.get());
//	orthoCamera->addChild(sceneNode);
//	viewer.setSceneData(orthoCamera);
//	viewer.realize();
//	viewer.getCamera()->setClearColor( osg::Vec4(0,0,0,0));
//	//viewer.setCameraManipulator(new osgGA::TrackballManipulator);
//	viewer.getCamera()->getOrCreateStateSet()->setMode( GL_BLEND, osg::StateAttribute::OFF|osg::StateAttribute::OVERRIDE);
//	//orthoCamera->getOrCreateStateSet()->setAttribute(clamp, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
//	STPager* databasePager = new STPager;
//	databasePager->cancel();
//	viewer.getScene()->setDatabasePager(databasePager);
//
//	//return 0;
//	GDALAllRegister();
//	const char *pszFormat = "GTiff";
//	GDALDriver *poDriver;
//	char **papszMetadata;
//
//	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
//	char **papszOptions = NULL;
//	double adfGeoTransform[6];
//	osg::Vec3d translate(0,0,0);
//	//float* tmpBuf = new float[texSize*texSize];
//	//for(int i=0;i<texSize*texSize;i++)
//	//{
//	//	tmpBuf[i] = 1000;
//	//}
//
//	//updateCameraMatrixOrtho(viewer.getCamera(), tileBB, 90);
//	//updateCameraMatrixOrtho(orthoCamera, tileBB, 90);
//	//viewer.frame();
//	//for (int n=0;n<50;n++)
//	//{
//	//	databasePager->frame();
//	//    viewer.frame();
//	//} 
//	setupObliqueCamera(orthoCamera, bb, 45);
//	viewer.frame();
//	databasePager->frame();
//
//	//for (int n=0;n<50;n++)
//	//{
//	//	databasePager->frame();
//	//    viewer.frame();
//	//} 
//	while (databasePager->getFileRequestListSize() > 0)
//	{
//		viewer.frame();
//		databasePager->frame();
//	}
//	viewer.frame();
//	osgDB::writeImageFile(*img, "oblique.png");
//	exit(0);
//	return 0;
//	for (int irow=0;irow<nrows;irow++)
//	{
//		for (int icol=0;icol<ncols;icol++)
//		{
//			curX = bb.xMin() + icol * tileW;
//			curY = bb.yMax() - irow * tileH;
//			tileBB = osg::BoundingBox(curX,curY-tileH,bb.zMin(),curX+tileW,curY,bb.zMax());
//			updateCameraMatrixOrtho(viewer.getCamera(),tileBB,90);
//			updateCameraMatrixOrtho(orthoCamera,tileBB,90);
//            //viewer.frame();
//			//for (int n=0;n<50;n++)
//			//{
//			//	databasePager->frame();
//			//    viewer.frame();
//			//} 
//			viewer.frame();
//			databasePager->frame();
//
//			//for (int n=0;n<50;n++)
//			//{
//			//	databasePager->frame();
//			//    viewer.frame();
//			//} 
//			while(databasePager->getFileRequestListSize() > 0)
//			{
//				viewer.frame();
//				databasePager->frame();
//			}
//	   		viewer.frame();     
//
//			std::stringstream ssSrc;
//			ssSrc << outdir.data() << irow << "_" << icol << ".png";
//			osgDB::writeImageFile(*img, ssSrc.str());
//
//
//
//			std::stringstream ssDest;
//			ssDest << outdir << irow << "_" << icol << ".tif";
//			//g_pWarpIntance->copyImageMemory();
//			//osgDB::writeImageFile(*img,ss.str());
//			GDALDataset *poSrcDS = 
//				(GDALDataset *) GDALOpen(ssSrc.str().data(), GA_ReadOnly );
//			if(!poSrcDS)
//				continue;
//
//			adfGeoTransform[0] = curX+translate.x();///* top left x */
//			adfGeoTransform[1] = (double)tileW/(double)texSizeX;///* w-e pixel resolution */
//			adfGeoTransform[2] = 0;///* 0 */
//			adfGeoTransform[3] = curY+translate.y();// /* top left y */
//			adfGeoTransform[4] = 0;///* 0 */
//			adfGeoTransform[5] =-(double)tileH/(double)texSizeY;///* n-s pixel resolution (negative value) */
//			//GDALDataset *poDstDS = poDriver->Create( ssDest.str().data(),texSize,texSize,1,GDT_Float32, NULL );
//			GDALDataset *poDstDS = poDriver->CreateCopy(ssDest.str().data(), poSrcDS, FALSE,
//				NULL, NULL, NULL);
//			poDstDS->SetGeoTransform(adfGeoTransform);
//			//poDstDS->SetProjection( pszSRS_WKT );
//			//double adfGeoTransform[6] = { 444720, 30, 0, 3751320, 0, -30 };
//			//osg::ref_ptr<osg::Image> newImg = new osg::Image(*img,osg::CopyOp::DEEP_COPY_ALL);
//			//newImg->flipVertical();
//			//float* buf = (float*)newImg->data();
//			/*		for (int k=0;k<texSize*texSize;k++)
//			{
//			if(buf[k] > 10)
//			printf("%f\n",buf[k]);
//			}*/
//			//poDstDS->GetRasterBand(1)->SetNoDataValue(0);
//			////poDstDS->RasterIO(GF_Write,0, 0,texSize, texSize, 
//				//buf, texSize, texSize, GDT_Float32,1,NULL,sizeof(float) , sizeof(float) * texSize,0);
//			////if( poDstDS != NULL )
//			//poDstDS->FlushCache();
//			GDALClose( (GDALDatasetH) poDstDS );
//			GDALClose( (GDALDatasetH) poSrcDS );
//			printf("%d/%d\n",irow,nrows);
//            //img->flipVertical();
//
//			//printf("%d/%d\n",irow,nrows);
//		}
//	}
//
//	//OGRSpatialReference oSRS;
//	//char *pszSRS_WKT = NULL;
//	//oSRS.SetUTM( 50, TRUE );
//	//oSRS.SetWellKnownGeogCS( "WGS84" );
//	//oSRS.exportToWkt( &pszSRS_WKT );
//
//	//for (int irow=0;irow<nrows;irow++)
//	//{
//	//	for (int icol=0;icol<ncols;icol++)
//	//	{
//	//		curX = bb.xMin() + icol * tileW;
//	//		curY = bb.yMax() - irow * tileH;
//
//	//		std::stringstream ssSrc;
//	//		ssSrc << outdir << irow << "_" << icol << ".jpg";
//	//		std::stringstream ssDest;
//	//		ssDest << outdir << irow << "_" << icol << ".tif";
//	//		//g_pWarpIntance->copyImageMemory();
//	//		//osgDB::writeImageFile(*g_pWarpedTexture->getImage(),ss.str());
//	//		GDALDataset *poSrcDS = 
//	//			(GDALDataset *) GDALOpen(ssSrc.str().data(), GA_ReadOnly );
//	//		if(!poSrcDS)
//	//			continue;
//
//	//		adfGeoTransform[0] = curX+translate.x();///* top left x */
//	//		adfGeoTransform[1] = (double)tileW/(double)texSize;///* w-e pixel resolution */
//	//		adfGeoTransform[2] = 0;///* 0 */
//	//		adfGeoTransform[3] = curY+translate.y();// /* top left y */
//	//		adfGeoTransform[4] = 0;///* 0 */
//	//		adfGeoTransform[5] =-(double)tileH/(double)texSize;///* n-s pixel resolution (negative value) */
//	//		GDALDataset *poDstDS = poDriver->CreateCopy( ssDest.str().data(), poSrcDS, FALSE, 
//	//			NULL, NULL, NULL );
//	//		poDstDS->SetGeoTransform(adfGeoTransform);
//	//		//poDstDS->SetProjection( pszSRS_WKT );
//	//		//double adfGeoTransform[6] = { 444720, 30, 0, 3751320, 0, -30 };
//	//		if( poDstDS != NULL )
//	//			GDALClose( (GDALDatasetH) poDstDS );
//	//		GDALClose( (GDALDatasetH) poSrcDS );
//	//		printf("%d/%d\n",irow,nrows);
//	//	}
//	//}
//	/*while(!viewer.done())
//	{
//		updateCameraMatrixOrtho(viewer.getCamera(),bb);
//		databasePager->frame();
//		viewer.frame();
//	}*/
//	return 1;
//
//}
bool DeleteDirectory(const QString &path)
{
	if (path.isEmpty())
		return false;

	QDir dir(path);
	if (!dir.exists())
		return true;

	dir.setFilter(QDir::AllEntries | QDir::NoDotAndDotDot);
	QFileInfoList fileList = dir.entryInfoList();
	foreach(QFileInfo fi, fileList)
	{
		if (fi.isFile())
			fi.dir().remove(fi.fileName());
		else
			DeleteDirectory(fi.absoluteFilePath());
	}
	return dir.rmpath(dir.absolutePath());
}

//int main(int argc, char** argv)
//{
//	osg::DisplaySettings::instance()->setNumOfDatabaseThreadsHint(1);
//	osg::DisplaySettings::instance()->setNumOfHttpDatabaseThreadsHint(0);
//
//
//	osgViewer::Viewer viewer;
//	//viewer.getScene()->setDatabasePager()
//	// add the state manipulator
//	viewer.addEventHandler(new osgGA::StateSetManipulator(viewer.getCamera()->getOrCreateStateSet()));
//
//	// add the thread model handler
//	viewer.addEventHandler(new osgViewer::ThreadingHandler);
//
//	// add the window size toggle handler
//	viewer.addEventHandler(new osgViewer::WindowSizeHandler);
//
//	// add the stats handler
//	viewer.addEventHandler(new osgViewer::StatsHandler);
//
//
//	// add the record camera path handler
//	//viewer.addEventHandler(new osgViewer::RecordCameraPathHandler);
//
//	// add the LOD Scale handler
//	//viewer.addEventHandler(new osgViewer::LODScaleHandler);
//	QDir qoutdir("./img/");
//	DeleteDirectory(qoutdir.absolutePath());
//	qoutdir.mkpath(".");
//	// add the screen capture handler
//	//viewer.addEventHandler(new osgViewer::ScreenCaptureHandler);
//	//osg::Node* node = osgDB::readNodeFile("F:/Oblique_Photogrammetry/weihai/Data/all.desc");
//	osg::Node* node = osgDB::readNodeFile("E:/efficient-sparse-voxel-octrees-1.4/scenes/kechuang.ive");
//	//osg::Node* node = osgDB::readNodeFile("F:/Oblique_Photogrammetry/xinji/Data/all.desc");
//	//osg::Node* node = osgDB::readNodeFile("F:/Oblique_Photogrammetry/weihai/Data/all.desc");
//	//osg::Node* node = osgDB::readNodeFile("F:/Oblique_Photogrammetry/hongpeng/Data/all.desc");
//	//osg::Node* node = osgDB::readNodeFile("F:/Oblique_Photogrammetry/suzhou_gaoxin/Data/all.desc");
//	/*osg::Geode* node = new osg::Geode();
//	node->addDrawable(new osg::ShapeDrawable(new osg::Box(osg::Vec3(0.0f, 0.0f, 0.0f), 2 * 1)));*/
//	//node->addDrawable(new osg::ShapeDrawable(new osg::Sphere(osg::Vec3(0.0f, 0.0f, 0.0f), 2)));
//	std::stringstream ss_outdir;
//	std::string outdir = "F:/Oblique_Photogrammetry/weihai/DOM025/";
//	//transformGCPs("../jiashan/georeference.txt",outdir + "georeference.txt",obliqueMat);
//	//osg::Node* sceneNode = 
//	//osg::Node* node = osgDB::readNodeFile("D:/Projects/Oblique_Photogrammetry/xinji/Data/all.desc");
//	node->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
//	//node->getOrCreateStateSet()->setMode( GL_BLEND, osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);
//	osg::MatrixTransform* sceneNode = new osg::MatrixTransform;
//	sceneNode->addChild(node);
//	//sceneNode->setMatrix(osg::Matrix::rotate(osg::DegreesToRadians(90.0), osg::Vec3(0,0,1)));
//	osg::ComputeBoundsVisitor visitor;
//	sceneNode->accept(visitor);
//	osg::BoundingBox bb = visitor.getBoundingBox();
//	//osg::BoundingBox tileBB;
//
//	//double resolution = 0.5;
//	//float bbW = bb.xMax() - bb.xMin();
//	//float bbH = bb.yMax() - bb.yMin();
//	//int texSizeX = 4096;
//	//int texSizeY = texSizeX * (bbH / bbW);
//	int texSizeX = 512;
//	int texSizeY = 512;
//	//int texSizeX = texSizeY * (bbW / bbH);
//
//	/*float tileW = texSizeX * 0.25;
//	float tileH = texSizeY * 0.25;
//
//	float curX = bb.xMin();
//	float curY = bb.yMax();
//	int nrows = (int)(bbH / tileH) + 1;
//	int ncols = (int)(bbW / tileW) + 1;*/
//	//osg::Camera* camera2 = createScreenQuadCamera(node,bb,rtTexture);
//	//camera2->addChild(screenQuad);
//
//	//orthoCamera->getOrCreateStateSet()->setRenderBinDetails(2,"RenderBin"); 
//	osg::ref_ptr<osg::Texture2D> tex = new osg::Texture2D;
//
//	tex->setTextureSize(texSizeX, texSizeY);
//	tex->setResizeNonPowerOfTwoHint(false);
//	tex->setFilter(osg::Texture2D::MIN_FILTER, osg::Texture2D::LINEAR);
//	tex->setFilter(osg::Texture2D::MAG_FILTER, osg::Texture2D::LINEAR);
//	tex->setWrap(osg::Texture2D::WRAP_S, osg::Texture2D::REPEAT);
//	tex->setWrap(osg::Texture2D::WRAP_T, osg::Texture2D::REPEAT);
//	//rtTexture->setDataVariance(osg::Object::DYNAMIC);
//	tex->setInternalFormat(GL_RGBA);
//	tex->setSourceFormat(GL_RGBA);
//	tex->setSourceType(GL_UNSIGNED_BYTE);
//	/*tex->setInternalFormat(GL_ALPHA32F_ARB);
//	tex->setSourceFormat(GL_ALPHA);
//	tex->setSourceType(GL_FLOAT);*/
//	//tex->setInternalFormat(GL_RGBA32F_ARB);
//	//tex->setSourceFormat(GL_RGBA);
//	//tex->setSourceType(GL_FLOAT);
//	osg::ref_ptr<osg::Image> img = new osg::Image;
//	//img->allocateImage(texSize, texSize,1,GL_ALPHA,GL_FLOAT);
//	img->allocateImage(texSizeX, texSizeY, 1, GL_RGBA, GL_UNSIGNED_BYTE);
//	//img->allocateImage(texSize, texSize,1,GL_RGBA,GL_FLOAT); 
//	tex->setImage(img.get());
//	osg::StateSet* stateset = sceneNode->getOrCreateStateSet();
//	osg::ref_ptr<osg::Program> program = new osg::Program;
//	char vertexShaderSource[] =
//		"varying vec4 pos;\n"
//		"uniform vec4 bound;\n"
//		"void main(void)\n"
//		"{\n"
//		"gl_TexCoord[0] = gl_MultiTexCoord0;\n"
//		"pos = gl_Vertex;\n"
//		"vec4 min = vec4(bound.x,bound.z,0,1);\n"
//		"vec4 max = vec4(bound.y,bound.w,0,1);\n"
//		"gl_TexCoord[1] = vec4((pos.x-min.x)/(max.x-min.x),(pos.y-min.y)/(max.y-min.y),0,1);\n"
//		"gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\n"
//		"}\n";
//	char fragmentShaderSource[] =
//
//		"uniform vec4 color;\n"
//		"uniform sampler2D tex;\n"
//		"uniform sampler2D texMask;\n"
//		"varying vec4 pos;\n"
//		"void main(void) \n"
//		"{\n"
//		"//gl_FragColor = texture2D(tex, gl_TexCoord[0].xy);\n"
//		"gl_FragColor = gl_TexCoord[1];\n"
//		"if(gl_TexCoord[1].y < 0.5)\n"
//		"   gl_FragColor = vec4(0,0,0,1); \n"
//
//
//		"}\n";
//
//	program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));
//	program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexShaderSource));
//	/*osg::ref_ptr<osg::Image> maskImg = osgEarth::URI(baseDir+"mask.png").getImage();
//	if(!maskImg || !maskImg.valid())
//	{
//	maskImg = new osg::Image;
//	maskImg->allocateImage(1,1,1,GL_RGBA,GL_UNSIGNED_BYTE);
//	unsigned char* data = maskImg->data();
//	data[0]=255;data[1]=255;data[2]=255;data[3]=255;
//
//	}*/
//	//stateset->setAttribute(program.get(), osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
//	visitor.reset();
//	node->accept(visitor);
//	//osg::BoundingBox bound = visitor.getBoundingBox();
//
//	//stateset->addUniform(new osg::Uniform("bound", osg::Vec4(bound.xMin(), bound.xMax(), bound.yMin(), bound.yMax())));
//	osg::ref_ptr<osg::Texture2D> maskTex = new osg::Texture2D;
//	maskTex->setName("maskTex");
//	maskTex->setFilter(osg::Texture::MIN_FILTER, osg::Texture::LINEAR_MIPMAP_LINEAR);//_MIPMAP_LINEAR
//	maskTex->setFilter(osg::Texture::MAG_FILTER, osg::Texture::LINEAR);
//	//maskTex->setImage(maskImg.get());
//	//stateset->addUniform( new osg::Uniform("osg_ViewMatrixInverse", osg::Matrix::identity()) );
//	stateset->addUniform(new osg::Uniform("tex", 0));
//	//stateset->addUniform( new osg::Uniform("texPolygon", 1) );
//	//stateset->addUniform( new osg::Uniform("texMask", 2) );
//	//stateset->addUniform( new osg::Uniform("bound", osg::Vec4(bound.xMin(),bound.xMax(),bound.yMin(),bound.yMax())) );
//	//stateset->setTextureAttributeAndModes( 2, maskTex.get(), osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);
//	//stateset->setTextureAttributeAndModes( 1, maskTex.get(), osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE );
//	//stateset->setTextureAttributeAndModes( 0, maskTex.get(), osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE );
//
//	//mat->setCullingActive(false);
//	//stateset->setAttribute(program.get(),osg::StateAttribute::ON |osg::StateAttribute::OVERRIDE);  
//	//osg::ClampColor* clamp = new osg::ClampColor();
//	//clamp->setClampVertexColor(GL_FALSE);
//	//clamp->setClampFragmentColor(GL_FALSE);
//	//clamp->setClampReadColor(GL_FALSE);
//	//stateset->setAttribute(clamp, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
//	ObliqueCamera* orthoCamera = new ObliqueCamera(bb,45);
//	orthoCamera->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
//	orthoCamera->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//	orthoCamera->setClearColor(osg::Vec4(0, 0, 0, 1));
//
//	//camera->setClearColor(osg::Vec4(0.53f, 0.85f, 1.0f, 0.9f));				// Background
//	orthoCamera->setReferenceFrame(osg::Transform::ABSOLUTE_RF_INHERIT_VIEWPOINT);
//	orthoCamera->setViewport(0, 0, texSizeX, texSizeY);
//	orthoCamera->setRenderOrder(osg::Camera::PRE_RENDER);
//	orthoCamera->setRenderTargetImplementation(osg::Camera::FRAME_BUFFER_OBJECT);
//	orthoCamera->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//	viewer.setUpViewInWindow(0, 0, texSizeX, texSizeY);
//	orthoCamera->attach(osg::Camera::COLOR_BUFFER0, tex.get());
//	orthoCamera->attach(osg::Camera::COLOR_BUFFER0, img.get());
//	orthoCamera->addChild(sceneNode);
//	viewer.setSceneData(orthoCamera);
//	viewer.realize();
//	viewer.getCamera()->setClearColor(osg::Vec4(0, 0, 0, 0));
//	//viewer.setCameraManipulator(new osgGA::TrackballManipulator);
//	viewer.getCamera()->getOrCreateStateSet()->setMode(GL_BLEND, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
//	//orthoCamera->getOrCreateStateSet()->setAttribute(clamp, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
//	STPager* databasePager = new STPager;
//	databasePager->cancel();
//	viewer.getScene()->setDatabasePager(databasePager);
//
//	//return 0;
//	GDALAllRegister();
//	const char *pszFormat = "GTiff";
//	GDALDriver *poDriver;
//	char **papszMetadata;
//
//	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
//	char **papszOptions = NULL;
//	double adfGeoTransform[6];
//	orthoCamera->setupGrid(500);
//	for (int irow = 0; irow < orthoCamera->_nrows; irow++)
//	{
//		for (int icol = 0; icol < orthoCamera->_ncols; icol++)
//		{
//			orthoCamera->setupCameraAtCell(irow, icol);
//			while (databasePager->getFileRequestListSize() > 0)
//			{
//				viewer.frame();
//				databasePager->frame();
//			}
//			viewer.frame();
//			viewer.frame();
//	/*		viewer.frame();
//			viewer.frame();*/
//			orthoCamera->produceImage(irow, icol, img);
//		}
//	}
//	exit(0);
//}

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

	osg::ref_ptr<osg::Node> node = osgDB::readNodeFile(infile.data());

	node->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
	//node->getOrCreateStateSet()->setMode( GL_BLEND, osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);
	osg::ref_ptr<osg::MatrixTransform> sceneNode = new osg::MatrixTransform;
	sceneNode->addChild(node.get());
	osg::Matrix mat = osg::Matrix::rotate(osg::DegreesToRadians(azimuthAngle), osg::Vec3(0, 0, 1));
	sceneNode->setMatrix(mat);
	osg::ComputeBoundsVisitor visitor;
	sceneNode->accept(visitor);
	osg::BoundingBox bb = visitor.getBoundingBox();
	
	if (argc > 9)
	{
		double xmin, ymin, xmax, ymax;
		xmin = atof(argv[6]); xmax = atof(argv[7]);
		ymin = atof(argv[8]); ymax = atof(argv[9]);
		std::vector<osg::Vec3d> corners;
		corners.push_back(osg::Vec3d(xmin, ymin, bb.zMax()));
		corners.push_back(osg::Vec3d(xmin, ymax, bb.zMin()));
		corners.push_back(osg::Vec3d(xmax, ymin, bb.zMin()));
		corners.push_back(osg::Vec3d(xmax, ymax, bb.zMax()));
		osg::BoundingBox newbb;
		for (size_t i = 0; i < corners.size(); i++)
		{
			osg::Vec3d newCorner = corners[i] * mat;
			newbb.expandBy(newCorner);
		}

		bb = newbb;
	}
	int texSizeX = 1024;
	int texSizeY = 1024;
	
	osg::ref_ptr<osg::Texture2D> tex = new osg::Texture2D;

	tex->setTextureSize(texSizeX, texSizeY);
	tex->setResizeNonPowerOfTwoHint(false);
	tex->setFilter(osg::Texture2D::MIN_FILTER, osg::Texture2D::LINEAR);
	tex->setFilter(osg::Texture2D::MAG_FILTER, osg::Texture2D::LINEAR);
	tex->setWrap(osg::Texture2D::WRAP_S, osg::Texture2D::REPEAT);
	tex->setWrap(osg::Texture2D::WRAP_T, osg::Texture2D::REPEAT);
	//rtTexture->setDataVariance(osg::Object::DYNAMIC);
	tex->setInternalFormat(GL_RGBA);
	tex->setSourceFormat(GL_RGBA);
	tex->setSourceType(GL_UNSIGNED_BYTE);

	osg::ref_ptr<osg::Image> img = new osg::Image;
	//img->allocateImage(texSize, texSize,1,GL_ALPHA,GL_FLOAT);
	img->allocateImage(texSizeX, texSizeY, 1, GL_RGBA, GL_UNSIGNED_BYTE);
	//img->allocateImage(texSize, texSize,1,GL_RGBA,GL_FLOAT); 
	tex->setImage(img.get());

	//ObliqueCameraWrapper* orthoCamera = new ObliqueCameraWrapper(viewer.getCamera(),bb, 45);
	viewer.getCamera()->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
	viewer.getCamera()->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//viewer.getCamera()->setClearColor(osg::Vec4(0, 0, 0, 0));

	osg::ref_ptr<ObliqueCamera> orthoCamera = new ObliqueCamera(bb,altAngle);
	orthoCamera->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
	orthoCamera->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	orthoCamera->setClearColor(osg::Vec4(0, 0, 0, 0));
	
	orthoCamera->setReferenceFrame(osg::Transform::ABSOLUTE_RF_INHERIT_VIEWPOINT);
	orthoCamera->setViewport(0, 0, texSizeX, texSizeY);
	orthoCamera->setRenderOrder(osg::Camera::PRE_RENDER);
	orthoCamera->setRenderTargetImplementation(osg::Camera::FRAME_BUFFER_OBJECT);
	orthoCamera->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	viewer.setUpViewInWindow(0, 0, texSizeX, texSizeY);
	orthoCamera->attach(osg::Camera::COLOR_BUFFER0, tex.get());
	orthoCamera->attach(osg::Camera::COLOR_BUFFER0, img.get());
	orthoCamera->addChild(sceneNode.get());
	viewer.setSceneData(orthoCamera.get());
	viewer.setUpViewInWindow(0, 0, 1, 1);

	viewer.realize();
	//viewer.setCameraManipulator(new osgGA::TrackballManipulator);
	//viewer.getCamera()->getOrCreateStateSet()->setMode(GL_BLEND, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
	//orthoCamera->getOrCreateStateSet()->setAttribute(clamp, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
	STPager* databasePager = new STPager;
	databasePager->cancel();
	viewer.getScene()->setDatabasePager(databasePager);

	//return 0;
	GDALAllRegister();
	std::vector<std::string> tilefiles;
	//char **papszOptions = NULL;
	//double adfGeoTransform[6];
	orthoCamera->setupGrid(texSizeX*resol, texSizeX);
	for (int irow = 0; irow < orthoCamera->_nrows; irow++)
	{
		for (int icol = 0; icol < orthoCamera->_ncols; icol++)
		{
			//if (orthoCamera->imagExists(irow, icol, img, outdir))
			//	continue;
			orthoCamera->setupCameraAtCell(irow, icol);
			viewer.frame();
			databasePager->frame();
			while (databasePager->getFileRequestListSize() > 0)
			{
				viewer.frame();
				databasePager->frame();
			}
			viewer.frame();
			tilefiles.push_back(orthoCamera->produceImage(irow, icol, img,outdir));
		}
	}

	const char *pszFormat = "GTiff";
	GDALDriver *poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	std::stringstream commandss;
	commandss << "gdalbuildvrt " << outdir + "oblique.vrt " << outdir + "*.tif";
	system(commandss.str().data());


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


	exit(0);
}

