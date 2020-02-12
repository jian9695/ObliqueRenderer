
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

using namespace osgDB;
using namespace OpenThreads;
#include "ogrsf_frmts.h"
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
		poDS = OGRSFDriverRegistrar::Open(filename.data(), update);
		poLayer = poDS->GetLayer(0);
	}

	void close()
	{
		if (poDS)
			OGRDataSource::DestroyDataSource(poDS);
		poDS = NULL;
	}
	void create(const char* filename, OGRSpatialReference* spatialRef = NULL, OGRFeatureDefn *poFDefn = NULL, OGRwkbGeometryType geotype = wkbPoint)
	{


		g_mFileName = filename;
		if (poDS)
			OGRDataSource::DestroyDataSource(poDS);
		const char *pszDriverName = "ESRI Shapefile";
		OGRSFDriver *poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(
			pszDriverName);

		if (QFileInfo(filename).exists())
		{
			poDriver->DeleteDataSource(filename);
			//poDS = OGRSFDriverRegistrar::Open(filename, TRUE);
			//poLayer = poDS->GetLayer(0);
		}

		poDS = poDriver->CreateDataSource(filename, NULL);
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
	OGRDataSource   *poDS;
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
	FindCompileableGLObjectsVisitor(const DatabasePager* pager):
	  osgUtil::StateToCompile(osgUtil::GLObjectsVisitor::COMPILE_DISPLAY_LISTS|osgUtil::GLObjectsVisitor::COMPILE_STATE_ATTRIBUTES),
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
				DatabasePager::FindCompileableGLObjectsVisitor stateToCompile(this);
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
void updateCameraMatrix(osg::Camera* camera,const osg::BoundingBox& bb,double altAngle, double azimuthAngle)
{
	
	////清除颜色和深度缓存，这意味着这个相机渲染子场景时将会覆盖之前任何相机的渲染数据
		camera->setClearMask( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	//
	////设置相机的坐标系，设置为ABSOLUTE_RF意味着相机的所有变换矩阵和观察/投影矩阵设置都是相对于世界坐标的，不会受到上级矩阵的影响
		camera->setReferenceFrame( osg::Transform::ABSOLUTE_RF );
		double radius = osg::BoundingSphere(bb).radius();
	//根据模型的包围球大小，设置相机的投影矩阵
		double viewDistance = 2.0 * radius;
		double znear = viewDistance - radius*2;
		double zfar = viewDistance + radius*2;
		float top = radius;
		float right = radius;
	//投影变换，设置正射投影矩阵的内容，6个参数分别表示平行视景体的6个面的位置（左，右，下，上，近，远）
	//还有一个透射投影。正射投影与透射投影的区别：
	//先说一下相机吧，包括3部分：1、物体在那儿，观察者所处的位置及观察角度。
	//2、将看到的投影到底片上，那么有两种方式正射---底片上的物体投影与实际一样大；透射---远小近大，
	//将视锥体之外的视景剪切掉，远大近小与参数far-near的值有关，视锥体外视景的裁剪与近裁切面
	//（除了far以外的参数有关）有关，根据本例，视景正好都在近裁切平面内，故不会有视景被裁掉。
	//物体在近裁剪面上（近裁剪面与物体之间还可能有距离，是这样，这两者之间为正射投影，能不能包含视
	//景全部，就看近裁剪平面的大小）.
	//3，将底片上的景物投射到屏幕上。
		camera->setProjectionMatrixAsOrtho( -right, right, -top, top, znear, zfar );
		//设置为正看,就是对应眼睛看的方向，人的头的向上的方向。
	   osg::Vec3d upDirection( 0.0,0.0,1.0 );
	  //设置眼睛所在位置，看物体的方向
	   osg::Vec3d viewDirection( 0.0,-1.0,0.0 );

	   osg::Vec3d referenceAxis(-1.0,0.0,0.0 ); 

	   osg::Matrixd rotateMatrix = osg::Matrixd::rotate(osg::DegreesToRadians(azimuthAngle), osg::Vec3d( 0.0,0.0,1.0 ));
	   viewDirection = viewDirection * rotateMatrix;
	   referenceAxis = referenceAxis * rotateMatrix;
	   rotateMatrix = osg::Matrixd::rotate(osg::DegreesToRadians(altAngle),referenceAxis);
	   viewDirection = viewDirection * rotateMatrix;
	   upDirection = -viewDirection ^ referenceAxis;

	   upDirection.normalize();
	   viewDirection.normalize();

	  //设置物体的位置
		osg::Vec3d center = bb.center();

	  //设置人眼所在的位置
		osg::Vec3d eyePoint = center + viewDirection * viewDistance;

	  //视点变换函数
		camera->setViewMatrixAsLookAt( eyePoint, center, upDirection );
	

}
void updateCameraMatrixOrtho(osg::Camera* camera, const osg::BoundingBox& bb,double altAngle=0)
{

	//

	//根据模型的包围球大小，设置相机的投影矩阵
	double viewDistance = 5.0 * bb.radius();
	double znear = viewDistance - bb.radius();
	double zfar = viewDistance + bb.radius();
	//float top   = (bb.yMax()-bb.yMin())*0.5;
	//float right = (bb.xMax()-bb.xMin());
	double xsize = (bb.yMax() - bb.yMin())*0.5;
	double ysize = ( bb.xMax() - bb.xMin())*0.5;
	double radius = xsize;
	if(radius < ysize)
		radius = ysize;
	float top = xsize * sin(osg::DegreesToRadians(altAngle));
	float right = ysize;
    camera->setReferenceFrame( osg::Transform::ABSOLUTE_RF );
	//camera->setProjectionMatrixAsOrtho( bb.yMin(), bb.yMax(),bb.xMin(), bb.xMax(),  znear, zfar );
    camera->setProjectionMatrixAsOrtho( -right, right, -top, top, znear, zfar );
	//设置为正看,就是对应眼睛看的方向，人的头的向上的方向。
	osg::Vec3d upDirection( 0.0,1.0,0.0 );
	//设置眼睛所在位置，看物体的方向
	osg::Vec3d viewDirection( 0.0,0.0,1.0 );

	//设置物体的位置
	osg::Vec3d center = bb.center();

	//设置人眼所在的位置
	osg::Vec3d eyePoint = center + viewDirection * viewDistance;

	//osg::Matrix mvp = osg::Matrix::lookAt( eyePoint, center, upDirection ) * osg::Matrix::ortho( -right, right, -top, top, znear, zfar );
    camera->setViewMatrixAsLookAt( eyePoint, center, upDirection );


	////设置为正看,就是对应眼睛看的方向，人的头的向上的方向。
	//osg::Vec3d upDirection( 0.0,0.0,1.0 );
	////设置眼睛所在位置，看物体的方向
	//osg::Vec3d viewDirection( 0.0,-1.0,0.0 );

	//osg::Vec3d referenceAxis(-1.0,0.0,0.0 ); 

	//osg::Matrixd rotateMatrix = osg::Matrixd::rotate(osg::DegreesToRadians(0.0), osg::Vec3d( 0.0,0.0,1.0 ));
	//viewDirection = viewDirection * rotateMatrix;
	//referenceAxis = referenceAxis * rotateMatrix;
	//rotateMatrix = osg::Matrixd::rotate(osg::DegreesToRadians(altAngle),referenceAxis);
	//viewDirection = viewDirection * rotateMatrix;
	//upDirection = -viewDirection ^ referenceAxis;

	//upDirection.normalize();
	//viewDirection.normalize();

	////设置物体的位置
	//osg::Vec3d center = bb.center();

	////设置人眼所在的位置
	//osg::Vec3d eyePoint = center + viewDirection * viewDistance;

	////视点变换函数
	//camera->setViewMatrixAsLookAt( eyePoint, center, upDirection );
}
#include <gdal_priv.h>
#include <ogr_spatialref.h>
void transformGCPs(const std::string infile,const std::string outfile,const osg::Matrix& mat)
{
	std::ifstream ifs(infile.data());

	std::ofstream ofs(outfile.data());	
	char buf[1024];
	while(ifs.getline(buf,1024))
	{
		std::stringstream str(buf);
		double x1,y1,x2,y2;
		str >> x1 >> y1 >> x2 >> y2;
		osg::Vec3d pos(x1,y1,0);
		pos = pos * mat;
		x1 = pos.x();
		y1 = pos.y();
		ofs << x1 << "	" << y1 <<  "	" << x2 <<  "	" <<  y2 << std::endl;

	}
	ifs.close();
	ofs.flush();
	ofs.close();
}

void printfVec3(osg::Vec3 p1, osg::Vec3 p2)
{
	printf("(%f,%f,%f),(%f,%f,%f)\n", p1.x(), p1.y(), p1.z(), p2.x(), p2.y(), p2.z());
}

#include "osg/LineSegment"
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
	double _left;
	double _right;
	double _bottom;
	double _top;
	double _znear;
	double _zfar;
	int _ncols;
	int _nrows;
	double _cellsize;
	osg::BoundingBox _bb;
	osg::BoundingBox _bbOrigial;
	GDALDriver *_poDriver;
	int _curRow;
	int _curCol;
	ObliqueCamera(const osg::BoundingBox& bb, double altAngle = 0,  double azimuthAngle = 0)
	{
		_bbOrigial = bb;
		float angle = osg::DegreesToRadians(altAngle);

		osg::Vec3d bbmin(bb.xMin(), bb.yMin(), bb.zMin());
		osg::Vec3d bbmax(bb.xMax(), bb.yMax(), bb.zMin());
		osg::Vec3d bbsize = bbmax - bbmin;

		float maxsize = bbsize.x() > bbsize.y() ? bbsize.x() : bbsize.y();
		float viewDistance = maxsize * 5;

		osg::Vec3d translation(0, -viewDistance * cos(angle), viewDistance * sin(angle));
		if (altAngle == 0.1)
		{
			translation = osg::Vec3d(0, 0, viewDistance);
		}
		osg::Vec3d bbminNear = bbmin + translation;
		osg::Vec3d bbmaxNear = bbmax + translation;
		osg::Vec3d bbminFar = bbmin - translation;
		osg::Vec3d bbmaxFar = bbmax - translation;
		double znear = viewDistance - bb.radius();
		double zfar = viewDistance + bb.radius();
		float  right = (bb.xMax() - bb.xMin())*0.5;

		this->setReferenceFrame(osg::Transform::ABSOLUTE_RF);

		osg::Vec3d upDirection(0.0, 1.0, 0.0);
		upDirection.normalize();
		osg::Vec3d viewDirection = -translation;
		viewDirection.normalize();
		osg::Vec3d center = (bbmin + bbmax) * 0.5;

		osg::Vec3d eyePoint = center - viewDirection * viewDistance;

		this->setViewMatrixAsLookAt(eyePoint, center, upDirection);


		float height = bb.zMax() - bb.zMin();
		float part1 = height / tan(angle);
		float part2 = (bb.yMax() - bb.yMin()) * 0.5;
		float top = sin(angle) * (part1 + part2);
		float bottom = sin(angle) * part2;

		this->setProjectionMatrixAsOrtho(-right, right, -bottom, top, znear, zfar);

		_azimuthAngle = osg::DegreesToRadians(azimuthAngle);
		_altAngle = osg::DegreesToRadians(altAngle);
		_viewDirection = viewDirection;
		_center = center;
		_left = -right;
		_right = right;
		_top = top;
		_bottom = -bottom;
		_znear = znear;
		_zfar = zfar;
		_viewDistance = viewDistance;
		_upDirection = upDirection;

		_bb = bb;
		_center.y() = _center.y() + part1 * 0.5;

		//return 0;
		GDALAllRegister();
		const char *pszFormat = "GTiff";
		
		char **papszMetadata;

		_poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
		char **papszOptions = NULL;
		//std::stringstream ssSrc;
		//ssSrc << nrow << "_" << ncol << ".tif";

	}

	void setupGrid(double cellsize)
	{
		//double width = _right - _left;
		//double height = _top - _bottom;
		double width = _bb.xMax() - _bb.xMin();
		double height = _bb.yMax() - _bb.yMin();
		float margin = (_bb.zMax() - _bb.zMin()) / tan(_altAngle) * sin(_altAngle);

		height += margin * 2;

		_nrows = height / cellsize;
		while (_nrows * cellsize < height)
			_nrows = _nrows + 1;

		_ncols = width / cellsize;
		while (_ncols * cellsize < width)
			_ncols = _ncols + 1;


		_cellsize = cellsize;

		//OGREnvelope newbound;
		_left = _left -  ((_ncols * cellsize) - width) * 0.5;
		_right = _right+ ((_ncols * cellsize) - width) * 0.5;
		_bottom = _bottom - ((_nrows * cellsize) - height) * 0.5;
		_top = _top + ((_nrows * cellsize) - height) * 0.5;

		_bb.xMin() = _center.x() - width * 0.5;
		_bb.xMax() = _center.x() + width * 0.5;
		_bb.yMin() = _center.y() - height * 0.5;
		_bb.yMax() = _center.y() + height * 0.5;


	}

	void setupCameraAtPosition(osg::Vec3d pos)
	{
		osg::Vec3d l1 = pos + _viewDirection * _viewDistance;
		osg::Vec3d l2 = pos - _viewDirection * _viewDistance;
		osg::ref_ptr<osg::LineSegment> line = new osg::LineSegment;
		line->set(l1, l2);
		osg::Vec3d center = _center;
		center.z() = _bb.zMax();
		osg::Vec3d p1 = center + (osg::Vec3d(_bb.xMin(), _bb.yMin(), _bb.zMin()) - _center) * 10;
		osg::Vec3d p2 = center + (osg::Vec3d(_center.x(), _bb.yMax(), _bb.zMin()) - _center) * 10;
		osg::Vec3d p3 = center + (osg::Vec3d(_bb.xMax(), _bb.yMin(), _bb.zMin()) - _center) * 10;
		double ratio;
	    line->intersect(p1, p2, p3,ratio);
		pos = line->start() + (line->end() - line->start()) * ratio;
		//center.x() = _bb.xMin() + ncol * _cellsize + _cellsize * 0.5;
		//center.y() = _bb.yMax() - nrow * _cellsize - _cellsize * 0.5;
		int nrow = (int)((_bb.yMax() - pos.y()) / _cellsize);
		int ncol = (int)((pos.x() - _bb.xMin()) / _cellsize);
		//osg::Vec3d p1 = _center + (_bb.xmin())
	    //line.intersect
	/*	osg::Vec3d center = _center;*/
		center.z() = _bb.zMin();
		center.x() = _bb.xMin() + ncol * _cellsize + _cellsize * 0.5;
		center.y() = _bb.yMax() - nrow * _cellsize - _cellsize * 0.5;
		osg::Vec3d eyePoint = center - _viewDirection * _viewDistance;

		//printf("%f,%f,%f: %f,%f,%f\n", center.x(), center.y(), center.z(), eyePoint.x(), eyePoint.y(), eyePoint.z());
		double halfsize = _cellsize * 0.5;
		float left = -halfsize;
		float right = halfsize;
		float bottom = -sin(_altAngle) * halfsize;
		float top = sin(_altAngle) * halfsize;
		float height = _bb.zMax() - _bb.zMin();

		this->setViewMatrixAsLookAt(eyePoint, center, _upDirection);
		double edge =  _cellsize * 2;
		this->setProjectionMatrixAsOrtho(left - edge, right + edge, bottom - edge, top + edge, _znear, _zfar);

		_curRow = nrow;
		_curCol = ncol;
	}
	void getRowCol(osg::Vec3d pos,int& nrow,int& ncol)
	{
		osg::Vec3d l1 = pos + _viewDirection * _viewDistance;
		osg::Vec3d l2 = pos - _viewDirection * _viewDistance;
		osg::ref_ptr<osg::LineSegment> line = new osg::LineSegment;
		line->set(l1, l2);
		osg::Vec3d center = _center;
		center.z() =  _bb.zMin();

		osg::Vec3d p1 = center + (osg::Vec3d(_bb.xMin(), _bb.yMin(),  _bb.zMin()) - _center) * 10;
		osg::Vec3d p2 = center + (osg::Vec3d(_center.x(), _bb.yMax(), _bb.zMin()) - _center) * 10;
		osg::Vec3d p3 = center + (osg::Vec3d(_bb.xMax(), _bb.yMin(),  _bb.zMin()) - _center) * 10;
		double ratio;
		line->intersect(p1, p2, p3, ratio);
		pos = line->start() + (line->end() - line->start()) * ratio;
		
		nrow = (int)((_bb.yMax() - pos.y()) / _cellsize);
		ncol = (int)((pos.x() - _bb.xMin()) / _cellsize);
	}
	void setupCameraAtCell(int nrow, int ncol)
	{
		osg::Vec3d center = _center;
		center.x() = _bb.xMin() + ncol * _cellsize + _cellsize * 0.5;
		center.y() = _bb.yMax() - nrow * _cellsize - _cellsize * 0.5;
		osg::Vec3d eyePoint = center - _viewDirection * _viewDistance;

		//printf("%f,%f,%f: %f,%f,%f\n", center.x(), center.y(), center.z(), eyePoint.x(), eyePoint.y(), eyePoint.z());
		double halfsize = _cellsize * 0.5;
		float left = -halfsize;
		float right = halfsize;
		float bottom = -sin(_altAngle) * halfsize;
		float top = sin(_altAngle) * halfsize;
		float height = _bb.zMax() - _bb.zMin();
		//float part1 = height / tan(_altAngle);
		//float part2 = halfsize;
		//if (nrow == 0)
		//{
		//	top = (part1 + part2) * sin(_altAngle);
		//}
		//if (nrow == _nrows - 1)
		//{
		//	bottom = -(part1 + part2) * sin(_altAngle);
		//}

		this->setViewMatrixAsLookAt(eyePoint, center, _upDirection);

		this->setProjectionMatrixAsOrtho(left, right, bottom, top, _znear, _zfar);
		_curRow = nrow;
		_curCol = ncol;
	}

	void produceImage(int nrow, int ncol,osg::Image* img,std::string dir = "./img/")
	{
		QDir qdir(dir.data());
		if (!qdir.exists())
			qdir.mkpath(".");
		std::stringstream ssSrc;
		ssSrc << dir << nrow << "_" << ncol << "_2.tif";
		std::stringstream ssDest;
		ssDest << dir << nrow << "_" << ncol << ".tif";

		osgDB::writeImageFile(*img, ssSrc.str().data());
		GDALDataset *poSrcDS = (GDALDataset *)GDALOpen(ssSrc.str().data(), GA_ReadOnly);


		//GDALDataset *poDstDS = (GDALDataset *)GDALOpen(ssDest.str().data(), GA_Update);
		double adfGeoTransform[6];
		adfGeoTransform[0] = _bb.xMin() + ncol * _cellsize;///* top left x */
		adfGeoTransform[1] = _cellsize / img->s();///* w-e pixel resolution */
		adfGeoTransform[2] = 0;///* 0 */
		adfGeoTransform[3] = _bb.yMax() - nrow * _cellsize;// /* top left y */
		adfGeoTransform[4] = 0;///* 0 */
		adfGeoTransform[5] = -_cellsize / img->t();///* n-s pixel resolution (negative value) */
														   


		char **papszOptions = NULL;

		papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");
		papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "LZW");
		GDALDataset * poDstDS = _poDriver->CreateCopy(ssDest.str().data(), poSrcDS, FALSE,
			papszOptions, GDALTermProgress, NULL);

		poDstDS->SetGeoTransform(adfGeoTransform);
		//GDALClose((GDALDatasetH)poDstDS);

		/* Once we're done, close properly the dataset */
		if (poDstDS != NULL)
			GDALClose((GDALDatasetH)poDstDS);
		CSLDestroy(papszOptions);
		GDALClose((GDALDatasetH)poSrcDS);
		QFile::remove(ssSrc.str().data());

	}
	void produceImage(int id, osg::Image* img)
	{

		std::stringstream ssSrc;
		ssSrc << "img/" << id  << "_2.tif";
		std::stringstream ssDest;
		ssDest << "img/" << id << ".tif";

		osgDB::writeImageFile(*img, ssSrc.str().data());
		GDALDataset *poSrcDS = (GDALDataset *)GDALOpen(ssSrc.str().data(), GA_ReadOnly);


		//GDALDataset *poDstDS = (GDALDataset *)GDALOpen(ssDest.str().data(), GA_Update);
		double adfGeoTransform[6];
		adfGeoTransform[0] = _bb.xMin() + _curCol * _cellsize;///* top left x */
		adfGeoTransform[1] = _cellsize / img->s();///* w-e pixel resolution */
		adfGeoTransform[2] = 0;///* 0 */
		adfGeoTransform[3] = _bb.yMax() - _curRow * _cellsize;// /* top left y */
		adfGeoTransform[4] = 0;///* 0 */
		adfGeoTransform[5] = -_cellsize / img->t();///* n-s pixel resolution (negative value) */



		char **papszOptions = NULL;

		papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");
		papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "LZW");
		GDALDataset * poDstDS = _poDriver->CreateCopy(ssDest.str().data(), poSrcDS, FALSE,
			papszOptions, GDALTermProgress, NULL);

		poDstDS->SetGeoTransform(adfGeoTransform);
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

		//GDALDataset *poDstDS = (GDALDataset *)GDALOpen(ssDest.str().data(), GA_Update);
		double adfGeoTransform[6];
		double halfsize = _cellsize * 0.5;
		//float left = -halfsize;
		//float right = halfsize;
		//float bottom = -sin(_altAngle) * halfsize;
		//float top = sin(_altAngle) * halfsize;
		//float height = _bb.zMax() - _bb.zMin();
		////float part1 = height / tan(_altAngle);
		////float part2 = halfsize;
		////if (nrow == 0)
		////{
		////	top = (part1 + part2) * sin(_altAngle);
		////}
		////if (nrow == _nrows - 1)
		////{
		////	bottom = -(part1 + part2) * sin(_altAngle);
		////}

		//this->setViewMatrixAsLookAt(eyePoint, center, _upDirection);

		//this->setProjectionMatrixAsOrtho(left, right, bottom, top, _znear, _zfar);

		adfGeoTransform[0] = _bb.xMin() + _curCol * _cellsize;///* top left x */
		adfGeoTransform[1] = _cellsize / img->s();///* w-e pixel resolution */
		adfGeoTransform[2] = 0;///* 0 */
		adfGeoTransform[3] = _bb.yMax() - _curRow * _cellsize;// /* top left y */
		adfGeoTransform[4] = 0;///* 0 */
		adfGeoTransform[5] = -_cellsize / img->t();///* n-s pixel resolution (negative value) */


		std::vector<osg::Vec3d> pixels;
		unsigned char* data = img->data();
		for (size_t i = 0; i < img->t(); i++)
		{
			double y = adfGeoTransform[3] + adfGeoTransform[5] * (img->t()-i-1) + adfGeoTransform[5] * 0.5;
			for (size_t j = 0; j < img->s(); j++)
			{
				double x = adfGeoTransform[0] + adfGeoTransform[1] * j + adfGeoTransform[1] * 0.5;
				if(*data > 0.5)
				   pixels.push_back(osg::Vec3d(x, y, 0));
				data+=4;
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
//class ObliqueCameraWrapper
//{
//public:
//	double _altAngle;
//	double _azimuthAngle;
//	double _viewDistance;
//	osg::Vec3d _viewDirection;
//	osg::Vec3d _centerEyePoint;
//	osg::Vec3d _center;
//	osg::Vec3d _upDirection;
//	double _left;
//	double _right;
//	double _bottom;
//	double _top;
//	double _znear;
//	double _zfar;
//	int _ncols;
//	int _nrows;
//	double _cellsize;
//	osg::BoundingBox _bb;
//	GDALDriver *_poDriver;
//	osg::Camera* _camera;
//	ObliqueCameraWrapper(osg::Camera* camera, const osg::BoundingBox& bb, double altAngle = 0, double azimuthAngle = 0)
//	{
//		_camera = camera;
//		float angle = osg::DegreesToRadians(altAngle);
//
//		osg::Vec3d bbmin(bb.xMin(), bb.yMin(), bb.zMin());
//		osg::Vec3d bbmax(bb.xMax(), bb.yMax(), bb.zMin());
//		osg::Vec3d bbsize = bbmax - bbmin;
//
//		float maxsize = bbsize.x() > bbsize.y() ? bbsize.x() : bbsize.y();
//		float viewDistance = maxsize * 5;
//
//		osg::Vec3d translation(0, -viewDistance * cos(angle), viewDistance * sin(angle));
//		if (altAngle == 0.1)
//		{
//			translation = osg::Vec3d(0, 0, viewDistance);
//		}
//		osg::Vec3d bbminNear = bbmin + translation;
//		osg::Vec3d bbmaxNear = bbmax + translation;
//		osg::Vec3d bbminFar = bbmin - translation;
//		osg::Vec3d bbmaxFar = bbmax - translation;
//		double znear = viewDistance - bb.radius();
//		double zfar = viewDistance + bb.radius();
//		float  right = (bb.xMax() - bb.xMin())*0.5;
//
//		_camera->setReferenceFrame(osg::Transform::ABSOLUTE_RF);
//
//		osg::Vec3d upDirection(0.0, 1.0, 0.0);
//		upDirection.normalize();
//		osg::Vec3d viewDirection = -translation;
//		viewDirection.normalize();
//		osg::Vec3d center = (bbmin + bbmax) * 0.5;
//
//		osg::Vec3d eyePoint = center - viewDirection * viewDistance;
//
//		_camera->setViewMatrixAsLookAt(eyePoint, center, upDirection);
//
//
//		float height = bb.zMax() - bb.zMin();
//		float part1 = height / tan(angle);
//		float part2 = (bb.yMax() - bb.yMin()) * 0.5;
//		float top = sin(angle) * (part1 + part2);
//		float bottom = sin(angle) * part2;
//
//		_camera->setProjectionMatrixAsOrtho(-right, right, -bottom, top, znear, zfar);
//
//		_azimuthAngle = osg::DegreesToRadians(azimuthAngle); ;
//		_altAngle = osg::DegreesToRadians(altAngle);
//		_viewDirection = viewDirection;
//		_center = center;
//		_left = -right;
//		_right = right;
//		_top = top;
//		_bottom = -bottom;
//		_znear = znear;
//		_zfar = zfar;
//		_viewDistance = viewDistance;
//		_upDirection = upDirection;
//
//		_bb = bb;
//		_center.y() = _center.y() + part1 * 0.5;
//
//		//return 0;
//		GDALAllRegister();
//		const char *pszFormat = "GTiff";
//
//		char **papszMetadata;
//
//		_poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
//		char **papszOptions = NULL;
//		//std::stringstream ssSrc;
//		//ssSrc << nrow << "_" << ncol << ".tif";
//
//	}
//
//	void setupGrid(double cellsize)
//	{
//		double width = _right - _left;
//		double height = _top - _bottom;
//
//
//
//
//		float margin = (_bb.zMax() - _bb.zMin()) / tan(_altAngle) * sin(_altAngle);
//
//		height += margin * 2;
//
//		_nrows = height / cellsize;
//		while (_nrows * cellsize < height)
//			_nrows = _nrows + 1;
//
//		_ncols = width / cellsize;
//		while (_ncols * cellsize < width)
//			_ncols = _ncols + 1;
//
//
//		_cellsize = cellsize;
//
//		//OGREnvelope newbound;
//		_left = _left - ((_ncols * cellsize) - width) * 0.5;
//		_right = _right + ((_ncols * cellsize) - width) * 0.5;
//		_bottom = _bottom - ((_nrows * cellsize) - height) * 0.5;
//		_top = _top + ((_nrows * cellsize) - height) * 0.5;
//
//		_bb.xMin() = _center.x() - width * 0.5;
//		_bb.xMax() = _center.x() + width * 0.5;
//		_bb.yMin() = _center.y() - height * 0.5;
//		_bb.yMax() = _center.y() + height * 0.5;
//
//
//	}
//
//
//	void setupCameraAtCell(int nrow, int ncol)
//	{
//		osg::Vec3d center = _center;
//		center.x() = _bb.xMin() + ncol * _cellsize + _cellsize * 0.5;
//		center.y() = _bb.yMax() - nrow * _cellsize - _cellsize * 0.5;
//		osg::Vec3d eyePoint = center - _viewDirection * _viewDistance;
//
//		//printf("%f,%f,%f: %f,%f,%f\n", center.x(), center.y(), center.z(), eyePoint.x(), eyePoint.y(), eyePoint.z());
//		double halfsize = _cellsize * 0.5;
//		float left = -halfsize;
//		float right = halfsize;
//		float bottom = -sin(_altAngle) * halfsize;
//		float top = sin(_altAngle) * halfsize;
//		float height = _bb.zMax() - _bb.zMin();
//		//float part1 = height / tan(_altAngle);
//		//float part2 = halfsize;
//		//if (nrow == 0)
//		//{
//		//	top = (part1 + part2) * sin(_altAngle);
//		//}
//		//if (nrow == _nrows - 1)
//		//{
//		//	bottom = -(part1 + part2) * sin(_altAngle);
//		//}
//
//		_camera->setViewMatrixAsLookAt(eyePoint, center, _upDirection);
//
//		_camera->setProjectionMatrixAsOrtho(left, right, bottom, top, _znear, _zfar);
//
//	}
//
//	void produceImage(int nrow, int ncol, osg::Image* img)
//	{
//
//
//		std::stringstream ssDest;
//		ssDest << "img/" << nrow << "_" << ncol << ".tif";
//		//		//g_pWarpIntance->copyImageMemory();
//		osgDB::writeImageFile(*img, ssDest.str().data());
//		GDALDataset *poDstDS = (GDALDataset *)GDALOpen(ssDest.str().data(), GA_Update);
//		double adfGeoTransform[6];
//		adfGeoTransform[0] = _bb.xMin() + ncol * _cellsize;///* top left x */
//		adfGeoTransform[1] = _cellsize / img->s();///* w-e pixel resolution */
//		adfGeoTransform[2] = 0;///* 0 */
//		adfGeoTransform[3] = _bb.yMax() - nrow * _cellsize;// /* top left y */
//		adfGeoTransform[4] = 0;///* 0 */
//		adfGeoTransform[5] = -_cellsize / img->t();///* n-s pixel resolution (negative value) */
//
//		poDstDS->SetGeoTransform(adfGeoTransform);
//		GDALClose((GDALDatasetH)poDstDS);
//
//	}
//
//};
void setupObliqueCamera(osg::Camera* camera, const osg::BoundingBox& bb, double altAngle = 0)
{
	altAngle = 45;
	float angle = osg::DegreesToRadians(altAngle);

	osg::Vec3d bbmin(bb.xMin(), bb.yMin(), bb.zMin());
	osg::Vec3d bbmax(bb.xMax(), bb.yMax(), bb.zMin());
	osg::Vec3d bbsize = bbmax - bbmin;

	float maxsize = bbsize.x() > bbsize.y() ? bbsize.x() : bbsize.y();
	float viewDistance = maxsize * 5;

	osg::Vec3d translation(0, -viewDistance * cos(angle), viewDistance * sin(angle));
	if (altAngle == 0.1)
	{
		translation = osg::Vec3d(0, 0, viewDistance);
	}
	osg::Vec3d bbminNear = bbmin + translation;
	osg::Vec3d bbmaxNear = bbmax + translation;
	osg::Vec3d bbminFar = bbmin - translation;
	osg::Vec3d bbmaxFar = bbmax - translation;
	double znear = viewDistance - bb.radius();
	double zfar = viewDistance + bb.radius();
	float  right= (bb.xMax() - bb.xMin())*0.5;

	camera->setReferenceFrame(osg::Transform::ABSOLUTE_RF);

	osg::Vec3d upDirection(0.0, 1.0, 0.0);
	upDirection.normalize();
	osg::Vec3d viewDirection = -translation;
	viewDirection.normalize();
	osg::Vec3d center = (bbmin + bbmax) * 0.5;

	osg::Vec3d eyePoint = center - viewDirection * viewDistance;

	camera->setViewMatrixAsLookAt(eyePoint, center, upDirection);


	float height = bb.zMax() - bb.zMin();
	float part1 = height / tan(angle);
	float part2 = (bb.yMax() - bb.yMin()) * 0.5;
	float halftop = sin(angle) * (part1 + part2);
	float halfbottom = sin(angle) * part2;

	camera->setProjectionMatrixAsOrtho(-right, right,-halfbottom, halftop,  znear, zfar);

	//for (size_t i = 0; i < 8; i++)
	//{
	//	osg::Vec3 pos = bb.corner(i) * camera->getViewMatrix() * camera->getProjectionMatrix();
	//	printfVec3(bb.corner(i), pos);
	//}

}
#include <osg/ShapeDrawable>
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


#include "osg/Point"
#include "osg/CullFace"
osg::Geode* createPoint(osg::Vec3d point)
{
	osg::Geode* pointgeode = new osg::Geode;
	osg::ref_ptr<osg::Program> program = new osg::Program;

	char vertexShaderSource[] =
		"void main(void)\n"
		"{\n"
		"   gl_Position    = gl_ModelViewProjectionMatrix *  gl_Vertex;\n"
		"}\n";
	char fragmentShaderSource[] =
		"void main(void) \n"
		"{\n"
		"     gl_FragColor = vec4(1,0,0,1);\n"
		"}\n";

	program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));
	program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexShaderSource));
	osg::StateSet* ss = pointgeode->getOrCreateStateSet();

	ss->setAttribute(program.get(), osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);

	pointgeode->setCullingActive(false);
	osg::ref_ptr<osg::Point> pointSize = new osg::Point(10);
	//point->setDistanceAttenuation(osg::Vec3(0.0, 0.0000, 0.05f));
	ss->setAttribute(pointSize.get());


	osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;
	osg::ref_ptr<osg::Geometry> polyGeom = new osg::Geometry();
	//point.z() +=0.01;
	vertices->push_back(point);
	
	polyGeom->setVertexArray(vertices.get());
	polyGeom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS, 0, vertices->size()));

	//osg::ref_ptr<osg::Vec3Array> normals = new osg::Vec3Array;
	//normals->push_back(osg::Vec3(0.0f, -1.0f, 0.0f));
	//pointsGeom->setNormalArray(normals);
	//pointsGeom->setNormalBinding(osg::Geometry::BIND_OVERALL);

	//osg::CullFace* cull = new osg::CullFace(osg::CullFace::BACK;
	pointgeode->addDrawable(polyGeom.get());
	polyGeom->getOrCreateStateSet()->setMode(GL_CULL_FACE, osg::StateAttribute::OFF| osg::StateAttribute::OVERRIDE);
	return pointgeode;
}




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

	//node->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
	//node->getOrCreateStateSet()->setMode( GL_BLEND, osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);
	osg::ref_ptr<osg::MatrixTransform> sceneNode = new osg::MatrixTransform;
	sceneNode->addChild(node.get());
	sceneNode->setMatrix(osg::Matrix::rotate(osg::DegreesToRadians(azimuthAngle), osg::Vec3(0,0,1)));
	osg::ComputeBoundsVisitor visitor;
	sceneNode->accept(visitor);
	osg::BoundingBox bb = visitor.getBoundingBox();

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
	viewer.getCamera()->setClearColor(osg::Vec4(0, 0, 0, 0));

	osg::ref_ptr<ObliqueCamera> orthoCamera = new ObliqueCamera(bb,altAngle);
	orthoCamera->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
	orthoCamera->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	orthoCamera->setClearColor(osg::Vec4(0, 0, 0, 0));
	
	//camera->setClearColor(osg::Vec4(0.53f, 0.85f, 1.0f, 0.9f));				// Background
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
	viewer.getCamera()->setClearColor(osg::Vec4(0, 0, 0, 0));
	//viewer.setCameraManipulator(new osgGA::TrackballManipulator);
	viewer.getCamera()->getOrCreateStateSet()->setMode(GL_BLEND, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
	//orthoCamera->getOrCreateStateSet()->setAttribute(clamp, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
	STPager* databasePager = new STPager;
	databasePager->cancel();
	viewer.getScene()->setDatabasePager(databasePager);

	//return 0;
	GDALAllRegister();
	const char *pszFormat = "GTiff";
	GDALDriver *poDriver;
	char **papszMetadata;

	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	char **papszOptions = NULL;
	double adfGeoTransform[6];
	orthoCamera->setupGrid(texSizeX*resol);
	

	ShapeFile shp("points.shp");
	OGRFeature *poFeature;
	shp.poLayer->ResetReading();
	int pointid = 0;
	osg::ref_ptr<osg::Node> pointNode;
	
	
	std::ofstream ofs( outdir + "georeference.txt");
	std::ofstream ofs2(outdir + "georeference2.txt");
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		/*if (pointid % 5 == 0)
		{*/
		OGRPoint* ogrp = (OGRPoint*)poFeature->GetGeometryRef();
		double x = poFeature->GetFieldAsDouble("x");
		double y = poFeature->GetFieldAsDouble("y");
		double z = poFeature->GetFieldAsDouble("z");
		double nx = poFeature->GetFieldAsDouble("nx");
		double ny = poFeature->GetFieldAsDouble("ny");
		double nz = poFeature->GetFieldAsDouble("nz");
		osg::Vec3d pos(x, y, z);
		osg::Vec3d normal(nx, ny, nz);
		pos = pos + normal * 0.01;
		osg::Vec3d pos2 = pos * sceneNode->getMatrix();
		int row;
		int col;
		orthoCamera->getRowCol(pos2, row, col);
		//if (row != 34 || col != 30)
		//	continue;
		//-3, 229.5
		bool found = false;
		if (pointNode && pointNode.valid())
			sceneNode->removeChild(pointNode);
		pointNode = createPoint(pos);
		sceneNode->addChild(pointNode.get());
		
	
		orthoCamera->setupCameraAtCell(row, col );
		viewer.frame();
		databasePager->frame();
		while (databasePager->getFileRequestListSize() > 0)
		{
			viewer.frame();
			databasePager->frame();
		}
		for (size_t n = 0; n < 10; n++)
		{
			viewer.frame();
		}

		osg::Vec3d destPos;
		if (row > -1 && col > -1)
		{
			if (orthoCamera->findPixelPosition(img, destPos))
			{
				ofs << destPos.x() << "	" << destPos.y() << "	" << pos.x() << "	" << pos.y() << std::endl;
				ofs2 << pointid << "	" << row << "	" << col << "	" << destPos.x() << "	" << destPos.y() << "	" << pos.x() << "	" << pos.y() << std::endl;
				//orthoCamera->produceImage(row, col, img, "./test/");
				//found = true;
				//orthoCamera->produceImage(row + i, col + j, img, "./test/");
				//break;
			}
		}
	
		OGRFeature::DestroyFeature(poFeature);
		pointid++;
		printf("%d\n", pointid);
		
	}
	ofs.close();
	ofs2.close();

	exit(0);
}





//

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
//	//QDir qoutdir("./img");
//	//DeleteDirectory(qoutdir.absolutePath());
//	//qoutdir.mkpath(".");
//	// add the screen capture handler
//	//viewer.addEventHandler(new osgViewer::ScreenCaptureHandler);
//	//osg::Node* node = osgDB::readNodeFile("F:/Oblique_Photogrammetry/weihai/Data/all.desc");
//	//osg::Node* node = osgDB::readNodeFile("E:/efficient-sparse-voxel-octrees-1.4/scenes/kechuang.ive");
//	//osg::Node* node = osgDB::readNodeFile("F:/Oblique_Photogrammetry/xinji/Data/all.desc");
//	osg::Node* node = osgDB::readNodeFile("F:/Oblique_Photogrammetry/weihai/Data/all.desc");
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
//
//	osg::ref_ptr<osg::PolygonMode> polymode = new osg::PolygonMode;
//	polymode->setMode(osg::PolygonMode::FRONT_AND_BACK, osg::PolygonMode::LINE);
//	//node->getOrCreateStateSet()->setAttributeAndModes(polymode, osg::StateAttribute::OVERRIDE | osg::StateAttribute::ON);
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
//		"     gl_FragColor = vec4(0,0,0,1);\n"
//		"}\n";
//
//	program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));
//	program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexShaderSource));
//
//	node->getOrCreateStateSet()->setAttribute(program.get(), osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
//
//	osg::MatrixTransform* sceneNode = new osg::MatrixTransform;
//	sceneNode->addChild(node);
//	sceneNode->setMatrix(osg::Matrix::rotate(osg::DegreesToRadians(90.0), osg::Vec3(0,0,1)));
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
//	int texSizeX = 1024;
//	int texSizeY = 1024;
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
//
//	osg::ref_ptr<osg::Image> img = new osg::Image;
//	//img->allocateImage(texSize, texSize,1,GL_ALPHA,GL_FLOAT);
//	img->allocateImage(texSizeX, texSizeY, 1, GL_RGBA, GL_UNSIGNED_BYTE);
//	//img->allocateImage(texSize, texSize,1,GL_RGBA,GL_FLOAT); 
//	tex->setImage(img.get());
//
//	//ObliqueCameraWrapper* orthoCamera = new ObliqueCameraWrapper(viewer.getCamera(),bb, 45);
//	viewer.getCamera()->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
//	viewer.getCamera()->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//	viewer.getCamera()->setClearColor(osg::Vec4(0, 0, 0, 1));
//
//	ObliqueCamera* orthoCamera = new ObliqueCamera(bb, 45);
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
//	viewer.setUpViewInWindow(0, 0, 1, 1);
//	//viewer.getCamera()->attach(osg::Camera::COLOR_BUFFER0, tex.get());
//	//viewer.getCamera()->attach(osg::Camera::COLOR_BUFFER0, img.get());
//	//orthoCamera->addChild(sceneNode);
//	//viewer.setSceneData(sceneNode);
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
//	OGRRegisterAll();
//
//	const char *pszFormat = "GTiff";
//	GDALDriver *poDriver;
//	char **papszMetadata;
//
//	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
//	char **papszOptions = NULL;
//	double adfGeoTransform[6];
//	orthoCamera->setupGrid(1024);
//	//orthoCamera->setupCameraAtCell(10, 10);
//	//viewer.frame();
//	//databasePager->frame();
//	//while (databasePager->getFileRequestListSize() > 0)
//	//{
//	//		viewer.frame();
//	//		databasePager->frame();
//	//}
//
//	//for (size_t i = 0; i < 20; i++)
//	//{
//	//	viewer.frame();
//	//}
//
//	//orthoCamera->produceImage(10, 10, img,"B:/VGE/DSMRenderer_2016_02_14/x64/Release/test/");
//	ShapeFile shp("points.shp");
//	OGRFeature *poFeature;
//	shp.poLayer->ResetReading();
//	int pointid = 0;
//	osg::ref_ptr<osg::Node> pointNode;
//
//
//	//void transformGCPs(const std::string infile, const std::string outfile, const osg::Matrix& mat)
//	//{
//	//	std::ifstream ifs(infile.data());
//
//	//	std::ofstream ofs(outfile.data());
//	//	char buf[1024];
//	//	while (ifs.getline(buf, 1024))
//	//	{
//	//		std::stringstream str(buf);
//	//		double x1, y1, x2, y2;
//	//		str >> x1 >> y1 >> x2 >> y2;
//	//		osg::Vec3d pos(x1, y1, 0);
//	//		pos = pos * mat;
//	//		x1 = pos.x();
//	//		y1 = pos.y();
//	//		ofs << x1 << "	" << y1 << "	" << x2 << "	" << y2 << std::endl;
//
//	//	}
//	//	ifs.close();
//	//	ofs.flush();
//	//	ofs.close();
//	//}
//
//	std::ofstream ofs("georeference_45_90.txt");
//	std::ofstream ofs2("georeference_45_90_2.txt");
//	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
//	{
//		/*if (pointid % 5 == 0)
//		{*/
//		OGRPoint* ogrp = (OGRPoint*)poFeature->GetGeometryRef();
//		double x = poFeature->GetFieldAsDouble("x");
//		double y = poFeature->GetFieldAsDouble("y");
//		double z = poFeature->GetFieldAsDouble("z");
//		double nx = poFeature->GetFieldAsDouble("nx");
//		double ny = poFeature->GetFieldAsDouble("ny");
//		double nz = poFeature->GetFieldAsDouble("nz");
//		osg::Vec3d pos(x, y, z);
//		osg::Vec3d normal(nx, ny, nz);
//		pos = pos + normal * 0.01;
//		osg::Vec3d pos2 = pos * sceneNode->getMatrix();
//		int row;
//		int col;
//		orthoCamera->getRowCol(pos2, row, col);
//		//if (row != 34 || col != 30)
//		//	continue;
//		//-3, 229.5
//		bool found = false;
//		if (pointNode && pointNode.valid())
//			sceneNode->removeChild(pointNode);
//		pointNode = createPoint(pos);
//		sceneNode->addChild(pointNode.get());
//		//for (int i = -1; i <= 1; i++)
//		//{
//		//	if (found)
//		//		break;
//		//	if (row + i < 0 || row + i > orthoCamera->_nrows - 1)
//		//		continue;
//		//	if (found)
//		//		break;
//		//	for (int j = -1; j <= 1; j++)
//		//	{
//		//		if (col + j < 0 || col + j > orthoCamera->_ncols - 1)
//		//			continue;
//		//		orthoCamera->setupCameraAtCell(row + i, col + j);
//		//		viewer.frame();
//		//		databasePager->frame();
//		//		while (databasePager->getFileRequestListSize() > 0)
//		//		{
//		//			viewer.frame();
//		//			databasePager->frame();
//		//		}
//		//		for (size_t n = 0; n < 5; n++)
//		//		{
//		//			viewer.frame();
//		//		}
//
//
//
//		//		osg::Vec3d destPos;
//
//		//		if (orthoCamera->findPixelPosition(img, destPos))
//		//		{
//		//			ofs << destPos.x() << "	" << destPos.y() << "	" << pos.x() << "	" << pos.y() << std::endl;
//		//			found = true;
//		//			//orthoCamera->produceImage(row + i, col + j, img, "./test/");
//		//			break;
//		//		}
//
//		//		//viewer.frame();
//		//		//std::stringstream ss;
//		//		//ss << pointid << ".tif";
//		//		//ss << orthoCamera->_curRow << "_" << orthoCamera->_curCol << ".tif";
//
//
//		//	}
//		//}
//
//				orthoCamera->setupCameraAtCell(row, col );
//				viewer.frame();
//				databasePager->frame();
//				while (databasePager->getFileRequestListSize() > 0)
//				{
//					viewer.frame();
//					databasePager->frame();
//				}
//				for (size_t n = 0; n < 10; n++)
//				{
//					viewer.frame();
//				}
//
//
//
//				osg::Vec3d destPos;
//				if (row > -1 && col > -1)
//				{
//					if (orthoCamera->findPixelPosition(img, destPos))
//					{
//						ofs << destPos.x() << "	" << destPos.y() << "	" << pos.x() << "	" << pos.y() << std::endl;
//						ofs2 << pointid << "	" << row << "	" << col << "	" << destPos.x() << "	" << destPos.y() << "	" << pos.x() << "	" << pos.y() << std::endl;
//						//orthoCamera->produceImage(row, col, img, "./test/");
//						//found = true;
//						//orthoCamera->produceImage(row + i, col + j, img, "./test/");
//						//break;
//					}
//				}
//
//				
//		
//
//		OGRFeature::DestroyFeature(poFeature);
//		pointid++;
//		printf("%d\n", pointid);
//	
//	}
//	ofs.close();
//	ofs2.close();
//	exit(0);
//}