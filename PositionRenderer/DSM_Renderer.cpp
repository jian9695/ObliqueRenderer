
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
//替换+字符，原因是含有该字符的文件无法在天启服务器实现发布

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
int main(int argc, char** argv)
{
	osg::DisplaySettings::instance()->setNumOfDatabaseThreadsHint(1);
	osg::DisplaySettings::instance()->setNumOfHttpDatabaseThreadsHint(0);


	osgViewer::Viewer viewer;
	//viewer.getScene()->setDatabasePager()
	// add the state manipulator
	viewer.addEventHandler( new osgGA::StateSetManipulator(viewer.getCamera()->getOrCreateStateSet()) );

	// add the thread model handler
	viewer.addEventHandler(new osgViewer::ThreadingHandler);

	// add the window size toggle handler
	viewer.addEventHandler(new osgViewer::WindowSizeHandler);

	// add the stats handler
	viewer.addEventHandler(new osgViewer::StatsHandler);


	// add the record camera path handler
	//viewer.addEventHandler(new osgViewer::RecordCameraPathHandler);

	// add the LOD Scale handler
	//viewer.addEventHandler(new osgViewer::LODScaleHandler);

	// add the screen capture handler
	//viewer.addEventHandler(new osgViewer::ScreenCaptureHandler);
	osg::Node* node = osgDB::readNodeFile("E:/efficient-sparse-voxel-octrees-1.4/scenes/kechuang.ive");
	std::stringstream ss_outdir;
	std::string outdir = "E:/ClimateSkyModel64/Kechuang/";


	//osg::Node* sceneNode = osgDB::readNodeFile("D:/Projects/Oblique_Photogrammetry/suzhou_gaoxin/Data/all.desc");
	//osg::Node* node = osgDB::readNodeFile("D:/Projects/Oblique_Photogrammetry/xinji/Data/all.desc");
	node->getOrCreateStateSet()->setMode(GL_LIGHTING,osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
	//node->getOrCreateStateSet()->setMode( GL_BLEND, osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);
	osg::MatrixTransform* sceneNode = new osg::MatrixTransform;
	sceneNode->addChild(node);
	/*osg::Matrix rot;
	double obliqueAngle = 90;
	double heading = 255.5;
	double pitch   = 1;
	double roll    = 359;
	rot.makeRotate( 
		osg::DegreesToRadians(heading), osg::Vec3(0,0,1),
		osg::DegreesToRadians(pitch),   osg::Vec3(1,0,0),
		osg::DegreesToRadians(roll),    osg::Vec3(0,1,0) );
	osg::Matrix mat = osg::Matrix::scale(24,24,24);
	mat = mat * rot;
	osg::Matrix obliqueMat = osg::Matrix::rotate(osg::DegreesToRadians(obliqueAngle), osg::Vec3(0,0,1));
	mat = mat * obliqueMat;*/
	//sceneNode->setMatrix(mat);
	/*<model name="嘉善城区模型">
		<range range_max="30000" range_min="1" range_on="true" />
		<position hat="271" lat="30.85100056359863" long="120.9086136235352" srs="+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs " />
		<scale x="24" y="24" z="24" />
		<local_rotation w="-0.6122172800344491" x="0" y="0" z="0.7906895737438435" />
		<style name="modelNode">
		<symbols>
		<model heading="255.5" pitch="1" roll="359" url="../../model/Data/all.desc" />
		</symbols>
		</style>
		</model>*/
	//osg::Node* sceneNode = osgDB::readNodeFile("cow.osg");

	//创建一个相机节点
	osg::ComputeBoundsVisitor visitor;
	visitor.apply(*sceneNode);
	osg::BoundingBox bb = visitor.getBoundingBox();
	osg::BoundingBox tileBB;
	int texSize = 512;
	double resolution = 0.2;
	float bbW = bb.xMax() - bb.xMin();
	float bbH = bb.yMax() - bb.yMin();

	float tileW = texSize * 0.25;
	float tileH = texSize * 0.25;

	float curX = bb.xMin();
	float curY = bb.yMax();
	int nrows = (int)(bbH / tileH)+1;
	int ncols = (int)(bbW / tileW)+1;
	//osg::Camera* camera2 = createScreenQuadCamera(node,bb,rtTexture);
	//camera2->addChild(screenQuad);

	//orthoCamera->getOrCreateStateSet()->setRenderBinDetails(2,"RenderBin"); 
	osg::ref_ptr<osg::Texture2D> tex = new osg::Texture2D;

	tex->setTextureSize(texSize, texSize);
	tex->setResizeNonPowerOfTwoHint(false);
	tex->setFilter(osg::Texture2D::MIN_FILTER,osg::Texture2D::LINEAR);
	tex->setFilter(osg::Texture2D::MAG_FILTER,osg::Texture2D::LINEAR);
	tex->setWrap(osg::Texture2D::WRAP_S,osg::Texture2D::REPEAT);
	tex->setWrap(osg::Texture2D::WRAP_T,osg::Texture2D::REPEAT);
	//rtTexture->setDataVariance(osg::Object::DYNAMIC);
	//tex->setInternalFormat(GL_RGBA);
	//tex->setSourceFormat(GL_RGBA);
	//tex->setSourceType(GL_UNSIGNED_BYTE);
	tex->setInternalFormat(GL_ALPHA32F_ARB);
	tex->setSourceFormat(GL_ALPHA);
	tex->setSourceType(GL_FLOAT);
	//tex->setInternalFormat(GL_RGBA32F_ARB);
	//tex->setSourceFormat(GL_RGBA);
	//tex->setSourceType(GL_FLOAT);
	osg::ref_ptr<osg::Image> img = new osg::Image;
	img->allocateImage(texSize, texSize,1,GL_ALPHA,GL_FLOAT);
    //img->allocateImage(texSize, texSize,1,GL_RGBA,GL_UNSIGNED_BYTE); 
    //img->allocateImage(texSize, texSize,1,GL_RGBA,GL_FLOAT); 
	tex->setImage(img.get());
	osg::StateSet* stateset = node->getOrCreateStateSet();
	osg::ref_ptr<osg::Program> program = new osg::Program;
	char vertexShaderSource[] = 
		"varying vec4 pos;\n"
		"void main(void)\n"
		"{\n"
		"gl_TexCoord[0] = gl_MultiTexCoord0;\n"
		"pos = gl_Vertex;\n"
		"gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\n"
		"}\n";
	char fragmentShaderSource[] = 

		"uniform vec4 color;\n"
		"uniform sampler2D tex;\n"
		"uniform sampler2D texMask;\n"
		"varying vec4 pos;\n"
		"void main(void) \n"
		"{\n"
		"//gl_FragColor = texture2D(tex, gl_TexCoord[0].xy);\n"
		"gl_FragData[0] = vec4(pos.z);\n"
		"}\n";

	program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));
	program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexShaderSource));
	/*osg::ref_ptr<osg::Image> maskImg = osgEarth::URI(baseDir+"mask.png").getImage();
	if(!maskImg || !maskImg.valid())
	{
		maskImg = new osg::Image;
		maskImg->allocateImage(1,1,1,GL_RGBA,GL_UNSIGNED_BYTE);
		unsigned char* data = maskImg->data();
		data[0]=255;data[1]=255;data[2]=255;data[3]=255;

	}*/
	node->accept(visitor);
	osg::BoundingBox bound = visitor.getBoundingBox();
	//构建地形图
	osg::ref_ptr<osg::Texture2D> maskTex = new osg::Texture2D;  
	maskTex->setName( "maskTex" );
	maskTex->setFilter(osg::Texture::MIN_FILTER, osg::Texture::LINEAR_MIPMAP_LINEAR);//_MIPMAP_LINEAR
	maskTex->setFilter(osg::Texture::MAG_FILTER, osg::Texture::LINEAR);
	//maskTex->setImage(maskImg.get());
	//stateset->addUniform( new osg::Uniform("osg_ViewMatrixInverse", osg::Matrix::identity()) );
	stateset->addUniform( new osg::Uniform("tex", 0) );
	//stateset->addUniform( new osg::Uniform("texPolygon", 1) );
	//stateset->addUniform( new osg::Uniform("texMask", 2) );
	//stateset->addUniform( new osg::Uniform("bound", osg::Vec4(bound.xMin(),bound.xMax(),bound.yMin(),bound.yMax())) );
	//stateset->setTextureAttributeAndModes( 2, maskTex.get(), osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);
	//stateset->setTextureAttributeAndModes( 1, maskTex.get(), osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE );
	//stateset->setTextureAttributeAndModes( 0, maskTex.get(), osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE );
	//stateset->setAttribute(program.get(),osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE); 
	//mat->setCullingActive(false);
	stateset->setAttribute(program.get(),osg::StateAttribute::ON |osg::StateAttribute::OVERRIDE);  
	osg::ClampColor* clamp = new osg::ClampColor();
	clamp->setClampVertexColor(GL_FALSE);
	clamp->setClampFragmentColor(GL_FALSE);
	clamp->setClampReadColor(GL_FALSE);
	stateset->setAttribute(clamp, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
	osg::Camera* orthoCamera = new osg::Camera;
	orthoCamera->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
	orthoCamera->setClearMask( GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT );
	orthoCamera->setClearColor( osg::Vec4(0,0,0,0));

	//camera->setClearColor(osg::Vec4(0.53f, 0.85f, 1.0f, 0.9f));				// Background
	orthoCamera->setReferenceFrame( osg::Transform::ABSOLUTE_RF_INHERIT_VIEWPOINT );
	orthoCamera->setViewport( 0,0, texSize, texSize);
	orthoCamera->setRenderOrder(osg::Camera::PRE_RENDER);
	orthoCamera->setRenderTargetImplementation(osg::Camera::FRAME_BUFFER_OBJECT); 
	orthoCamera->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	viewer.setUpViewInWindow(0,0,texSize,texSize);
	orthoCamera->attach(osg::Camera::COLOR_BUFFER0,tex.get());
	orthoCamera->attach(osg::Camera::COLOR_BUFFER0,img.get());
	orthoCamera->addChild(sceneNode);
	viewer.setSceneData(orthoCamera);
	viewer.realize();
	viewer.getCamera()->setClearColor( osg::Vec4(0,0,0,0));
	//viewer.setCameraManipulator(new osgGA::TrackballManipulator);
	viewer.getCamera()->getOrCreateStateSet()->setMode( GL_BLEND, osg::StateAttribute::OFF|osg::StateAttribute::OVERRIDE);
	orthoCamera->getOrCreateStateSet()->setAttribute(clamp, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
	STPager* databasePager = new STPager;
	databasePager->cancel();
	viewer.getScene()->setDatabasePager(databasePager);

	//transformGCPs("../jiashan/georeference.txt",outdir + "georeference.txt",obliqueMat);
	//return 0;
	GDALAllRegister();
	const char *pszFormat = "GTiff";
	GDALDriver *poDriver;
	char **papszMetadata;

	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	char **papszOptions = NULL;
	double adfGeoTransform[6];
	osg::Vec3d translate(0,0,0);
	//float* tmpBuf = new float[texSize*texSize];
	//for(int i=0;i<texSize*texSize;i++)
	//{
	//	tmpBuf[i] = 1000;
	//}
	for (int irow=0;irow<nrows;irow++)
	{
		for (int icol=0;icol<ncols;icol++)
		{
			curX = bb.xMin() + icol * tileW;
			curY = bb.yMax() - irow * tileH;
			tileBB = osg::BoundingBox(curX,curY-tileH,bb.zMin(),curX+tileW,curY,bb.zMax());
			updateCameraMatrixOrtho(viewer.getCamera(),tileBB,90);
			updateCameraMatrixOrtho(orthoCamera,tileBB,90);
            //viewer.frame();
			//for (int n=0;n<50;n++)
			//{
			//	databasePager->frame();
			//    viewer.frame();
			//} 
			viewer.frame();
			databasePager->frame();

			//for (int n=0;n<50;n++)
			//{
			//	databasePager->frame();
			//    viewer.frame();
			//} 
			while(databasePager->getFileRequestListSize() > 0)
			{
				viewer.frame();
				databasePager->frame();
			}
	   		viewer.frame();     

	/*		std::stringstream ss;;
			ss << outdir.data() << irow << "_" << icol << ".png";*/
			//osgDB::writeImageFile(*img,ss.str());



			std::stringstream ssDest;
			ssDest << outdir << irow << "_" << icol << ".tif";
			//g_pWarpIntance->copyImageMemory();
			//osgDB::writeImageFile(*img,ss.str());
		/*	GDALDataset *poSrcDS = 
				(GDALDataset *) GDALOpen(ssSrc.str().data(), GA_ReadOnly );
			if(!poSrcDS)
				continue;*/

			adfGeoTransform[0] = curX+translate.x();///* top left x */
			adfGeoTransform[1] = (double)tileW/(double)texSize;///* w-e pixel resolution */
			adfGeoTransform[2] = 0;///* 0 */
			adfGeoTransform[3] = curY+translate.y();// /* top left y */
			adfGeoTransform[4] = 0;///* 0 */
			adfGeoTransform[5] =-(double)tileH/(double)texSize;///* n-s pixel resolution (negative value) */
			GDALDataset *poDstDS = poDriver->Create( ssDest.str().data(),texSize,texSize,1,GDT_Float32, NULL );
			poDstDS->SetGeoTransform(adfGeoTransform);
			//poDstDS->SetProjection( pszSRS_WKT );
			//double adfGeoTransform[6] = { 444720, 30, 0, 3751320, 0, -30 };
			osg::ref_ptr<osg::Image> newImg = new osg::Image(*img,osg::CopyOp::DEEP_COPY_ALL);
			newImg->flipVertical();
			float* buf = (float*)newImg->data();
			/*		for (int k=0;k<texSize*texSize;k++)
			{
			if(buf[k] > 10)
			printf("%f\n",buf[k]);
			}*/
			poDstDS->GetRasterBand(1)->SetNoDataValue(0);
			poDstDS->RasterIO(GF_Write,0, 0,texSize, texSize, 
				buf, texSize, texSize, GDT_Float32,1,NULL,sizeof(float) , sizeof(float) * texSize,0);
			////if( poDstDS != NULL )
			poDstDS->FlushCache();
			GDALClose( (GDALDatasetH) poDstDS );
			//GDALClose( (GDALDatasetH) poSrcDS );
			printf("%d/%d\n",irow,nrows);
            //img->flipVertical();

			//printf("%d/%d\n",irow,nrows);
		}
	}

	//OGRSpatialReference oSRS;
	//char *pszSRS_WKT = NULL;
	//oSRS.SetUTM( 50, TRUE );
	//oSRS.SetWellKnownGeogCS( "WGS84" );
	//oSRS.exportToWkt( &pszSRS_WKT );

	//for (int irow=0;irow<nrows;irow++)
	//{
	//	for (int icol=0;icol<ncols;icol++)
	//	{
	//		curX = bb.xMin() + icol * tileW;
	//		curY = bb.yMax() - irow * tileH;

	//		std::stringstream ssSrc;
	//		ssSrc << outdir << irow << "_" << icol << ".jpg";
	//		std::stringstream ssDest;
	//		ssDest << outdir << irow << "_" << icol << ".tif";
	//		//g_pWarpIntance->copyImageMemory();
	//		//osgDB::writeImageFile(*g_pWarpedTexture->getImage(),ss.str());
	//		GDALDataset *poSrcDS = 
	//			(GDALDataset *) GDALOpen(ssSrc.str().data(), GA_ReadOnly );
	//		if(!poSrcDS)
	//			continue;

	//		adfGeoTransform[0] = curX+translate.x();///* top left x */
	//		adfGeoTransform[1] = (double)tileW/(double)texSize;///* w-e pixel resolution */
	//		adfGeoTransform[2] = 0;///* 0 */
	//		adfGeoTransform[3] = curY+translate.y();// /* top left y */
	//		adfGeoTransform[4] = 0;///* 0 */
	//		adfGeoTransform[5] =-(double)tileH/(double)texSize;///* n-s pixel resolution (negative value) */
	//		GDALDataset *poDstDS = poDriver->CreateCopy( ssDest.str().data(), poSrcDS, FALSE, 
	//			NULL, NULL, NULL );
	//		poDstDS->SetGeoTransform(adfGeoTransform);
	//		//poDstDS->SetProjection( pszSRS_WKT );
	//		//double adfGeoTransform[6] = { 444720, 30, 0, 3751320, 0, -30 };
	//		if( poDstDS != NULL )
	//			GDALClose( (GDALDatasetH) poDstDS );
	//		GDALClose( (GDALDatasetH) poSrcDS );
	//		printf("%d/%d\n",irow,nrows);
	//	}
	//}
	/*while(!viewer.done())
	{
		updateCameraMatrixOrtho(viewer.getCamera(),bb);
		databasePager->frame();
		viewer.frame();
	}*/
	return 1;

}
