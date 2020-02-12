
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
#include <osg/ShapeDrawable>
#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include "osg/LineWidth"
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

void printfVec3(osg::Vec3 p1, osg::Vec3 p2)
{
	printf("(%f,%f,%f),(%f,%f,%f)\n", p1.x(), p1.y(), p1.z(), p2.x(), p2.y(), p2.z());
}

#include "osg/LineSegment"

class ObliqueCamera : public osg::Camera, public osgGA::GUIEventHandler
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
	double _srcTileSizeY;
	double _srcTileSizeX;
	osg::BoundingBox _bb;
	osg::BoundingBox _oribb;

	GDALDriver *_poDriver;
	int _curRow;
	int _curCol;
	osg::Group* _childFrustrumNode;
	ObliqueCamera(const osg::BoundingBox& bb, double altAngle = 0, double azimuthAngle = 0)
	{
		_altAngle = osg::DegreesToRadians(altAngle);

		osg::Vec3d bbmin(bb.xMin(), bb.yMin(), bb.zMin());
		osg::Vec3d bbmax(bb.xMax(), bb.yMax(), bb.zMin());
		osg::Vec3d bbsize = bbmax - bbmin;

		float maxsize = bbsize.x() > bbsize.y() ? bbsize.x() : bbsize.y();
		_viewDistance = maxsize * 5;

		osg::Vec3d translation(0, -_viewDistance * cos(_altAngle), _viewDistance * sin(_altAngle));
		if (altAngle == 0.1)
		{
			translation = osg::Vec3d(0, 0, _viewDistance);
		}
		osg::Vec3d bbminNear = bbmin + translation;
		osg::Vec3d bbmaxNear = bbmax + translation;
		osg::Vec3d bbminFar = bbmin - translation;
		osg::Vec3d bbmaxFar = bbmax - translation;
		_znear = _viewDistance - bb.radius();
		_zfar = _viewDistance + bb.radius();
		float  right = (bb.xMax() - bb.xMin())*0.5;

		this->setReferenceFrame(osg::Transform::ABSOLUTE_RF);

		_upDirection = osg::Vec3d(0.0, 1.0, 0.0);
		_upDirection.normalize();
		_viewDirection = -translation;
		_viewDirection.normalize();
		_center = (bbmin + bbmax) * 0.5;

		osg::Vec3d eyePoint = _center - _viewDirection * _viewDistance;
		this->setViewMatrixAsLookAt(eyePoint, _center, _upDirection);

		_bb = bb;
		_oribb = bb;
		GDALAllRegister();
		const char *pszFormat = "GTiff";

		char **papszMetadata;

		_poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
		char **papszOptions = NULL;
		//std::stringstream ssSrc;
		//ssSrc << nrow << "_" << ncol << ".tif";
	}

	bool handle(const osgGA::GUIEventAdapter& ea, osgGA::GUIActionAdapter& aa)
	{
		if (ea.getEventType() == osgGA::GUIEventAdapter::RELEASE && ea.getButton() == osgGA::GUIEventAdapter::RIGHT_MOUSE_BUTTON)
		{

		}
		else if (ea.getKey() == osgGA::GUIEventAdapter::KEY_B  &&  ea.getEventType() == osgGA::GUIEventAdapter::KEYUP)
		{
			_childFrustrumNode->setNodeMask(!_childFrustrumNode->getNodeMask());
		}

		return false;
	}
	void bindShader(osg::StateSet* stateset,osg::Vec4 color)
	{
		osg::ref_ptr<osg::Program> program = new osg::Program;

		/*	char vertexShaderSource[] =
		"void main(void)\n"
		"{\n"
		"   gl_Position    = gl_ModelViewProjectionMatrix *  gl_Vertex;\n"
		"}\n"*/;
		char fragmentShaderSource[] =
			"uniform vec4 color; \n"
			"void main(void) \n"
			"{\n"
			"     gl_FragColor = color;\n"
			"}\n";

		program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));
		//program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexShaderSource));
		stateset->getOrCreateUniform("color", osg::Uniform::FLOAT_VEC4)->set(color);
		stateset->setAttribute(program.get(), osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
	}
	osg::Geometry* createLineBox(osg::Vec3d center,double xsizehalf,double ysizehalf,double zsizehalf)
	{
		//osg::Geode* geode = new osg::Geode;
		osg::Geometry* geom = new osg::Geometry;
		osg::ref_ptr<osg::Vec3dArray> vertices = new osg::Vec3dArray;
		vertices->push_back(osg::Vec3(center.x() - xsizehalf, center.y() - ysizehalf, center.z() - zsizehalf));
		vertices->push_back(osg::Vec3(center.x() - xsizehalf, center.y() + ysizehalf, center.z() - zsizehalf));
		vertices->push_back(osg::Vec3(center.x() - xsizehalf, center.y() + ysizehalf, center.z() - zsizehalf));
		vertices->push_back(osg::Vec3(center.x() + xsizehalf, center.y() + ysizehalf, center.z() - zsizehalf));
		vertices->push_back(osg::Vec3(center.x() + xsizehalf, center.y() + ysizehalf, center.z() - zsizehalf));
		vertices->push_back(osg::Vec3(center.x() + xsizehalf, center.y() - ysizehalf, center.z() - zsizehalf));
		vertices->push_back(osg::Vec3(center.x() + xsizehalf, center.y() - ysizehalf, center.z() - zsizehalf));
		vertices->push_back(osg::Vec3(center.x() - xsizehalf, center.y() - ysizehalf, center.z() - zsizehalf));


		vertices->push_back(osg::Vec3(center.x() - xsizehalf, center.y() - ysizehalf, center.z() + zsizehalf));
		vertices->push_back(osg::Vec3(center.x() - xsizehalf, center.y() + ysizehalf, center.z() + zsizehalf));
		vertices->push_back(osg::Vec3(center.x() - xsizehalf, center.y() + ysizehalf, center.z() + zsizehalf));
		vertices->push_back(osg::Vec3(center.x() + xsizehalf, center.y() + ysizehalf, center.z() + zsizehalf));
		vertices->push_back(osg::Vec3(center.x() + xsizehalf, center.y() + ysizehalf, center.z() + zsizehalf));
		vertices->push_back(osg::Vec3(center.x() + xsizehalf, center.y() - ysizehalf, center.z() + zsizehalf));
		vertices->push_back(osg::Vec3(center.x() + xsizehalf, center.y() - ysizehalf, center.z() + zsizehalf));
		vertices->push_back(osg::Vec3(center.x() - xsizehalf, center.y() - ysizehalf, center.z() + zsizehalf));


		vertices->push_back(osg::Vec3(center.x() - xsizehalf, center.y() - ysizehalf, center.z() - zsizehalf));
		vertices->push_back(osg::Vec3(center.x() - xsizehalf, center.y() - ysizehalf, center.z() + zsizehalf));
		vertices->push_back(osg::Vec3(center.x() - xsizehalf, center.y() + ysizehalf, center.z() - zsizehalf));
		vertices->push_back(osg::Vec3(center.x() - xsizehalf, center.y() + ysizehalf, center.z() + zsizehalf));
		vertices->push_back(osg::Vec3(center.x() + xsizehalf, center.y() + ysizehalf, center.z() - zsizehalf));
		vertices->push_back(osg::Vec3(center.x() + xsizehalf, center.y() + ysizehalf, center.z() + zsizehalf));
		vertices->push_back(osg::Vec3(center.x() + xsizehalf, center.y() - ysizehalf, center.z() - zsizehalf));
		vertices->push_back(osg::Vec3(center.x() + xsizehalf, center.y() - ysizehalf, center.z() + zsizehalf));
		geom->setVertexArray(vertices.get());
		geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES, 0, vertices->size()));
		//geode->addChild(geom.get());
		return geom;
	}

	osg::Node* createLine(osg::Vec3d center, double translate,osg::Vec3d dir, double len, int width,osg::Vec4 color)
	{
		osg::ref_ptr<osg::Geode> geode = new osg::Geode;
		osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
		osg::ref_ptr<osg::Vec3dArray> vertices = new osg::Vec3dArray;
		vertices->push_back(osg::Vec3(0, 0, -len * 0.5 + translate));
		vertices->push_back(osg::Vec3(0, 0, len * 0.5 + translate));

		geom->setVertexArray(vertices.get());
		geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES, 0, vertices->size()));
		
		bindShader(geom->getOrCreateStateSet(), color);

		osg::LineWidth* linewidth = new osg::LineWidth();
		linewidth->setWidth(width);
		geom->getOrCreateStateSet()->setAttributeAndModes(linewidth, osg::StateAttribute::OVERRIDE | osg::StateAttribute::ON);
		geom->getOrCreateStateSet()->setMode(GL_CULL_FACE, osg::StateAttribute::OFF);
		geode->addDrawable(geom.get());

		osg::MatrixTransform* transform = new osg::MatrixTransform;

		transform->setMatrix(osg::Matrix::rotate(osg::Vec3d(0, 0, 1), dir) * osg::Matrix::translate(center));

		transform->addChild(geode.get());
		geode->getOrCreateStateSet()->setMode(GL_BLEND, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
		return transform;
	}

	osg::Node* createFrustrumGrid(osg::Vec3d& center, double translate,osg::Vec3d& viewDirection, int nrows, int ncols, double left2right, double top2bottom, double near2far, osg::Vec4 linecolor)
	{
		
		osg::ref_ptr<osg::Geode> geode = new osg::Geode();


		osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
		osg::ref_ptr<osg::Vec3dArray> vertices = new osg::Vec3dArray;
		double cellsizeY = top2bottom / nrows;
		double cellsizeX = left2right / ncols;
		double minx =  - left2right * 0.5;
		double miny =  - top2bottom * 0.5;

		for (size_t i = 1; i < nrows; i++)
		{
			vertices->push_back(osg::Vec3(minx + 0, miny + i*cellsizeY, translate));
			vertices->push_back(osg::Vec3(minx + left2right, miny + i*cellsizeY, translate));
		}
		for (size_t i = 1; i < ncols; i++)
		{
			vertices->push_back(osg::Vec3(minx + i*cellsizeX, miny + 0, translate));
			vertices->push_back(osg::Vec3(minx + i*cellsizeX, miny + top2bottom, translate));
		}

		//for (size_t i = 1; i < nrows; i++)
		//{
		//	vertices->push_back(osg::Vec3(minx + 0, miny + i*cellsizeY, center.z() + near2far*0.5));
		//	vertices->push_back(osg::Vec3(minx + left2right, miny + i*cellsizeY, center.z() + near2far*0.5));
		//}
		//for (size_t i = 1; i < ncols; i++)
		//{
		//	vertices->push_back(osg::Vec3(minx + i*cellsizeX, miny + 0, center.z() + near2far*0.5));
		//	vertices->push_back(osg::Vec3(minx + i*cellsizeX, miny + top2bottom, center.z() + near2far*0.5));
		//}
		geom->setVertexArray(vertices.get());
		geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES, 0, vertices->size()));



		osg::LineWidth* linewidth = new osg::LineWidth();
		linewidth->setWidth(2.0f);
		geom->getOrCreateStateSet()->setAttributeAndModes(linewidth, osg::StateAttribute::OVERRIDE | osg::StateAttribute::ON);
		geom->getOrCreateStateSet()->setMode(GL_CULL_FACE, osg::StateAttribute::OFF);
		geom->getOrCreateStateSet()->setMode(GL_BLEND, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);

		geode->addChild(geom.get());

		osg::MatrixTransform* transform = new osg::MatrixTransform;

		transform->setMatrix(osg::Matrix::rotate(osg::Vec3d(0, 0, 1), viewDirection) * osg::Matrix::translate(center.x(),center.y(),center.z()));

		transform->addChild(geode.get());

		return transform;

	}

	osg::Node* createFrustrumGeo()
	{
		osg::Vec3d center = _center;
		double left2right = _bb.xMax() - _bb.xMin();
		double top2bottom = (_bb.yMax() - _bb.yMin()) * sin(_altAngle);
		double near2far = top2bottom * 3;
		osg::Vec3d viewDirection = _viewDirection;
		osg::Group* root = new osg::Group;

		osg::Vec3d up(0, 0, 1);

		osg::Vec3d bbmin(_oribb.xMin(), _oribb.yMin(), _oribb.zMin());
		osg::Vec3d bbmax(_oribb.xMax(), _oribb.yMax(), _oribb.zMax());
		osg::Vec3d bbcenter = (bbmin+ bbmax) * 0.5;
		double height = _oribb.zMax() - _oribb.zMin();
		double axislen = near2far * 5;

		double cellsizeY = top2bottom / _nrows;
		double cellsizeX = left2right / _ncols;
		_childFrustrumNode = new osg::Group;
		_childFrustrumNode->addChild(createFrustrumGrid(center, near2far* 0.5, _viewDirection, _nrows,_ncols, left2right, top2bottom, near2far, osg::Vec4(0, 0, 0, 1)));
		_childFrustrumNode->addChild(createFrustrumGrid(center, -near2far* 0.5, _viewDirection, _nrows, _ncols, left2right, top2bottom, near2far, osg::Vec4(0, 0, 0, 1)));
		_childFrustrumNode->addChild(createFrustrumGeo2(center, osg::Vec3d(0, -cellsizeY*0.5, 0), _viewDirection, cellsizeX, cellsizeY, near2far+10, osg::Vec4(1, 0, 0, 0.5), osg::Vec4(0, 0, 0, 1)));
		_childFrustrumNode->setNodeMask(false);
		//root->addChild(_childFrustrumNode);

		//root->addChild(createLine(bbcenter, -near2far * 0.5,_viewDirection, near2far*2, 2, osg::Vec4(0, 0, 0, 1.0)));
		root->addChild(createLine(bbcenter, axislen * 0.5, osg::Vec3d(1, 0, 0), axislen, 2, osg::Vec4(1, 0, 0, 1.0)));
		root->addChild(createLine(bbcenter, axislen * 0.5, osg::Vec3d(0, 1, 0), axislen, 2, osg::Vec4(0, 1, 0, 1.0)));
		root->addChild(createLine(bbcenter, axislen * 0.5, osg::Vec3d(0, 0, 1), axislen, 2, osg::Vec4(0, 0, 1, 1.0)));
		root->addChild(createFrustrumGeo2(bbcenter, osg::Vec3d(0, 0, 0), up, _oribb.xMax() - _oribb.xMin(), _oribb.yMax() - _oribb.yMin(), height, osg::Vec4(0, 1, 0, 0.25), osg::Vec4(0, 0, 0, 1)));
		//root->addChild(createFrustrumGeo2(center, osg::Vec3d(0, 0, 0), _viewDirection, left2right, top2bottom, near2far, osg::Vec4(1, 1, 0, 0.5), osg::Vec4(0, 0, 0, 1)));
		//root->getOrCreateStateSet()->setMode(GL_BLEND, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
		//root->getOrCreateStateSet()->setRenderingHint(osg::StateSet::TRANSPARENT_BIN);
		return root;

	}
	osg::Node* createFrustrumGeo2(osg::Vec3d& center, osg::Vec3d translate, osg::Vec3d& viewDirection, double left2right, double top2bottom, double near2far,osg::Vec4 bbcolor, osg::Vec4 linecolor)
	{
		/*c
		osg::Vec3d viewDirection = _viewDirection;
		double left2right = _bb.xMax() - _bb.xMin();
		double top2bottom = _bb.yMax() - _bb.yMin();
		double near2far = _viewDistance;*/
		osg::ref_ptr<osg::Geode> geode = new osg::Geode();
		//osg::StateSet* stateset = new osg::StateSet();
		//geode->setStateSet(stateset);
		osg::TessellationHints* hints = new osg::TessellationHints;
		hints->setDetailRatio(0.5f);

		osg::ref_ptr<osg::ShapeDrawable> boxPolygon = new osg::ShapeDrawable(new osg::Box(osg::Vec3(0.0f, 0.0f, 0.0f), left2right, top2bottom, near2far));
		
		osg::ref_ptr<osg::PolygonMode> polymode = new osg::PolygonMode;
		polymode->setMode(osg::PolygonMode::FRONT, osg::PolygonMode::FILL);
		boxPolygon->getOrCreateStateSet()->setAttributeAndModes(polymode.get(), osg::StateAttribute::OVERRIDE | osg::StateAttribute::ON);
		boxPolygon->getOrCreateStateSet()->setMode(GL_BLEND, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
		boxPolygon->getOrCreateStateSet()->setRenderingHint(osg::StateSet::TRANSPARENT_BIN);
		bindShader(boxPolygon->getOrCreateStateSet(), bbcolor);
		osg::ref_ptr<osg::Geometry> boxPolyline = createLineBox(osg::Vec3d(0,0,0), left2right*0.5, top2bottom*0.5, near2far*0.5);
		bindShader(boxPolyline->getOrCreateStateSet(), linecolor);


		osg::LineWidth* linewidth = new osg::LineWidth();
		linewidth->setWidth(3.0f);
		boxPolyline->getOrCreateStateSet()->setAttributeAndModes(linewidth,osg::StateAttribute::OVERRIDE | osg::StateAttribute::ON);
		boxPolyline->getOrCreateStateSet()->setMode(GL_CULL_FACE, osg::StateAttribute::OFF);
		boxPolyline->getOrCreateStateSet()->setMode(GL_BLEND, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
		
		
		geode->addDrawable(boxPolyline.get());
		geode->addDrawable(boxPolygon.get());
		//osg::ref_ptr<osg::ShapeDrawable> boxLine = new osg::ShapeDrawable(new osg::Box(osg::Vec3(0.0f, 0.0f, 0.0f), left2right, top2bottom, near2far));
		//geode->addDrawable(boxLine.get());

		//osg::ref_ptr<osg::PolygonMode> polymode = new osg::PolygonMode;
		//polymode->setMode(osg::PolygonMode::FRONT_AND_BACK, osg::PolygonMode::LINE);
		//boxLine->getOrCreateStateSet()->setAttributeAndModes(polymode.get(), osg::StateAttribute::OVERRIDE | osg::StateAttribute::ON);
		//bindShader(boxLine->getOrCreateStateSet(), osg::Vec4(1, 0, 0, 1));

		//stateset->setMode(GL_LIGHTING, osg::StateAttribute::ON);
		osg::MatrixTransform* transform = new osg::MatrixTransform;

		transform->setMatrix(osg::Matrix::translate(translate) * osg::Matrix::rotate(osg::Vec3d(0, 0, 1), viewDirection) * osg::Matrix::translate(center));

		transform->addChild(geode.get());

		return transform;

	}





	void setupGrid(double destTileSize)
	{
		_srcTileSizeY = destTileSize * tan(_altAngle);
		_srcTileSizeX = destTileSize;
		double zspan = _bb.zMax() - _bb.zMin();
		double xspan = _bb.xMax() - _bb.xMin();
		double yspan = _bb.yMax() - _bb.yMin();

		double yspanProjected = yspan + zspan / tan(_altAngle);


		_nrows = (int)(yspanProjected / _srcTileSizeY);
		while (_nrows * _srcTileSizeY < yspanProjected)
			_nrows = _nrows + 1;

		_ncols = (int)(xspan / _srcTileSizeX);
		while (_ncols * _srcTileSizeX < xspan)
			_ncols = _ncols + 1;
		_center.y() = _center.y() + zspan / tan(_altAngle) * 0.5;
		xspan = _ncols * _srcTileSizeX;
		osg::BoundingBox bb = _bb;
		bb.xMin() = _center.x() - _ncols * _srcTileSizeX * 0.5;
		bb.xMax() = _center.x() + _ncols * _srcTileSizeX * 0.5;
		bb.yMin() = _center.y() - _nrows * _srcTileSizeY * 0.5;
		bb.yMax() = _center.y() + _nrows * _srcTileSizeY * 0.5;

		_bb = bb;
	}

	void setupGrid(double destTileSizeX, double destTileSizeY)
	{
		_srcTileSizeY = destTileSizeY * tan(_altAngle);
		_srcTileSizeX = destTileSizeX;
		double zspan = _bb.zMax() - _bb.zMin();
		double xspan = _bb.xMax() - _bb.xMin();
		double yspan = _bb.yMax() - _bb.yMin();

		double yspanProjected = yspan + zspan / tan(_altAngle);


		_nrows = (int)(yspanProjected / _srcTileSizeY);
		while (_nrows * _srcTileSizeY < yspanProjected)
			_nrows = _nrows + 1;

		_ncols = (int)(xspan / _srcTileSizeX);
		while (_ncols * _srcTileSizeX < xspan)
			_ncols = _ncols + 1;


		//_nrows = _ncols = 1;
		_center.y() = _center.y() + zspan / tan(_altAngle) * 0.5;
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
		osg::Vec3d center = _center;
		center.x() = _bb.xMin() + ncol * _srcTileSizeX + _srcTileSizeX * 0.5;
		center.y() = _bb.yMax() - nrow * _srcTileSizeY - _srcTileSizeY * 0.5;
		osg::Vec3d eyePoint = center - _viewDirection * _viewDistance;
		float left = -_srcTileSizeX * 0.5;
		float right = _srcTileSizeX * 0.5;

		float bottom = -_srcTileSizeY * sin(_altAngle) * 0.5;
		float top = _srcTileSizeY * sin(_altAngle) * 0.5;
		this->setViewMatrixAsLookAt(eyePoint, center, _upDirection);

		this->setProjectionMatrixAsOrtho(left, right, bottom, top, _znear, _zfar);
		_curRow = nrow;
		_curCol = ncol;
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


		//GDALDataset *poDstDS = (GDALDataset *)GDALOpen(ssDest.str().data(), GA_Update);
		double adfGeoTransform[6];
		adfGeoTransform[0] = _bb.xMin() + ncol * _srcTileSizeX;///* top left x */
		adfGeoTransform[1] = _srcTileSizeX / img->s();///* w-e pixel resolution */
		adfGeoTransform[2] = 0;///* 0 */
		adfGeoTransform[3] = _bb.yMax() - nrow * _srcTileSizeY / tan(_altAngle);// /* top left y */
		adfGeoTransform[4] = 0;///* 0 */
		adfGeoTransform[5] = -_srcTileSizeX / img->t();///* n-s pixel resolution (negative value) */



		char **papszOptions = NULL;

		papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");
		papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "LZW");
		GDALDataset * poDstDS = _poDriver->CreateCopy(ssDest.str().data(), poSrcDS, FALSE,
			papszOptions, GDALTermProgress, NULL);

		poDstDS->SetGeoTransform(adfGeoTransform);

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

		//GDALDataset *poDstDS = (GDALDataset *)GDALOpen(ssDest.str().data(), GA_Update);
		double adfGeoTransform[6];
		adfGeoTransform[0] = _bb.xMin() + _curCol * _srcTileSizeX;///* top left x */
		adfGeoTransform[1] = _srcTileSizeX / img->s();///* w-e pixel resolution */
		adfGeoTransform[2] = 0;///* 0 */
		adfGeoTransform[3] = _bb.yMax() - _curRow * _srcTileSizeY / tan(_altAngle);// /* top left y */
		adfGeoTransform[4] = 0;///* 0 */
		adfGeoTransform[5] = -_srcTileSizeX / img->t();///* n-s pixel resolution (negative value) */


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
		double adfGeoTransform[6];
		double halfsize = _srcTileSizeY * 0.5;
		
		adfGeoTransform[0] = _bb.xMin() + _curCol * _srcTileSizeY;///* top left x */
		adfGeoTransform[1] = _srcTileSizeY / img->s();///* w-e pixel resolution */
		adfGeoTransform[2] = 0;///* 0 */
		adfGeoTransform[3] = _bb.yMax() - _curRow * _srcTileSizeY;// /* top left y */
		adfGeoTransform[4] = 0;///* 0 */
		adfGeoTransform[5] = -_srcTileSizeY / img->t();///* n-s pixel resolution (negative value) */
		std::vector<osg::Vec3d> pixels;
		unsigned char* data = img->data();
		for (size_t i = 0; i < img->t(); i++)
		{
			double y = adfGeoTransform[3] + adfGeoTransform[5] * (img->t() - i - 1) + adfGeoTransform[5] * 0.5;
			for (size_t j = 0; j < img->s(); j++)
			{
				double x = adfGeoTransform[0] + adfGeoTransform[1] * j + adfGeoTransform[1] * 0.5;
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

void lookat(osg::Vec3d eye, osg::Vec3d center,osg::Vec3d up)
{
	//Vec3d f(center - eye);
	//f.normalize();
	//Vec3d s(f^up);
	//s.normalize();
	//Vec3d u(s^f);
	//u.normalize();

	osg::Vec3d z(center - eye);
	z.normalize();
	osg::Vec3d y(z^up);
	y.normalize();
	osg::Vec3d x(y^z);
	x.normalize();
	osg::Vec3d translate(x*eye,y*eye,z*eye);
	printf("%f,%f,%f", translate.x(), translate.y(), translate.z());
	//osg::Vec3d v;
	//double mat[4][4];
	//mat[0] = {};
	////	s[0], u[0], -f[0], 0.0,
	////	s[1], u[1], -f[1], 0.0,
	////	s[2], u[2], -f[2], 0.0,
	////	0.0, 0.0, 0.0, 1.0);

	//for (unsigned i = 0; i < 3; ++i)
	//{
	//	double tmp = v[i];
	//	mat[3][0] += tmp*mat[i][0];
	//	mat[3][1] += tmp*mat[i][1];
	//	mat[3][2] += tmp*mat[i][2];
	//	mat[3][3] += tmp*mat[i][3];
	//}
}



//
//set(
//	s[0], u[0], -f[0], 0.0,
//	s[1], u[1], -f[1], 0.0,
//	s[2], u[2], -f[2], 0.0,
//	0.0, 0.0, 0.0, 1.0);
//
//preMultTranslate(-eye);



int main(int argc, char** argv)
{
	//osg::Matrix viewMat = osg::Matrix::lookAt(osg::Vec3(100, 0, 0), osg::Vec3(0, 0, 0), osg::Vec3(0, 0, 1));
	//lookat(osg::Vec3(100, 0, 100), osg::Vec3(0, 0, 0), osg::Vec3(0, 0, 1));
	//double len = osg::Vec3(100, 100, 100).length();
	//osg::Vec3 point = osg::Vec3(0, 1000, 0) * viewMat;

	osg::Node* scene = osgDB::readNodeFile("F:/Work/VGE/Jiangshan/JSmoxing/wangou.osg");
	//scene->getOrCreateStateSet()->setMode(GL_CULL_FACE,osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
	scene->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
	
	//osg::ComputeBoundsVisitor cbv;
	//scene->accept(cbv);
	//osg::MatrixTransform* matnode = new osg::MatrixTransform;
	//osg::Vec3d center = -cbv.getBoundingBox().center();
	//center.z() = -cbv.getBoundingBox().center().z() + (cbv.getBoundingBox().zMax() - cbv.getBoundingBox().zMin()) * 0.5;
	//matnode->setMatrix(osg::Matrix::translate(center));
	//matnode->addChild(scene);
	//osgDB::writeNodeFile(*matnode, "F:/Work/VGE/Jiangshan/JSmoxing/wangou.ive");
	osg::Vec3d right(1, 0, 0);
	osg::Vec3d forward(0, 1, 0);
	osg::Vec3d up = forward ^ right ;
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
	node->getOrCreateStateSet()->setMode(GL_BLEND, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
	node->getOrCreateStateSet()->setRenderingHint(osg::StateSet::OPAQUE_BIN);
	osg::ref_ptr<osg::MatrixTransform> sceneNode = new osg::MatrixTransform;
	sceneNode->addChild(node.get());
	osg::Matrix mat = osg::Matrix::rotate(osg::DegreesToRadians(azimuthAngle), osg::Vec3(0, 0, 1));
	sceneNode->setMatrix(mat);
	osg::ComputeBoundsVisitor visitor;
	sceneNode->accept(visitor);
	osg::BoundingBox bb = visitor.getBoundingBox();
	
 //   osg::MatrixTransform* trans = new osg::MatrixTransform;
	//osg::ComputeBoundsVisitor cbv;
	//osg::Vec3 bbdim = osg::Vec3(bb.xMax() - bb.xMin(), bb.yMax() - bb.yMin(), bb.zMax() - bb.zMin());
	//osg::Vec3 scale = osg::Vec3(100/bbdim.x(), 100 / bbdim.y(), 100 / bbdim.z());
	//trans->addChild(node);
	//trans->setMatrix(osg::Matrix::scale(scale) * osg::Matrix::rotate(osg::DegreesToRadians(90.0), osg::Vec3(0, 0, 1)));
	//osgDB::writeNodeFile(*trans, "F:/Work/VGE/Jiangshan/JSmoxing/building2.ive");
	//	osg::MatrixTransform* trans = new osg::MatrixTransform;
	//osg::ComputeBoundsVisitor cbv;
	//trans->addChild(root);
	//osg::Vec3 bbdim = osg::Vec3(bb.xMax() - bb.xMin(), bb.yMax() - bb.yMin(), bb.zMax() - bb.zMin());
	//osg::Vec3 scale = osg::Vec3(100/bbdim.x(), 100 / bbdim.y(), 100 / bbdim.z());
	//trans->setMatrix(osg::Matrix::scale(scale) * osg::Matrix::rotate(osg::DegreesToRadians(90.0), osg::Vec3(0, 0, 1)));
	//trans->accept(cbv);
	//osg::BoundingBox newbb = cbv.getBoundingBox();


	//osgDB::writeNodeFile(*trans, "F:/Work/VGE/Jiangshan/JSmoxing/building.ive");

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

		//bb = newbb;
	}
	//int texSizeX = (int)((bb.xMax() - bb.xMin()) / resol);
	//double zspan = bb.zMax() - bb.zMin();
	//double projectedLengthY = (bb.yMax() - bb.yMin()) + zspan / tan(osg::DegreesToRadians(altAngle));
	//int texSizeY = (int)(projectedLengthY  / resol);
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
	//viewer.getCamera()->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
	viewer.getCamera()->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	viewer.getCamera()->setClearColor(osg::Vec4(1, 1, 1, 0));

	ObliqueCamera* orthoCamera = new ObliqueCamera(bb,altAngle);
	orthoCamera->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
	orthoCamera->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	orthoCamera->setClearColor(osg::Vec4(1, 0, 0, 1));
	
	orthoCamera->setReferenceFrame(osg::Transform::ABSOLUTE_RF_INHERIT_VIEWPOINT);
	orthoCamera->setViewport(0, 0, texSizeX, texSizeY);
	orthoCamera->setRenderOrder(osg::Camera::PRE_RENDER);
	orthoCamera->setRenderTargetImplementation(osg::Camera::FRAME_BUFFER_OBJECT);
	orthoCamera->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	viewer.setUpViewInWindow(0, 0, texSizeX, texSizeY);
	orthoCamera->attach(osg::Camera::COLOR_BUFFER0, tex.get());
	orthoCamera->attach(osg::Camera::COLOR_BUFFER0, img.get());
	orthoCamera->addChild(sceneNode.get());
	//viewer.setSceneData(orthoCamera.get());
	orthoCamera->setupGrid(texSizeX*resol, texSizeY*resol);
	viewer.addEventHandler((osgGA::GUIEventHandler*)orthoCamera);
	osg::ref_ptr<osg::Group> root = new osg::Group;
	root->addChild(sceneNode.get());
	root->addChild(orthoCamera->createFrustrumGeo());
	//osgDB::writeNodeFile(*root, "F:/Work/VGE/Jiangshan/JSmoxing/building.ive");
	viewer.setSceneData(root.get());
	viewer.setUpViewInWindow(100, 100, 1024, 768);
	viewer.setCameraManipulator(new osgGA::TrackballManipulator);
	viewer.realize();
	//viewer.getCamera()->getOrCreateStateSet()->setMode(GL_BLEND, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
	//STPager* databasePager = new STPager;
	//databasePager->cancel();
	//viewer.getScene()->setDatabasePager(databasePager);

	return viewer.run();



	//////return 0;
	//GDALAllRegister();
	//std::vector<std::string> tilefiles;
	////char **papszOptions = NULL;
	////double adfGeoTransform[6];
	//for (int irow = 0; irow < orthoCamera->_nrows; irow++)
	//{
	//	for (int icol = 0; icol < orthoCamera->_ncols; icol++)
	//	{
	//		//if (orthoCamera->imagExists(irow, icol, img, outdir))
	//		//	continue;
	//		orthoCamera->setupCameraAtCell(irow, icol);
	//		viewer.frame();
	//		databasePager->frame();
	//		while (databasePager->getFileRequestListSize() > 0)
	//		{
	//			viewer.frame();
	//			databasePager->frame();
	//		}
	//		viewer.frame();
	//		tilefiles.push_back(orthoCamera->produceImage(irow, icol, img,outdir));
	//	}
	//}

	//const char *pszFormat = "GTiff";
	//GDALDriver *poDriver;
	//poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	//std::stringstream commandss;
	//commandss << "gdalbuildvrt " << outdir + "oblique.vrt " << outdir + "*.tif";
	//system(commandss.str().data());


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

