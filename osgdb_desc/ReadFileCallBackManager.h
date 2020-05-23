#pragma once
//#include  "vge_foundation_global.h"
#include <Windows.h>
#include <map>
#include "CallBackGroup.h"
#include "osgDB\ReaderWriter"
#include "osgDB\Callbacks"
#include "osgDB\Registry"
#include "osgDB\Options"
using namespace osgDB;
using namespace osg;

namespace VGE
{

	class OnReadFileCallBack : public CallBackObject
	{
	public:

		virtual void beforeOpenArchive(ReaderWriter::ReadResult& rr,const std::string& filename,ReaderWriter::ArchiveStatus status, unsigned int indexBlockSizeHint, const Options* useObjectCache){}

		virtual void beforeReadObject(ReaderWriter::ReadResult& rr,const std::string& filename, const Options* options){}

		virtual void beforeReadImage(ReaderWriter::ReadResult& rr,const std::string& filename, const Options* options){}

		virtual void beforeReadHeightField(ReaderWriter::ReadResult& rr,const std::string& filename, const Options* options){}

		virtual void beforeReadNode(ReaderWriter::ReadResult& rr,const std::string& filename, const Options* options){}

		virtual void beforeReadShader(ReaderWriter::ReadResult& rr,const std::string& filename, const Options* options){}

		virtual void afterOpenArchive(ReaderWriter::ReadResult& rr,const std::string& filename,ReaderWriter::ArchiveStatus status, unsigned int indexBlockSizeHint, const Options* useObjectCache){}

		virtual void afterReadObject(ReaderWriter::ReadResult& rr,const std::string& filename, const Options* options){}

		virtual void afterReadImage(ReaderWriter::ReadResult& rr,const std::string& filename, const Options* options){}

		virtual void afterReadHeightField(ReaderWriter::ReadResult& rr,const std::string& filename, const Options* options){}

		virtual void afterReadNode(ReaderWriter::ReadResult& rr,const std::string& filename, const Options* options){}

		virtual void afterReadShader(ReaderWriter::ReadResult& rr,const std::string& filename, const Options* options){}

	};


	class /*VGE_FOUNDATION_EXPORT*/ ReadFileCallBackManager : public CallBackGroup,public ReadFileCallback
	{
	private:
		static ReadFileCallBackManager* g_pInstance;
		ReadFileCallBackManager(void){}
	public:
		static ReadFileCallBackManager* instance();

		virtual ReaderWriter::ReadResult openArchive(const std::string& filename,ReaderWriter::ArchiveStatus status, unsigned int indexBlockSizeHint, const Options* useObjectCache);

		virtual ReaderWriter::ReadResult readObject(const std::string& filename, const Options* options);

		virtual ReaderWriter::ReadResult readImage(const std::string& filename, const Options* options);

		virtual ReaderWriter::ReadResult readHeightField(const std::string& filename, const Options* options);

		virtual ReaderWriter::ReadResult readNode(const std::string& filename, const Options* options);

		virtual ReaderWriter::ReadResult readShader(const std::string& filename, const Options* options);

	};

}

