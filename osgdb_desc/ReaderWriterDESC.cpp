/* -*-c++-*- OpenSceneGraph - Copyright (C) 1998-2006 Robert Osfield
 *
 * This application is open source and may be redistributed and/or modified
 * freely and without restriction, both in commercial and non commercial
 * applications, as long as this copyright notice is maintained.
 *
 * This application is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
*/

#include <sstream>
#include <memory>

#include <osg/Notify>
#include <osg/NodeVisitor>
#include <osgDB/ReaderWriter>
#include <osgDB/FileNameUtils>
#include <osgDB/Registry>
#include <osgDB/ConvertUTF>

#include <OpenThreads/ScopedLock>

#include "ReaderWriterDESC.h"
#include <osgDB/ReadFile>
#include <osg/MatrixTransform>
#ifdef WIN32
#include "windows.h"
#endif
#include "osg/MatrixTransform"
#include <osgDB/FileUtils>
#include "osg/ComputeBoundsVisitor"
#include "osgEarth/URI"
#include "osg/Texture2D"
#include "osg/PagedLOD"
#include "ReadFileCallBackManager.h"
#include "qfileinfo.h"
#include "GeometryVisitor.h"
#include "osgDB\ReaderWriter"
#define SERIALIZER() OpenThreads::ScopedLock<OpenThreads::ReentrantMutex> lock(_serializerMutex)
using namespace osg;
using namespace osgDB;
using namespace VGE;


std::vector<std::string> ReaderWriterDESC::split(const char& delimiter,const std::string& line) const
{
	std::string str = line;
	std::vector<std::string> splits;
	if(str.length() < 1)
		return splits;
	std::string::iterator iter = str.begin();
	int index = 0;
	int lastIndex = 0;
	while(iter != str.end() && index < str.length())
	{
		char& c = *iter;
		if(c ==delimiter || index == str.length()-1)
		{
			if(index-lastIndex>0)
			   splits.push_back(str.substr(lastIndex,index-lastIndex));
			lastIndex = index+1;
		}
		index++;
		iter++;
	}
	return splits;
}
std::string ReaderWriterDESC::getDirectory(const std::string& filename) const
{
	std::string  tmpStr = filename;
	 std::string::iterator iter = tmpStr.end();
	 iter--;
	 int lastIndex = -1;
	 int index = tmpStr.length()-1;
	 while(iter != tmpStr.begin())
	 {
		 char& c = *iter;
		 if(c =='/' || c == '\\')
		 {
			 lastIndex = index;
			 break;
		 }
		 index--;
		 iter--;
	 }
	std::string dir = "";
	if(lastIndex > -1)
	{
		dir = tmpStr.substr(0,lastIndex+1);
	}
	return dir;
}
std::string ReaderWriterDESC::replace(const std::string& name) const
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
			else if(*it == '\\')
			{
				*it = '/';
			}
		}
		return str;
	}
	
	return name;
}

void bindShader(osg::Node* node,std::string baseDir)
{
	osg::BoundingBox bound;
	osg::ComputeBoundsVisitor visitor;
	node->accept(visitor);
	osg::BoundingBox bb = visitor.getBoundingBox();
	bound.expandBy(bb);
	osg::StateSet* stateset = node->getOrCreateStateSet();
	stateset->setMode(GL_NORMALIZE, osg::StateAttribute::OVERRIDE | osg::StateAttribute::ON);
	stateset->setMode(GL_LIGHTING, osg::StateAttribute::OVERRIDE | osg::StateAttribute::OFF);
	osg::ref_ptr<osg::Image> maskImg = osgEarth::URI(baseDir+"mask.png").getImage();
	if(!maskImg || !maskImg.valid())
	{
		maskImg = new osg::Image;
		maskImg->allocateImage(1,1,1,GL_RGBA,GL_UNSIGNED_BYTE);
		unsigned char* data = maskImg->data();
		data[0]=255;data[1]=255;data[2]=255;data[3]=255;

	}

	osg::ref_ptr<osg::Program> program = new osg::Program;
	char vertexShaderSource[] = 
		"uniform vec4 bound;\n"
		"uniform vec4 bound2;\n"
		"uniform mat4 osg_ViewMatrixInverse;\n"
		"uniform sampler2D texDSM;\n"
		"uniform bool hasDSM;\n"
		"void main(void)\n"
		"{\n"
		"gl_TexCoord[0] = gl_MultiTexCoord0;\n"
		"vec4 pos = gl_Vertex;\n"
		"vec4 min = vec4(bound.x,bound.z,0,1);\n"
		"vec4 max = vec4(bound.y,bound.w,0,1);\n"
		"gl_TexCoord[1] = vec4((pos.x-min.x)/(max.x-min.x),(pos.y-min.y)/(max.y-min.y),0,1);\n"
		"min = vec4(bound2.x,bound2.z,0,1);\n"
		"max = vec4(bound2.y,bound2.w,0,1);\n"
		"gl_TexCoord[2] = vec4((pos.x-min.x)/(max.x-min.x),(pos.y-min.y)/(max.y-min.y),0,1);\n"
		"float h = texture2D(texDSM, gl_TexCoord[2].xy).a;\n"
		"//if(h > -1000 && h < 9000)\n"
		"  //pos.z = h;\n"
		"gl_Position = gl_ModelViewProjectionMatrix * pos;\n"
		"}\n";
	char fragmentShaderSource[] = 
		"uniform vec4 color;\n"
		"uniform sampler2D tex;\n"
		"uniform sampler2D texPolygon;\n"
		"uniform sampler2D texMask;\n"
		"uniform sampler2D texDOM;\n"
		"uniform sampler2D texDSM;\n"
		"uniform bool hasDSM;\n"
		"void main(void) \n"
		"{\n"
		"vec4 maskColor = texture2D(texMask, gl_TexCoord[1].xy);\n"
		"if(maskColor.r < 0.1)\n"
		"{\n"
		"   discard;\n"
		"}\n"
		"vec4 texColor = texture2D(tex, gl_TexCoord[0].xy);\n"
		"if(texture2D(texDSM, gl_TexCoord[2].xy).r > 0.5)\n"
		"   texColor = texture2D(texDOM, gl_TexCoord[2].xy);\n"
		"gl_FragColor = texColor;\n"
		"}\n";

	program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));
	program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexShaderSource));

	//构建地形图
	osg::ref_ptr<osg::Texture2D> maskTex = new osg::Texture2D;  
//	maskTex->setName( "maskTex" );
	maskTex->setFilter(osg::Texture::MIN_FILTER, osg::Texture::LINEAR_MIPMAP_LINEAR);//_MIPMAP_LINEAR
	maskTex->setFilter(osg::Texture::MAG_FILTER, osg::Texture::LINEAR);
	maskTex->setImage(maskImg.get());
	//stateset->addUniform( new osg::Uniform("osg_ViewMatrixInverse", osg::Matrix::identity()) );
	if(!stateset->getUniform("tex"))
		stateset->addUniform( new osg::Uniform("tex", 0) );
	if(!stateset->getUniform("texPolygon"))
		stateset->addUniform( new osg::Uniform("texPolygon", 1) );
	if(!stateset->getUniform("texMask"))
		stateset->addUniform( new osg::Uniform("texMask", 2) );
	if(!stateset->getUniform("texDOM"))
		stateset->addUniform( new osg::Uniform("texDOM", 3) );
	if(!stateset->getUniform("texDSM"))
		stateset->addUniform( new osg::Uniform("texDSM", 4) );
	if(!stateset->getUniform("bound"))
		stateset->addUniform( new osg::Uniform("bound", osg::Vec4(bound.xMin(),bound.xMax(),bound.yMin(),bound.yMax())) );
	if(!stateset->getUniform("bound2"))
		stateset->addUniform( new osg::Uniform("bound2", osg::Vec4(bound.xMin(),bound.xMax(),bound.yMin(),bound.yMax())) );
	if(!stateset->getUniform("hasDSM"))
		stateset->addUniform( new osg::Uniform("hasDSM", false) );
	//stateset->setTextureAttributeAndModes( 2, maskTex.get(), osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);
	//stateset->setTextureAttributeAndModes( 1, maskTex.get(), osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE );
	//stateset->setTextureAttributeAndModes( 0, maskTex.get(), osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE );
	//stateset->setAttribute(program.get(),osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE); 
	//osg::ClampColor* clamp = new osg::ClampColor();
	//clamp->setClampVertexColor(GL_FALSE);
	//clamp->setClampFragmentColor(GL_FALSE);
	//clamp->setClampReadColor(GL_FALSE);
	//stateset->setAttribute(clamp, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
}

osgDB::ReaderWriter::ReadResult
	ReaderWriterDESC::readNodeFile(const std::vector<std::string>& lines,const std::string& dir) const
{
	/* SERIALIZER();*/

	std::string baseDir = dir;
	if(baseDir[baseDir.length()-1] != '/' && baseDir[baseDir.length()-1] != '\\')
	{
		baseDir += "/";
	}
    osg::Group* mat = new osg::Group;
	osgDB::ReaderWriter* reader = osgDB::Registry::instance()->getReaderWriterForExtension("osgb");
	bool isHttp = false;
	if(dir.length() > 5 && dir.substr(0,4) == "http")
	{
		isHttp = true;
	}
	osg::BoundingBox bound;
	for (int i=0;i<lines.size();i++)
	{    
		std::string line = lines[i];
		if(line.length() < 5)
			continue;
		std::vector<std::string> splits = split(',',line);
		if(splits.size()>1)
		{
			line  = splits[0];
		}
		std::string fullname = baseDir+line;//
		if(isHttp)
		{
			fullname = replace(fullname);
			if(fullname[fullname.length()-1] != 'c' && fullname[fullname.length()-1] != 'e' && fullname[fullname.length()-1] != 'b')
			   fullname = fullname.substr(0,fullname.length()-1);
		}
	
		/*if(splits.size()>5)
		{

		osg::BoundingBox bb(atof(splits[1].data()),atof(splits[2].data()),atof(splits[3].data()),
		atof(splits[4].data()),atof(splits[5].data()),atof(splits[6].data()));
		osg::ref_ptr<osg::PagedLOD> pageLOD = new osg::PagedLOD();
		pageLOD->setCenter(bb.center());
		pageLOD->setRadius(bb.radius());
		pageLOD->setFileName(0,fullname);
		pageLOD->setCenterMode(osg::LOD::USER_DEFINED_CENTER);
		pageLOD->setRangeMode(osg::LOD::PIXEL_SIZE_ON_SCREEN);
		pageLOD->setRange(0,10,100000000);
		mat->addChild(pageLOD.get());
		bound.expandBy(bb);
		}
		else
		{*/
			osg::ref_ptr<osg::Node> node = osgEarth::URI(fullname).getNode();
			if(!node|| !node.valid())
				continue;
			mat->addChild(node.get());
			osg::ComputeBoundsVisitor visitor;
			mat->traverse(visitor);
			osg::BoundingBox bb = visitor.getBoundingBox();
			bound.expandBy(bb);
		//}

	}
	
	//bindShader(mat,dir);
	//static ReadModelCallBack* readModelCallBack = NULL;
	//if(!readModelCallBack)
	//{
	//	readModelCallBack = new ReadModelCallBack;
	//    ReadFileCallBackManager::instance()->addCallBack(readModelCallBack);
	//}

	return mat;
}
osgDB::ReaderWriter::ReadResult
	ReaderWriterDESC::readNode(std::istream& fin,
	const osgDB::ReaderWriter::Options* options) const
{
	/* SERIALIZER();*/
	if(options->getDatabasePathList().size() < 1)
		return ReadResult::FILE_NOT_FOUND;
	std::string dir = options->getDatabasePathList()[0];
	if(dir.empty())
		return ReadResult::FILE_NOT_FOUND;
	char str[1024];  
	std::vector<std::string> lines;
	while( fin.getline(str,1024) )
	{    
		std::string line = str;
		if(line.length() < 5)
			continue;
		lines.push_back(str);
	}
	return readNodeFile(lines, dir);
}
osgDB::ReaderWriter::ReadResult
ReaderWriterDESC::readNode(const std::string& fname,
        const osgDB::ReaderWriter::Options* options) const
{
   /* SERIALIZER();*/

	std::string ext = osgDB::getLowerCaseFileExtension(fname);
	if (!acceptsExtension(ext)) return ReadResult::FILE_NOT_HANDLED;
	std::string filename = osgDB::findDataFile( fname, options );
	if(fname.substr(0,4) == "http")
	   filename = fname;
	if (filename.empty()) return ReadResult::FILE_NOT_FOUND;

	std::string dir = getDirectory(filename);
	if(dir == "")
		return ReadResult::FILE_NOT_FOUND;

	osgEarth::ReadResult r = osgEarth::URI(filename).readString();
	if ( r.failed() )
	{
		//OE_WARN << LC << "Failed to read TMS tile map file from " << metadataXML << std::endl;
		return 0L;
	}
	std::stringstream buf( r.getString() );
	//std::ifstream fin(buf);
	char str[1024];  
	std::vector<std::string> lines;
	osgDB::ReaderWriter* reader = osgDB::Registry::instance()->getReaderWriterForExtension("osgb");
	while( buf.getline(str,1024) )
	{    
		std::string line = str;
		if(line.length() < 5)
		   continue;
		lines.push_back(line);

	}
	return readNodeFile(lines, dir);
}

///////////////////////////////////////////////////////////////////////////
// Add ourself to the Registry to instantiate the reader/writer.

REGISTER_OSGPLUGIN(desc, ReaderWriterDESC)

// vim: set sw=4 ts=8 et ic ai:
