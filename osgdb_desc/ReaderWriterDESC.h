#ifndef _READERWRITERDESC_H_
#define _READERWRITERDESC_H_


#include <OpenThreads/ReentrantMutex>

///////////////////////////////////////////////////////////////////////////
// OSG reader/writer plugin for the COLLADA 1.4.x ".dae" format.
// See http://collada.org/ and http://khronos.org/collada/

class ReaderWriterDESC : public osgDB::ReaderWriter
{
public:
    ReaderWriterDESC()
    {
        // Collada document
        supportsExtension("desc","Oblique City Model Description");
        // Collada zip archive (contains one or more dae files and a manifest.xml)
    }
    const char* className() const { return "Oblique City Model Description reader/writer"; }

   // ReadResult readNode(std::istream&, const Options* = NULL) const;
    ReadResult readNode(const std::string&, const Options* = NULL) const;
	ReadResult readNode(std::istream& fin, const Options* options) const;
	ReadResult readNodeFile(const std::vector<std::string>& lines, const std::string& dir) const;
private:
    mutable OpenThreads::ReentrantMutex _serializerMutex;
	std::string getDirectory(const std::string& filename) const;
	std::string replace(const std::string& name) const;
	std::vector<std::string> split(const char& delimiter,const std::string& line) const;
};


///////////////////////////////////////////////////////////////////////////

#endif

