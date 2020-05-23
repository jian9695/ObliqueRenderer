#pragma once
#include "ogrsf_frmts.h"

class ShapeFile
{
public:
	ShapeFile();
	ShapeFile(std::string filename, int update = 0);
	~ShapeFile();
	void close();
	int getOrCreateField(const char* name, OGRFieldType tp);
	int getField(const char* name);
	void create(std::string filename, OGRSpatialReference* spatialRef = NULL, OGRFeatureDefn *poFDefn = NULL, OGRwkbGeometryType geotype = wkbPolygon);
	static void move(std::string src, std::string dest);
	static void copy(std::string src, std::string dest);
	static void remove(std::string filename);
	static void copyDropFields(std::string src, std::string dest, std::vector<std::string> fields2drop);
	static void copyDirDropFields(std::string srcDir, std::string destDir, std::vector<std::string> fields2drop);
	static double getTotal(std::string shapefile, std::string fieldname);
	static double getTotal(std::vector<std::string> shapefiles, std::string fieldname);
	static void copyField(std::string shapefile, std::string oldfield, std::string newfield);
	static void reproject(std::string inputfile, std::string outputfile, std::string destwktfile);
	static void reprojectDir(std::string indir, std::string outdir, std::string destwktfile);
	static void changeFieldNames(std::string filename, std::vector<std::string> oldnames, std::vector<std::string> newnames);
	static void copy(std::string src, std::string dest, std::vector<std::string> fields2keep);
	static void copyDir(std::string srcDir, std::string destDir, std::vector<std::string> fields2keep);
	static void copyDropGeometry(std::string src, std::string dest, std::vector<std::string> fields2keep);
	static void copyDirDropGeometry(std::string srcDir, std::string destDir, std::vector<std::string> fields2keep);
	static void scaleCAFields(std::string srcDir);
	static bool isNULL(std::string filename);
	static void mergeFile(std::vector<std::string> inputfiles, std::vector<std::string> fields2keep, std::string outputfile);
	double sum(std::string fieldname);
	std::vector<double> sum(std::vector<std::string> fieldnames);
public:
	OGRLayer       *poLayer;
	GDALDataset    *poDS;
  std::string proj;
private:
	std::string g_mFileName;
};

