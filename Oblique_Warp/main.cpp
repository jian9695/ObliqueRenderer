// PolygonStatistics.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "ogrsf_frmts.h"
//#include "gdal_priv.h"
#include <QFileinfo>
#include <map>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>      // std::setprecision
#include <qdir.h>
#include "qfileinfo.h"
double FOOTPRINT_SCALE_FACTOR = 1;
#include <gdal_priv.h>


/*
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# intersect.py
# Created on : 2016 - 02 - 23 10 : 09:07.00000
#   (generated by ArcGIS/ModelBuilder)
# Description:
# ---------------------------------------------------------------------------

# Import arcpy module
import arcpy


# Local variables :
fishnet_shp__2_ = "B:\\Baltimore\\gridPrep_SHP_master\\WGS84\\result\\fishnet.shp"
CMVUnderway_shp = "B:\\Baltimore\\gridPrep_SHP_master\\WGS84\\CMVUnderway.shp"
sdfs_shp = "B:\\Baltimore\\gridPrep_SHP_master\\sdfs.shp"

# Process : Intersect
arcpy.Intersect_analysis("B:\\Baltimore\\gridPrep_SHP_master\\WGS84\\result\\fishnet.shp;B:\\Baltimore\\gridPrep_SHP_master\\WGS84\\CMVUnderway.shp", "B:\\Baltimore\\gridPrep_SHP_master\\sdfs.shp", "ALL", "", "INPUT")
*/

//# ---------------------------------------------------------------------------
//# warp.py
//# Created on : 2016 - 03 - 10 00 : 32 : 11.00000
//#   (generated by ArcGIS/ModelBuilder)
//# Description :
//# ---------------------------------------------------------------------------
//
//	# Import arcpy module
//	import arcpy
//
//
//	# Local variables :
//oblique_vrt = "B:\\Oblique_Renderer\\25CM\\25CM_45_0\\oblique.vrt"
//oblique_Warp = "C:\\Users\\Administrator\\Documents\\ArcGIS\\Default.gdb\\oblique_Warp"
//
//# Process : Warp
//arcpy.Warp_management(oblique_vrt, "'23 32';'453 454';'34 656';'2332 43'", "'2343 53';'12 67';'45 65';'645 65'", oblique_Warp, "POLYORDER3", "NEAREST")
////GDALDataset *poDstDS = (GDALDataset *)GDALOpen(ssDest.str().data(), GA_Update);
//double adfGeoTransform[6];
//adfGeoTransform[0] = _bb.xMin() + ncol * _cellsize;///* top left x */
//adfGeoTransform[1] = _cellsize / img->s();///* w-e pixel resolution */
//adfGeoTransform[2] = 0;///* 0 */
//adfGeoTransform[3] = _bb.yMax() - nrow * _cellsize;// /* top left y */
//adfGeoTransform[4] = 0;///* 0 */
//adfGeoTransform[5] = -_cellsize / img->t();///* n-s pixel resolution (negative value) */
struct GCP
{
	double xsource;
	double ysource;
	double xmap;
	double ymap;
	double x;//pixel in image space
	double y;//line in image space
	void getPixelLine(double* adfGeoTransform)
	{
		x = (xsource - adfGeoTransform[0]) / adfGeoTransform[1];
		y = (ysource - adfGeoTransform[3]) / adfGeoTransform[5];
	}
};
void warpFromFile(std::string inputfile, std::string outputfile, std::string gcpfile)
{

	std::ifstream ifs(gcpfile);
	std::string line;
	std::vector<GCP> gcps;
	while (getline(ifs, line)) {
		std::stringstream ss(line);
		GCP gcp;
		ss >> gcp.xsource >> gcp.ysource >> gcp.xmap >> gcp.ymap;
		gcps.push_back(gcp);
	}
	ifs.close();


	QFileInfo info(outputfile.data());
	if (info.exists()) {
		const char *pszFormat = "GTiff";
		GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
		poDriver->Delete(outputfile.data());
	}
	std::ofstream ofs;
	std::string scriptFile = (info.absoluteDir().absolutePath() + "/" + info.completeBaseName() + ".py").toLocal8Bit().data();
	ofs.open(scriptFile.data());
	std::stringstream script;
	script << "import arcpy" << "\n";
	//script << "inputfile = " << "\"" << inputfile << "\"" << "\n";
	//script << "outputfile = " << "\"" << outputfile << "\"" << "\n";
	////arcpy.Warp_management(oblique_vrt, "'23 32';'453 454';'34 656';'2332 43'", "'2343 53';'12 67';'45 65';'645 65'", oblique_Warp, "POLYORDER3", "NEAREST")
	//script << "arcpy.Warp_management(";
	////script << "inputfile1" << ";" << "inputfile2" << "," << "outputfile" << "," << "\"ALL\"" << "," << "\"\"" << ","  "\"INPUT\")";
	////script << "\"" << inputfile1 << ";" << inputfile2 << "\"" << "," << "outputfile" << "," << "\"NO_FID\"" << "," << "\"\"" << ","  "\"INPUT\")";
	//script << "inputfile" << ",";
	//script << "\"";
	//for (size_t i = 0; i < gcps.size(); i++)
	//{
	//	GCP& gcp = gcps[i];
	//	script << "'" << gcp.xsource << " " << gcp.ysource << "'";
	//	if (i != gcps.size() - 1)
	//		script << ";";
	//}
	//script << "\"" << ",";
	//script << "\"";
	//for (size_t i = 0; i < gcps.size(); i++)
	//{
	//	GCP& gcp = gcps[i];
	//	script << "'" << gcp.xmap << " " << gcp.ymap << "'";
	//	if (i != gcps.size() - 1)
	//		script << ";";
	//}
	//script << "\"" << ",";
	//script << "outputfile" << ",";
	//script << "\"" << "POLYORDER3" << "\"" << ",";
	//script << "\"" << "NEAREST" << "\"" << ")";

	script << "inputfile = " << "\"" << inputfile << "\"" << "\n";
	script << "outputfile = " << "\"" << outputfile << "\"" << "\n";
	script << "gcpfile = " << "\"" << gcpfile << "\"" << "\n";

	//arcpy.Warp_management(oblique_vrt, "'23 32';'453 454';'34 656';'2332 43'", "'2343 53';'12 67';'45 65';'645 65'", oblique_Warp, "POLYORDER3", "NEAREST")
	script << "arcpy.WarpFromFile_management(";
	//script << "inputfile1" << ";" << "inputfile2" << "," << "outputfile" << "," << "\"ALL\"" << "," << "\"\"" << ","  "\"INPUT\")";
	//script << "\"" << inputfile1 << ";" << inputfile2 << "\"" << "," << "outputfile" << "," << "\"NO_FID\"" << "," << "\"\"" << ","  "\"INPUT\")";
	script << "inputfile" << ",";
	script << "outputfile" << ",";
	script << "gcpfile" << ",";
	script << "\"" << "POLYORDER1" << "\"" << ",";
	script << "\"" << "NEAREST" << "\"" << ")";
	//arcpy.WarpFromFile_management
	ofs << script.str();
	ofs.close();
	printf("%s\n", script.str().data());
	system(scriptFile.data());

	QFile::remove(scriptFile.data());
}
void warp(std::string inputfile, std::string outputfile, std::string gcpfile)
{

	std::ifstream ifs(gcpfile);
	std::string line;
	std::vector<GCP> gcps;
	while (getline(ifs, line)) {
		std::stringstream ss(line);
		GCP gcp;
		ss >> gcp.xsource >> gcp.ysource >> gcp.xmap >> gcp.ymap;
		gcps.push_back(gcp);
	}
	ifs.close();


	QFileInfo info(outputfile.data());
	if (info.exists()) {
		const char *pszFormat = "GTiff";
		GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
		poDriver->Delete(outputfile.data());
	}
	std::ofstream ofs;
	std::string scriptFile = (info.absoluteDir().absolutePath() + "/" + info.completeBaseName() + ".py").toLocal8Bit().data();
	ofs.open(scriptFile.data());
	std::stringstream script;
	script << "import arcpy" << "\n";
	script << "inputfile = " << "\"" << inputfile << "\"" << "\n";
	script << "outputfile = " << "\"" << outputfile << "\"" << "\n";
	//arcpy.Warp_management(oblique_vrt, "'23 32';'453 454';'34 656';'2332 43'", "'2343 53';'12 67';'45 65';'645 65'", oblique_Warp, "POLYORDER3", "NEAREST")
	script << "arcpy.Warp_management(";
	//script << "inputfile1" << ";" << "inputfile2" << "," << "outputfile" << "," << "\"ALL\"" << "," << "\"\"" << ","  "\"INPUT\")";
	//script << "\"" << inputfile1 << ";" << inputfile2 << "\"" << "," << "outputfile" << "," << "\"NO_FID\"" << "," << "\"\"" << ","  "\"INPUT\")";
	script << "inputfile" << ",";
	script << "\"";
	for (size_t i = 0; i < gcps.size(); i++)
	{
		GCP& gcp = gcps[i];
		script << "'" << gcp.xsource << " " << gcp.ysource << "'";
		if (i != gcps.size() - 1)
			script << ";";
	}
	script << "\"" << ",";
	script << "\"";
	for (size_t i = 0; i < gcps.size(); i++)
	{
		GCP& gcp = gcps[i];
		script << "'" << gcp.xmap << " " << gcp.ymap << "'";
		if (i != gcps.size() - 1)
			script << ";";
	}
	script << "\"" << ",";
	script << "outputfile" << ",";
	script << "\"" << "POLYORDER3" << "\"" << ",";
	script << "\"" << "NEAREST" << "\"" << ")";

	//script << "inputfile = " << "\"" << inputfile << "\"" << "\n";
	//script << "outputfile = " << "\"" << outputfile << "\"" << "\n";
	//script << "gcpfile = " << "\"" << gcpfile << "\"" << "\n";

	////arcpy.Warp_management(oblique_vrt, "'23 32';'453 454';'34 656';'2332 43'", "'2343 53';'12 67';'45 65';'645 65'", oblique_Warp, "POLYORDER3", "NEAREST")
	//script << "arcpy.WarpFromFile_management(";
	////script << "inputfile1" << ";" << "inputfile2" << "," << "outputfile" << "," << "\"ALL\"" << "," << "\"\"" << ","  "\"INPUT\")";
	////script << "\"" << inputfile1 << ";" << inputfile2 << "\"" << "," << "outputfile" << "," << "\"NO_FID\"" << "," << "\"\"" << ","  "\"INPUT\")";
	//script << "inputfile" << ",";
	//script << "outputfile" << ",";
	//script << "gcpfile" << ",";
	//script << "\"" << "POLYORDER1" << "\"" << ",";
	//script << "\"" << "NEAREST" << "\"" << ")";
	//arcpy.WarpFromFile_management
	ofs << script.str();
	ofs.close();
	printf("%s\n", script.str().data());
	system(scriptFile.data());

	QFile::remove(scriptFile.data());
}


void readTransform(const char* filename, double* adfGeoTransform)
{
	GDALDataset *poSrcDS = (GDALDataset*)GDALOpen(filename, GA_ReadOnly);
	poSrcDS->GetGeoTransform(adfGeoTransform);
	GDALClose(poSrcDS);
}
void walpGDAL(std::string inputfile, std::string outputfile, std::string gcpfile)
{
	std::string outdir = (QFileInfo(outputfile.data()).absoluteDir().absolutePath() + "/").toLocal8Bit().data();

	std::string outputfile_tmp = outdir + (QFileInfo(outputfile.data()).completeBaseName() + "_2" + ".tif").toLocal8Bit().data();

	const char *pszFormat = "GTiff";
	GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);

	std::ifstream ifs(gcpfile);
	std::string line;
	std::vector<GCP> gcps;

	double adfGeoTransform[6];
	readTransform(inputfile.data(), adfGeoTransform);



	while (getline(ifs, line)) {
		std::stringstream ss(line);
		GCP gcp;
		ss >> gcp.xsource >> gcp.ysource >> gcp.xmap >> gcp.ymap;
		gcp.getPixelLine(adfGeoTransform);
		gcps.push_back(gcp);
	}
	ifs.close();

	if (QFileInfo(outputfile.data()).exists()) {
		poDriver->Delete(outputfile.data());
	}


	std::stringstream gdal_translate_script;
	gdal_translate_script << "gdal_translate -of GTiff ";
	for (size_t i = 0; i < gcps.size(); i++)
	{
		GCP& gcp = gcps[i];
		//if (i % 2 != 0)
		//	continue;
		gdal_translate_script << "-gcp " << gcp.x << " " << gcp.y << " " << gcp.xmap << " " << gcp.ymap << " ";
	}

	gdal_translate_script << "\"" << inputfile << "\"" << " ";
	gdal_translate_script << "\"" << outputfile_tmp << "\"";
	printf("%s\n", gdal_translate_script.str().data());
	system(gdal_translate_script.str().data());

	std::stringstream gdal_warp_script;
	//gdal_warp_script << "gdalwarp -r cubic -order 3 -co COMPRESS=LZW ";
	//gdal_warp_script << "\"" << outputfile_tmp << "\"" << " ";
	//gdal_warp_script << "\"" << outputfile << "\"";
	//printf("%s\n",gdal_warp_script.str().data());

	gdal_warp_script << "gdalwarp -r near -order 3 ";
	gdal_warp_script << "\"" << outputfile_tmp << "\"" << " ";
	gdal_warp_script << "\"" << outputfile << "\"";
	//gdal_warp_script << " -co COMPRESS=LZW";

	
	printf("%s\n", gdal_warp_script.str().data());


	system(gdal_warp_script.str().data());

	if (QFileInfo(outputfile_tmp.data()).exists()) {
		poDriver->Delete(outputfile_tmp.data());
	}
	//QFile::remove(outputfile_tmp.data());
/*gdal_translate -of GTiff -gcp -4019.01 7059.03 32 323 -gcp -3971.75 7051.42 23 53 -gcp -3987.5 7066.91 3443 43 -gcp -3988.82 7048.79 56 45 -gcp -4004.31 7040.39 45 34 -gcp -4024 7039.08 435 6 -gcp -3990.66 7071.64 45 66 -gcp -3988.29 7101.57 53 434 -gcp -4015.86 7044.85 54 434 -gcp -4002.73 7041.44 45 54 -gcp -3986.72 7075.84 65 75 -gcp -4027.41 7064.28 56 54 "B:/VGE/DSMRenderer_2016_02_14/bin/10_45_0/11_13.tif" "C:/Users/jliang41/AppData/Local/Temp/11_13.tif"
gdalwarp -r bilinear -order 3 -co COMPRESS=LZW  "C:/Users/jliang41/AppData/Local/Temp/11_13.tif" "B:/VGE/DSMRenderer_2016_02_14/bin/10_45_0/11_13_modified.tif"*/
}
int main(int argc, char** argv)
{
	std::string infile = argv[1];
	std::string outfile = argv[2];
	std::string gcpfile = argv[3];
	OGRRegisterAll();
	GDALAllRegister();
	//warp(infile, outfile, gcpfile);

	walpGDAL(infile, outfile, gcpfile);
}
