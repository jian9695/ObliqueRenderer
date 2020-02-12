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
#include <gdal_priv.h>

struct RasterFileTransform
{
	std::string name;
	double adfGeoTransform[6];
	int xsize;
	int ysize;
};
void readTransform(const char* filename, double* adfGeoTransform)
{
	GDALDataset *poSrcDS = (GDALDataset*)GDALOpen(filename, GA_ReadOnly);
	poSrcDS->GetGeoTransform(adfGeoTransform);
	GDALClose(poSrcDS);
}

OGREnvelope getBoundFromFile(const char* rasterfile,double* adfGeoTransform,int& xsize,int& ysize)
{

	GDALDataset  *poDataset = (GDALDataset *)GDALOpen(rasterfile, GA_ReadOnly);
	GDALRasterBand *poBand = poDataset->GetRasterBand(4);
	xsize = poBand->GetXSize();
	ysize = poBand->GetYSize();
	int   ncells = xsize*ysize;
	unsigned char* data = new unsigned char[ncells];
	poBand->RasterIO(GF_Read, 0, 0, xsize, ysize,
		data, xsize, ysize, GDT_Byte,
		0, 0);

	poDataset->GetGeoTransform(adfGeoTransform);

	GDALClose((GDALDatasetH)poDataset);
	OGREnvelope bound;
	
	for (long yidx = 0; yidx < ysize; yidx++)
	{
		long num = 0;
		for (long xidx = 0; xidx < xsize; xidx++)
		{
			unsigned char alpha = data[xidx + yidx * xsize];
			if (alpha > 125)
				num++;
		}
		if (num > xsize / 2)
		{
			bound.MaxY = adfGeoTransform[3] + (yidx + 1) * adfGeoTransform[5];
			break;
		}

	}

	for (long yidx = ysize - 1; yidx >= 0; yidx--)
	{
		long num = 0;
		for (long xidx = 0; xidx < xsize; xidx++)
		{
			unsigned char alpha = data[xidx + yidx * xsize];
			if (alpha > 125)
				num++;
		}
		if (num > xsize / 2)
		{
			bound.MinY = adfGeoTransform[3] + (yidx + 1) * adfGeoTransform[5];
			break;
		}

	}

	for (long xidx = 0; xidx < xsize; xidx++)
	{
		long num = 0;
		for (long yidx = 0; yidx < ysize; yidx++)
		{
			unsigned char alpha = data[xidx + yidx * xsize];
			if (alpha > 125)
				num++;
		}
		if (num > ysize / 2)
		{
			bound.MinX = adfGeoTransform[0] + (xidx + 1) * adfGeoTransform[1];
			break;
		}

	}

	for (long xidx = xsize - 1; xidx >= 0; xidx--)
	{
		long num = 0;
		for (long yidx = 0; yidx < ysize; yidx++)
		{
			unsigned char alpha = data[xidx + yidx * xsize];
			if (alpha > 125)
				num++;
		}
		if (num > ysize / 2)
		{
			bound.MaxX = adfGeoTransform[0] + (xidx + 1) * adfGeoTransform[1];
			break;
		}

	}
	delete[] data;
	return bound;
}

OGREnvelope getBoundFromDir(std::string& inputdir, std::vector<RasterFileTransform>& filenames)
{
	OGREnvelope bound;
	bound.MinX = 9000000000;
	bound.MaxX = -9000000000;
	bound.MinY = 9000000000;
	bound.MaxY = -9000000000;

	std::vector<std::string> files;
	QDir input_dir(inputdir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	inputdir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".tif"))
			continue;
		RasterFileTransform rasterft;
		OGREnvelope subbound = getBoundFromFile(input_file.data(), rasterft.adfGeoTransform,rasterft.xsize,rasterft.ysize);
		if (subbound.MinX < bound.MinX)
			bound.MinX = subbound.MinX;
		if (subbound.MinY < bound.MinY)
			bound.MinY = subbound.MinY;
		if (subbound.MaxX > bound.MaxX)
			bound.MaxX = subbound.MaxX;
		if (subbound.MaxY > bound.MaxY)
			bound.MaxY = subbound.MaxY;
		rasterft.name = fileInfo.fileName().toLocal8Bit().data();
		filenames.push_back(rasterft);
	}

	return bound;
}
void crop(std::string inputdir, std::string outputdir)
{
	std::vector<RasterFileTransform> filenames;
	OGREnvelope bound = getBoundFromDir(inputdir, filenames);
	getBoundFromDir(inputdir, filenames);
	QDir qoutdir(outputdir.data());
	outputdir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();
	if (!qoutdir.exists())
		qoutdir.mkpath(".");
//	gdal_translate - projwin 23.0 53.0 34.0 23.0 NETCDF:"E:\FFDAS\input\d001.nc" : flux_h05
	for (size_t i = 0; i < filenames.size(); i++)
	{
		std::stringstream ss;
		ss << "gdal_translate -srcwin ";

		int x1 = (int)((bound.MinX - filenames[i].adfGeoTransform[0]) / filenames[i].adfGeoTransform[1]);
		if (x1 < 0)
			x1 = 0;
		int y1 = (int)((bound.MaxY - filenames[i].adfGeoTransform[3]) / filenames[i].adfGeoTransform[5]);
		if (y1 < 0)
			y1 = 0;
		
		int x2 = (int)((bound.MaxX - filenames[i].adfGeoTransform[0]) / filenames[i].adfGeoTransform[1]);
		if (x2 > filenames[i].xsize - 1)
			x2 = filenames[i].xsize - 1;
		int y2 = (int)((bound.MinY - filenames[i].adfGeoTransform[3]) / filenames[i].adfGeoTransform[5]);
		if (y2 > filenames[i].ysize - 1)
			y2 = filenames[i].ysize - 1;
		int xoffset = x1;
		int yoffset = y1;
		int xsize = x2 - x1;
		int ysize = y2 - y1;

		ss << xoffset << " " << yoffset << " " << xsize << " " << ysize << " -of GTiff ";
		ss << inputdir + filenames[i].name << " " << outputdir + filenames[i].name;
		printf("%s", ss.str().data());
		system(ss.str().data());
	}
}

struct ObliqueAngle
{
	double alt;
	double azimuth;
	ObliqueAngle(double _alt, double _azimuth)
	{
		alt = _alt;
		azimuth = _azimuth;
	}
};

bool removeDir(const QString &dirName)
{
	bool result = true;
	QDir dir(dirName);

	if (dir.exists(dirName)) {
		Q_FOREACH(QFileInfo info, dir.entryInfoList(QDir::NoDotAndDotDot | QDir::System | QDir::Hidden | QDir::AllDirs | QDir::Files, QDir::DirsFirst)) {
			if (info.isDir()) {
				result = removeDir(info.absoluteFilePath());
			}
			else {
				result = QFile::remove(info.absoluteFilePath());
			}

			if (!result) {
				return result;
			}
		}
		//result = dir.rmdir(dirName);
	}

	return result;
}
void compress(std::string infile, std::string outfile)
{
	const char *pszFormat = "GTiff";
	GDALDriver *poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);

	//std::string outfile = (QFileInfo(infile.data()).absoluteDir().absolutePath() + "/").toLocal8Bit().data() + std::string("tmp.tif");
	GDALDataset *poSrcDS = (GDALDataset *)GDALOpen(infile.data(), GA_ReadOnly);
	char **papszOptions = NULL;
	papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");
	papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "LZW");
	GDALDataset * poDstDS = poDriver->CreateCopy(outfile.data(), poSrcDS, FALSE,
		papszOptions, GDALTermProgress, NULL);
	/* Once we're done, close properly the dataset */
	if (poDstDS != NULL)
		GDALClose((GDALDatasetH)poDstDS);
	CSLDestroy(papszOptions);
	GDALClose((GDALDatasetH)poSrcDS);
	if (!QFileInfo(outfile.data()).exists())
		return;
	//poDriver->Delete(infile.data());
	//rename(outfile.data(), infile.data());
}

int main(int argc, char** argv)
{
	GDALAllRegister();


	std::string modelfile = argv[1];
	std::string shapefile = argv[2];
	std::string outdir = argv[3];

	double resol = atof(argv[4]);
	bool boundSpecified = false;
	double xmin, ymin, xmax, ymax;
	if (argc > 7)
	{
		xmin = atof(argv[7]); xmax = atof(argv[8]);
		ymin = atof(argv[9]); ymax = atof(argv[10]);
		boundSpecified = true;
	}

	double altAngle = 45;
	QDir qoutdir(outdir.data());
	if(!qoutdir.exists())
	{
		qoutdir.mkpath(".");
	}
	outdir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();


	std::string temp_outdir = (qoutdir.absolutePath() + "/" + "tmp" + "/").toLocal8Bit().data();
	QDir(temp_outdir.data()).mkpath(".");

	ObliqueAngle obliqueAngle(atof(argv[5]), atof(argv[6]));

	std::stringstream ssName;
	ssName << obliqueAngle.alt << "_" << obliqueAngle.azimuth;
	std::string name = ssName.str();

	std::stringstream ssTileDir;
	ssTileDir << outdir << name << "/";

	//Oblique_Renderer "F:/Oblique_Photogrammetry/xinji/Data/all.desc" "xinji/25cm/45_90" 0.25 45 90
	std::stringstream Oblique_Renderer_Command;
	Oblique_Renderer_Command << "Oblique_Renderer " << modelfile << " " << ssTileDir.str() << " ";
	Oblique_Renderer_Command << resol << " " << obliqueAngle.alt << " " << obliqueAngle.azimuth << " ";

	std::stringstream Oblique_GCP_Command;
	Oblique_GCP_Command << "Oblique_GCP " << modelfile << " " << ssTileDir.str() << " ";
	Oblique_GCP_Command << resol << " " << obliqueAngle.alt << " " << obliqueAngle.azimuth << " " << shapefile << " ";

	if (boundSpecified)
	{
		Oblique_Renderer_Command << xmin << " " << xmax << " " << ymin << " " << ymax;
		Oblique_GCP_Command << xmin << " " << xmax << " " << ymin << " " << ymax;
	}

	std::string vrtFile = ssTileDir.str() + "oblique.vrt";
	if (!QFileInfo(vrtFile.data()).exists())
	{
		printf("%s\n", Oblique_Renderer_Command.str().data());
		system(Oblique_Renderer_Command.str().data());
	}

	if (!QFileInfo(vrtFile.data()).exists())
		return 0;

	std::string gcpFile = ssTileDir.str() + "georeference.txt";
	if (!QFileInfo(gcpFile.data()).exists())
	{
		printf("%s\n", Oblique_GCP_Command.str().data());
		system(Oblique_GCP_Command.str().data());
	}
	if (!QFileInfo(gcpFile.data()).exists())
		return 0;
	std::string warppedFile = temp_outdir + name + ".tif";

	std::stringstream Oblique_Warp_Command;
	Oblique_Warp_Command << "Oblique_Warp " << vrtFile << " " << warppedFile << " " << gcpFile;
	std::string outFile = outdir + name + ".tif";
	if (!QFileInfo(warppedFile.data()).exists())
	{
		printf("%s\n", Oblique_Warp_Command.str().data());
		system(Oblique_Warp_Command.str().data());
	}
	if (!QFileInfo(outFile.data()).exists())
	{
		compress(warppedFile, outFile);
	}
	//std::stringstream Oblique_Crop_Command;
	//Oblique_Crop_Command << "Oblique_Crop " << temp_outdir << " " << outdir;
	//printf("%s", Oblique_Crop_Command.str().data());
	//system(Oblique_Crop_Command.str().data());

	//removeDir(temp_outdir.data());

	/*std::stringstream ssTileDir;
	ssTileDir << outdir << obliqueAngle.alt << "_" << obliqueAngle.azimuth << "/";*/
	//removeDir(ssTileDir.str().data());

}
