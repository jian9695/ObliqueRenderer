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

void compress(std::string infile)
{

	const char *pszFormat = "GTiff";
	GDALDriver *poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);

	std::string outfile = (QFileInfo(infile.data()).absoluteDir().absolutePath() + "/").toLocal8Bit().data() + std::string("tmp.tif");
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
	poDriver->Delete(infile.data());
	rename(outfile.data(), infile.data());
	//std::string commandline = "move " + outfile + " " + infile;
	//printf("%s\n", commandline.data());
	//system(commandline.data());
}

int main(int argc, char** argv)
{
	std::string indir = QDir::currentPath().toLocal8Bit().data();
	if (argc > 1)
	{
		indir = (QDir(argv[1]).absolutePath() + "/").toLocal8Bit().data();
	}
	printf("%s\n", indir.data());
	OGRRegisterAll();
	GDALAllRegister();
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		if (!fileInfo.fileName().endsWith(".tif"))
			continue;
		printf("%s\n", fileInfo.absoluteFilePath().toLocal8Bit().data());
		compress(fileInfo.absoluteFilePath().toLocal8Bit().data());
	}

}
