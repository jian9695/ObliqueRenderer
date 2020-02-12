
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "gdal_priv.h"
#include <string>
#include "ogrsf_frmts.h"
#include "qfileinfo.h"
#include "qdir.h"
//adfGeoTransform[0] = _bb.xMin() + ncol * _cellsize;///* top left x */
//adfGeoTransform[1] = _cellsize / img->s();///* w-e pixel resolution */
//adfGeoTransform[2] = 0;///* 0 */
//adfGeoTransform[3] = _bb.yMax() - nrow * _cellsize;// /* top left y */
//adfGeoTransform[4] = 0;///* 0 */
//adfGeoTransform[5] = -_cellsize / img->t();///* n-s pixel resolution (negative value) */

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
	void create(const char* filename, OGRSpatialReference* spatialRef = NULL, OGRFeatureDefn *poFDefn = NULL, OGRwkbGeometryType geotype = wkbPolygon)
	{


		g_mFileName = filename;
		if (poDS)
			GDALClose(poDS);
		const char *pszDriverName = "ESRI Shapefile";
		GDALDriver *poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(
			pszDriverName);


		if (QFileInfo(filename).exists())
		{
			poDriver->Delete(filename);
			//poDS = OGRSFDriverRegistrar::Open(filename, TRUE);
			//poLayer = poDS->GetLayer(0);
		}


		poDS = poDriver->Create(filename, 0, 0, 0, GDT_Unknown, NULL);
		poLayer = poDS->CreateLayer(QFileInfo(filename).baseName().toLocal8Bit().data(), spatialRef, geotype, NULL);

		/*OGRSpatialReference oSRS;
		oSRS.SetWellKnownGeogCS("WGS84");
		poLayer = poDS->CreateLayer(QFileInfo(filename).baseName().toLocal8Bit().data(), &oSRS, tp, NULL);
		OGRSpatialReference* spatialRef = poLayer->GetSpatialRef();*/


		//else
		//{
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

void georeferenceRaster(std::string infile, std::string reffile, double shift_x, double shift_y)
{
	GDALDataset *srccDS = (GDALDataset *)GDALOpen(infile.data(), GA_Update);
	GDALDataset *reffDS = (GDALDataset *)GDALOpen(reffile.data(), GA_ReadOnly);

	//"F:/Oblique_Photogrammetry/weihai/Weihai_DOM_025.tif";
	double adfGeoTransform[6];
	srccDS->GetGeoTransform(adfGeoTransform);
	adfGeoTransform[0] += shift_x;
	adfGeoTransform[3] += shift_y;
	srccDS->SetGeoTransform(adfGeoTransform);
	if (reffDS)
	{
		srccDS->SetProjection(reffDS->GetProjectionRef());
		GDALClose((GDALDatasetH)reffDS);
	}

	GDALClose((GDALDatasetH)srccDS);
}
void georeferenceShape(std::string infile, std::string reffile, double shift_x, double shift_y)
{


	ShapeFile shp(infile,1);

	OGRFeature *poFeature;
	shp.poLayer->ResetReading();

	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		OGRPoint* point = (OGRPoint*)poFeature->GetGeometryRef();
		OGRPoint po;
		po.setX(point->getX() + shift_x);
		po.setY(point->getY() + shift_y);
		poFeature->SetGeometry(&po);
		shp.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}





}
void georeferenceRaster2()
{
	GDALDataset *srccDS = (GDALDataset *)GDALOpen("F:/Oblique_Photogrammetry/weihai/Weihai_DOM_025.tif", GA_Update);

	double adfGeoTransform[6];
	srccDS->GetGeoTransform(adfGeoTransform);
	adfGeoTransform[0] = 0;
	adfGeoTransform[3] = 0;
	srccDS->SetGeoTransform(adfGeoTransform);

	GDALClose((GDALDatasetH)srccDS);
}
int main(int argc, char** argv)
{
	GDALAllRegister();
	OGRRegisterAll();
	//georeferenceRaster2();
	//return 0;
	std::string infile = argv[1];
	std::string reffile = argv[2];
	double shift_x = atof(argv[3]);
	double shift_y = atof(argv[4]);

	if (QFileInfo(infile.data()).isFile() && QFileInfo(infile.data()).exists())
	{
		if (QFileInfo(infile.data()).fileName().endsWith(".shp"))
		{
			georeferenceShape(QFileInfo(infile.data()).absoluteFilePath().toLocal8Bit().data(), reffile, shift_x, shift_y);
			return 0;
		}
		georeferenceRaster(infile, reffile, shift_x, shift_y);
		return 0;
	}

	QDir input_dir(infile.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".tif"))
		{
			continue;
		}
		printf("%s\n", fileInfo.absoluteFilePath().toLocal8Bit().data());
		georeferenceRaster(fileInfo.absoluteFilePath().toLocal8Bit().data(), reffile, shift_x, shift_y);

	}

}
