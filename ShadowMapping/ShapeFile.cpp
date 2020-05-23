#include "ShapeFile.h"
#include "qfileinfo.h"
#include <fstream>
#include <sstream>
#include "qdir.h"

void ShapeFile::close()
{
	if (poDS)
		GDALClose(poDS);
	poDS = NULL;
}
void ShapeFile::create(std::string filename, OGRSpatialReference* spatialRef, OGRFeatureDefn *poFDefn, OGRwkbGeometryType geotype)
{


	g_mFileName = filename;
	if (poDS)
		GDALClose(poDS);
	const char *pszDriverName = "ESRI Shapefile";
	GDALDriver *poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(
		pszDriverName);

	if (QFileInfo(filename.data()).exists())
	{
		poDriver->Delete(filename.data());
	}

	poDS = poDriver->Create(filename.data(), 0, 0, 0, GDT_Unknown, NULL);

	//if (!spatialRef)
	//{
	//	OGRSpatialReference oSRS;
	//	oSRS.SetWellKnownGeogCS("WGS84");
	//	poLayer = poDS->CreateLayer(QFileInfo(filename.data()).baseName().toLocal8Bit().data(), &oSRS, geotype, NULL);
	//}
	//else
	//{
	//	poLayer = poDS->CreateLayer(QFileInfo(filename.data()).baseName().toLocal8Bit().data(), spatialRef, geotype, NULL);
	//}

	poLayer = poDS->CreateLayer(QFileInfo(filename.data()).baseName().toLocal8Bit().data(), spatialRef, geotype, NULL);

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
int ShapeFile::getOrCreateField(const char* name, OGRFieldType tp) {
	int idx = poLayer->GetLayerDefn()->GetFieldIndex(name);
	if (idx > -1)
		return idx;
	OGRFieldDefn field(name, tp);
	poLayer->CreateField(&field);
	return poLayer->GetLayerDefn()->GetFieldIndex(name);
}

int ShapeFile::getField(const char* name) {
	return poLayer->GetLayerDefn()->GetFieldIndex(name);
}
double ShapeFile::getTotal(std::string shapefile, std::string fieldname)
{
	ShapeFile shp(shapefile);
	OGRFeature *poFeature;
	double sum = 0;
	int idx = shp.poLayer->GetLayerDefn()->GetFieldIndex(fieldname.data());
	if (idx < 0)
		return sum;
	
	shp.poLayer->ResetReading();
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		double value = poFeature->GetFieldAsDouble(idx);
		sum += value;
		OGRFeature::DestroyFeature(poFeature);
	}
	return sum;
}

double ShapeFile::getTotal(std::vector<std::string> shapefiles,std::string fieldname)
{
	double total = 0;
	for (size_t i = 0; i < shapefiles.size(); i++)
	{
		total += getTotal(shapefiles[i], fieldname);
	}
	return total;
}


void ShapeFile::copyField(std::string shapefile, std::string oldfield, std::string newfield)
{
	ShapeFile shp(shapefile,1);
	OGRFeature *poFeature;
	double sum = 0;
	int oldidx = shp.poLayer->GetLayerDefn()->GetFieldIndex(oldfield.data());
	if (oldidx < 0)
		return ;
	int newidx = shp.poLayer->GetLayerDefn()->GetFieldIndex(newfield.data());
	if (newidx < 0)
	{
		OGRFieldDefn* oldfielddef = shp.poLayer->GetLayerDefn()->GetFieldDefn(oldidx);
		OGRFieldDefn newfielddef(newfield.data(),oldfielddef->GetType());
		shp.poLayer->CreateField(&newfielddef);
		newidx = shp.poLayer->GetLayerDefn()->GetFieldIndex(newfield.data());
	}
	OGRFieldType tp = shp.poLayer->GetLayerDefn()->GetFieldDefn(oldidx)->GetType();

	shp.poLayer->ResetReading();
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		poFeature->SetField(newidx, poFeature->GetRawFieldRef(oldidx));
		shp.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}

}

ShapeFile::ShapeFile()
{
	poDS = NULL;
}

ShapeFile::ShapeFile(std::string filename, int update) {
	g_mFileName = filename;
	poDS = NULL;
	//poDS = OGRSFDriverRegistrar::Open(filename.data(), update);
	poDS = (GDALDataset*)GDALOpenEx(filename.data(), GDAL_OF_VECTOR | update, NULL, NULL, NULL);
	poLayer = poDS->GetLayer(0);
}

void ShapeFile::move(std::string src, std::string dest)
{
	const char *pszDriverName = "ESRI Shapefile";
	GDALDriver *poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(
		pszDriverName);
	poDriver->Delete(dest.data());
	GDALDataset* poDSSrc = (GDALDataset*)GDALOpenEx(src.data(), GDAL_OF_VECTOR, NULL, NULL, NULL);
	GDALDataset* poDSDest = (GDALDataset*)poDriver->CreateCopy(dest.data(), poDSSrc,0,0,0,0);
	GDALClose(poDSSrc);
	GDALClose(poDSDest);
	poDriver->Delete(src.data());
}
void ShapeFile::copy(std::string src, std::string dest)
{
	//const char *pszDriverName = "ESRI Shapefile";
	//GDALDriver *poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(
	//	pszDriverName);
	//poDriver->Delete(dest.data());
	//GDALDataset* poDSSrc = (GDALDataset*)GDALOpenEx(src.data(), GDAL_OF_VECTOR, NULL, NULL, NULL);
	//GDALDataset* poDSDest = (GDALDataset*)poDriver->CreateCopy(dest.data(), poDSSrc, 0, 0, 0, 0);
	//GDALClose(poDSSrc);
	//GDALClose(poDSDest);
	QFileInfo fileinfo(src.data());
	std::string srcname = fileinfo.completeBaseName().toLocal8Bit().data();

	std::string destname = QFileInfo(dest.data()).completeBaseName().toLocal8Bit().data();
	std::string exts[] = { ".shp",".dbf",".prj",".shx",".shp.xml",".cpg" };
	std::string srcdir = fileinfo.absoluteDir().absolutePath().toLocal8Bit().data();
	std::string destdir = QFileInfo(dest.data()).absoluteDir().absolutePath().toLocal8Bit().data();


	for (size_t i = 0; i < 6; i++)
	{
		std::string srcfile = srcdir + "/" + srcname + exts[i];
		std::string destfile = destdir + "/" + destname + exts[i];

		if (!QFileInfo(srcfile.data()).exists())
			continue;
		if (!QFileInfo(destfile.data()).exists())
		{
			QFile::remove(destfile.data());
		}
		QFile::copy(srcfile.data(), destfile.data());
	}


}
void ShapeFile::remove(std::string filename)
{
	const char *pszDriverName = "ESRI Shapefile";
	GDALDriver *poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(
		pszDriverName);
	poDriver->Delete(filename.data());
}
void ShapeFile::copyDropFields(std::string src, std::string dest, std::vector<std::string> fields2drop)
{
	ShapeFile srcshp(src, 0);
	std::vector<int> oldfields;
	std::vector<int> newfields;
	bool found = false;
	for (size_t i = 0; i < srcshp.poLayer->GetLayerDefn()->GetFieldCount(); i++)
	{
		OGRFieldDefn* fielddef = srcshp.poLayer->GetLayerDefn()->GetFieldDefn(i);
		for (size_t j = 0; j < fields2drop.size(); j++)
		{
			if (fielddef->GetNameRef() == fields2drop[j])
			{
				found = true;
				break;
			}

		}
	}
	if (!found)
	{
		srcshp.close();
		copy(src, dest);
		return;
	}

	ShapeFile destshp;
	destshp.create(dest, srcshp.poLayer->GetSpatialRef(), 0, srcshp.poLayer->GetGeomType());



	for (size_t i = 0; i < srcshp.poLayer->GetLayerDefn()->GetFieldCount(); i++)
	{
		OGRFieldDefn* fielddef = srcshp.poLayer->GetLayerDefn()->GetFieldDefn(i);
		oldfields.push_back(srcshp.poLayer->GetLayerDefn()->GetFieldIndex(fielddef->GetNameRef()));
		bool found = false;
		for (size_t j = 0; j < fields2drop.size(); j++)
		{
			if (fielddef->GetNameRef() == fields2drop[j])
			{
				found = true;
				break;
			}

		}

		if (!found)
		{
			destshp.poLayer->CreateField(fielddef);
			newfields.push_back(destshp.poLayer->GetLayerDefn()->GetFieldIndex(fielddef->GetNameRef()));
		}
		else
		{
			newfields.push_back(-1);
		}
	}


	srcshp.poLayer->ResetReading();
	OGRFeature *poFeatureOld;
	while ((poFeatureOld = srcshp.poLayer->GetNextFeature()) != NULL)
	{
		OGRFeature *poFeatureNew = OGRFeature::CreateFeature(destshp.poLayer->GetLayerDefn());
		poFeatureNew->SetGeometry(poFeatureOld->GetGeometryRef());
		for (size_t i = 0; i < oldfields.size(); i++)
		{
			int oldfieldidx = oldfields[i];
			int newfieldidx = newfields[i];
			if (newfieldidx == -1)
				continue;

			poFeatureNew->SetField(newfieldidx, poFeatureOld->GetRawFieldRef(oldfieldidx));

		}

		destshp.poLayer->CreateFeature(poFeatureNew);
		OGRFeature::DestroyFeature(poFeatureNew);
		OGRFeature::DestroyFeature(poFeatureOld);
	}

}
void ShapeFile::copyDirDropFields(std::string srcDir, std::string destDir, std::vector<std::string> fields2drop)
{
	if (!QDir(srcDir.data()).exists())
		return;

	QDir qoutdir(destDir.data());
	if (!qoutdir.exists())
		qoutdir.mkpath(".");
	destDir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();
	std::vector<std::string> files;

	QDir input_dir(srcDir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	srcDir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = fileInfo.absoluteFilePath().toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		std::string output_file = destDir + fileInfo.fileName().toLocal8Bit().data();
		copyDropFields(input_file, output_file, fields2drop);
		printf("%s\n", output_file.data());
	}

}

void ShapeFile::changeFieldNames(std::string filename, std::vector<std::string> oldnames, std::vector<std::string> newnames)
{
	ShapeFile shapefile(filename.data(), 1);
	std::vector<int> fieldindices;
	OGRFeatureDefn* layerdef = shapefile.poLayer->GetLayerDefn();
	for (size_t i = 0; i < oldnames.size(); i++)
	{
		int idx = shapefile.poLayer->GetLayerDefn()->GetFieldIndex(oldnames[i].data());
		if (idx < 0)
			continue;
		layerdef->GetFieldDefn(idx)->SetName(newnames[i].data());
		std::string newname = layerdef->GetFieldDefn(idx)->GetNameRef();
		printf("%s\n,", newname.data());
		shapefile.poLayer->AlterFieldDefn(idx, layerdef->GetFieldDefn(idx), 1);
	}

	shapefile.close();
}

void ShapeFile::copy(std::string src, std::string dest,std::vector<std::string> fields2keep)
{
	ShapeFile srcshp(src, 0);
	ShapeFile destshp;
	destshp.create(dest, srcshp.poLayer->GetSpatialRef(), 0, srcshp.poLayer->GetGeomType());

	std::vector<int> oldfields;
	std::vector<int> newfields;


	if (fields2keep.size() == 0)
	{
		for (size_t i = 0; i < srcshp.poLayer->GetLayerDefn()->GetFieldCount(); i++)
		{
			int oldidx = i;
			destshp.poLayer->CreateField(srcshp.poLayer->GetLayerDefn()->GetFieldDefn(oldidx));
			fields2keep.push_back(srcshp.poLayer->GetLayerDefn()->GetFieldDefn(oldidx)->GetNameRef());
			newfields.push_back(oldidx);
			oldfields.push_back(oldidx);
		}
	}
	else
	{
		for (size_t i = 0; i < fields2keep.size(); i++)
		{
			int oldidx = srcshp.poLayer->GetLayerDefn()->GetFieldIndex(fields2keep[i].data());
			if (oldidx < 0)
			{
				oldidx = -1;
				newfields.push_back(-1);
			}
			else
			{
				destshp.poLayer->CreateField(srcshp.poLayer->GetLayerDefn()->GetFieldDefn(oldidx));
				newfields.push_back(destshp.poLayer->GetLayerDefn()->GetFieldIndex(fields2keep[i].data()));
			}
			oldfields.push_back(oldidx);
		}
	}

	srcshp.poLayer->ResetReading();
	OGRFeature *poFeatureOld;
	while ((poFeatureOld = srcshp.poLayer->GetNextFeature()) != NULL)
	{
		OGRFeature *poFeatureNew = OGRFeature::CreateFeature(destshp.poLayer->GetLayerDefn());
		poFeatureNew->SetGeometry(poFeatureOld->GetGeometryRef());
		for (size_t i = 0; i < fields2keep.size(); i++)
		{
			int oldfieldidx = oldfields[i];
			int newfieldidx = newfields[i];
			if (oldfieldidx == -1 || newfieldidx == -1)
				continue;

			poFeatureNew->SetField(newfieldidx,poFeatureOld->GetRawFieldRef(oldfieldidx));

		}

		destshp.poLayer->CreateFeature(poFeatureNew);
		OGRFeature::DestroyFeature(poFeatureNew);
		OGRFeature::DestroyFeature(poFeatureOld);
	}

}
void ShapeFile::copyDir(std::string srcDir, std::string destDir, std::vector<std::string> fields2keep)
{
	if (!QDir(srcDir.data()).exists())
		return;

	QDir qoutdir(destDir.data());
	if (!qoutdir.exists())
		qoutdir.mkpath(".");
	destDir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();
	std::vector<std::string> files;

	QDir input_dir(srcDir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	srcDir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = fileInfo.absoluteFilePath().toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		std::string output_file = destDir + fileInfo.fileName().toLocal8Bit().data();
		if (QFileInfo(output_file.data()).exists())
			continue;
		//if (!fileInfo.fileName().contains("NonRoad", Qt::CaseInsensitive))
		//	continue;
		//reprojectWithOGR(input_file, output_file, pythonFile);
		//if (isNULL(input_file.data()))
		//{
		//	remove(input_file);
		//	continue;
		//}
		copy(input_file, output_file, fields2keep);
		printf("%s\n", output_file.data());
	}

}

void ShapeFile::mergeFile(std::vector<std::string> inputfiles, std::vector<std::string> fields2keep, std::string outputfile)
{
	ShapeFile firstshp(inputfiles[0], 0);
	ShapeFile destshp;
	destshp.create(outputfile, firstshp.poLayer->GetSpatialRef(), 0, firstshp.poLayer->GetGeomType());
	std::vector<int> newfields;
	std::vector<int> oldfields;
	for (size_t i = 0; i < fields2keep.size(); i++)
	{
		int oldidx = firstshp.poLayer->GetLayerDefn()->GetFieldIndex(fields2keep[i].data());
		oldfields.push_back(oldidx);
		int newidx = -1;
		if (oldidx > -1)
		{
			destshp.poLayer->CreateField(firstshp.poLayer->GetLayerDefn()->GetFieldDefn(oldidx));
			newidx=destshp.poLayer->GetLayerDefn()->GetFieldIndex(fields2keep[i].data());
		}
		newfields.push_back(newidx);
	}
	firstshp.close();
	for (size_t nfile = 0; nfile < inputfiles.size(); nfile++)
	{
		ShapeFile srcshp(inputfiles[nfile], 0);
		OGRFeature *poFeatureOld;
		while ((poFeatureOld = srcshp.poLayer->GetNextFeature()) != NULL)
		{
			OGRFeature *poFeatureNew = OGRFeature::CreateFeature(destshp.poLayer->GetLayerDefn());
			poFeatureNew->SetGeometry(poFeatureOld->GetGeometryRef());
			for (size_t i = 0; i < oldfields.size(); i++)
			{
				int oldfieldidx = oldfields[i];
				int newfieldidx = newfields[i];
				if(oldfieldidx > -1)
				   poFeatureNew->SetField(newfieldidx, poFeatureOld->GetRawFieldRef(oldfieldidx));

			}
			destshp.poLayer->CreateFeature(poFeatureNew);
			OGRFeature::DestroyFeature(poFeatureNew);
			OGRFeature::DestroyFeature(poFeatureOld);
		}
		srcshp.close();
	}
	destshp.close();
}

void ShapeFile::copyDropGeometry(std::string src, std::string dest, std::vector<std::string> fields2keep)
{
	ShapeFile srcshp(src, 0);
	ShapeFile destshp;
	destshp.create(dest, srcshp.poLayer->GetSpatialRef(), 0, OGRwkbGeometryType::wkbPoint);

	std::vector<int> oldfields;
	std::vector<int> newfields;
	if (fields2keep.size() == 0)
	{
		for (size_t i = 0; i < srcshp.poLayer->GetLayerDefn()->GetFieldCount(); i++)
		{
			int oldidx = i;
			destshp.poLayer->CreateField(srcshp.poLayer->GetLayerDefn()->GetFieldDefn(oldidx));
			fields2keep.push_back(srcshp.poLayer->GetLayerDefn()->GetFieldDefn(oldidx)->GetNameRef());
			newfields.push_back(oldidx);
			oldfields.push_back(oldidx);
		}
	}
	else
	{
		for (size_t i = 0; i < fields2keep.size(); i++)
		{
			int oldidx = srcshp.poLayer->GetLayerDefn()->GetFieldIndex(fields2keep[i].data());
			if (oldidx < 0)
			{
				oldidx = -1;
				newfields.push_back(-1);
			}
			else
			{
				destshp.poLayer->CreateField(srcshp.poLayer->GetLayerDefn()->GetFieldDefn(oldidx));
				newfields.push_back(destshp.poLayer->GetLayerDefn()->GetFieldIndex(fields2keep[i].data()));
			}
			oldfields.push_back(oldidx);
		}
	}

	srcshp.poLayer->ResetReading();
	OGRFeature *poFeatureOld;
	while ((poFeatureOld = srcshp.poLayer->GetNextFeature()) != NULL)
	{
		OGRFeature *poFeatureNew = OGRFeature::CreateFeature(destshp.poLayer->GetLayerDefn());
		OGRPoint pt;
		OGRwkbGeometryType geomt = poFeatureOld->GetGeometryRef()->getGeometryType();
		if (geomt == OGRwkbGeometryType::wkbPoint || geomt == OGRwkbGeometryType::wkbPoint25D || geomt == OGRwkbGeometryType::wkbPointZM || geomt == OGRwkbGeometryType::wkbPointM)
		{
			OGRPoint* oldpt = (OGRPoint*)poFeatureOld->GetGeometryRef();
			pt.setX(oldpt->getX());
			pt.setY(oldpt->getY());
		}
		else
		{
			OGREnvelope env;
			poFeatureOld->GetGeometryRef()->getEnvelope(&env);
			pt.setX((env.MinX + env.MaxX)*0.5);
			pt.setY((env.MinY + env.MaxY)*0.5);
		}
		poFeatureNew->SetGeometry(&pt);
		for (size_t i = 0; i < fields2keep.size(); i++)
		{
			int oldfieldidx = oldfields[i];
			int newfieldidx = newfields[i];
			if (oldfieldidx == -1 || newfieldidx == -1)
				continue;

			poFeatureNew->SetField(newfieldidx, poFeatureOld->GetRawFieldRef(oldfieldidx));

		}

		destshp.poLayer->CreateFeature(poFeatureNew);
		OGRFeature::DestroyFeature(poFeatureNew);
		OGRFeature::DestroyFeature(poFeatureOld);
	}

}
void ShapeFile::copyDirDropGeometry(std::string srcDir, std::string destDir, std::vector<std::string> fields2keep)
{
	if (!QDir(srcDir.data()).exists())
		return;

	QDir qoutdir(destDir.data());
	if (!qoutdir.exists())
		qoutdir.mkpath(".");
	destDir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();
	std::vector<std::string> files;

	QDir input_dir(srcDir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	srcDir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = fileInfo.absoluteFilePath().toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		std::string output_file = destDir + fileInfo.fileName().toLocal8Bit().data();
		if (QFileInfo(output_file.data()).exists())
		{
			//remove(output_file);
			continue;
		}

		//if (!fileInfo.fileName().contains("NonRoad", Qt::CaseInsensitive))
		//	continue;
		//reprojectWithOGR(input_file, output_file, pythonFile);
		//if (isNULL(input_file.data()))
		//{
		//	remove(input_file);
		//	continue;
		//}
		copyDropGeometry(input_file, output_file, fields2keep);
		printf("%s\n", output_file.data());
	}

}
void scaleField(std::string shapefile, double scalefactor)
{
	ShapeFile shp(shapefile, 1);
	OGRFeature *poFeature;
	std::vector<int> fields;
	for (size_t i = 0; i < shp.poLayer->GetLayerDefn()->GetFieldCount(); i++)
	{
		std::string fieldname = shp.poLayer->GetLayerDefn()->GetFieldDefn(i)->GetNameRef();
		if (fieldname.substr(0, 2) == "ca")
			fields.push_back(shp.poLayer->GetLayerDefn()->GetFieldIndex(fieldname.data()));
	}
	shp.poLayer->ResetReading();
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		for (size_t i = 0; i < fields.size(); i++)
		{
			int idx = fields[i];
			double val = poFeature->GetFieldAsDouble(idx);
			poFeature->SetField(idx, val*scalefactor);
		}
		shp.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}

}


void ShapeFile::scaleCAFields(std::string srcDir)
{
	if (!QDir(srcDir.data()).exists())
		return;

	std::vector<std::string> files;

	QDir input_dir(srcDir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	srcDir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = fileInfo.absoluteFilePath().toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		double scalingfactor = 1000;
		scaleField(input_file, scalingfactor);
	}

}
bool ShapeFile::isNULL(std::string filename)
{
	OGRFeature *poFeature;
	ShapeFile shp(filename.data());
	shp.poLayer->ResetReading();
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		return false;
	}
	return true;
}

double ShapeFile::sum(std::string fieldname)
{
	double total = 0;
	OGRFeature *poFeature;
	int idx = poLayer->GetLayerDefn()->GetFieldIndex(fieldname.data());
	if (idx < 0)
		return total;

	poLayer->ResetReading();
	while ((poFeature = poLayer->GetNextFeature()) != NULL)
	{
		double value = poFeature->GetFieldAsDouble(idx);
		total += value;
		OGRFeature::DestroyFeature(poFeature);
	}
	return total;
}
std::vector<double> ShapeFile::sum(std::vector<std::string> fieldnames)
{
	std::vector<double> totals;
	for (size_t i = 0; i < fieldnames.size(); i++)
	{
		totals.push_back(sum(fieldnames[i]));
	}
	return totals;
}
void readText2Buf(std::string filename, char*& buf, int& length)
{
	std::ifstream is(filename.data(), std::ifstream::binary);
	is.seekg(0, is.end);
	length = is.tellg();
	is.seekg(0, is.beg);

	buf = new char[length];
	// read data as a block:
	is.read(buf, length);

	is.close();
}
void ShapeFile::reproject(std::string inputfile, std::string outputfile, std::string destwktfile)
{

	if (QFileInfo(outputfile.data()).exists())
		return;
	GDALDataset* poDS = (GDALDataset*)GDALOpenEx(inputfile.data(), GDAL_OF_VECTOR, NULL, NULL, NULL);
	OGRLayer* poLayer = poDS->GetLayer(0);
	int srcEPSG = 4326;
	char* buf = NULL;
	int length = 0;
	readText2Buf(destwktfile, buf, length);
	OGRSpatialReference oSRS;
	oSRS.importFromWkt(&buf);

	bool isWGS84 = false;
	std::string srcwktfile = QFileInfo(destwktfile.data()).absoluteDir().absolutePath().toLocal8Bit().data() + std::string("/srcwkt.txt");

	if (!poLayer->GetSpatialRef() || !poLayer->GetSpatialRef()->IsSame(&oSRS))
	{
		static char* srcwktbuf = new char[10000];
		if (!poLayer->GetSpatialRef())
		{
			OGRSpatialReference oSRS;
			oSRS.SetWellKnownGeogCS("WGS84");
			oSRS.exportToWkt(&srcwktbuf);
		}
		else
		{
			poLayer->GetSpatialRef()->exportToWkt(&srcwktbuf);
		}
		//QFileInfo fileinfo(inputfile.)
		std::string srcwktstr = srcwktbuf;
		std::ofstream ofs;
		ofs.open(srcwktfile);
		ofs << srcwktstr;
		ofs.close();
	}
	else if (poLayer->GetSpatialRef()->IsSame(&oSRS))
	{
		//copy(inputfile, outputfile);
		return;
	}
	//if (!poLayer->GetSpatialRef() || poLayer->GetSpatialRef()->GetEPSGGeogCS() == srcEPSG)
	//{
	//	//OGRSpatialReference oSRS;
	//	//oSRS.SetWellKnownGeogCS("WGS84");
	//	printf("%s: spatial reference not found.\n", inputfile.data());
	//	isWGS84 = true;
	//}
	//else if (poLayer->GetSpatialRef()->IsSame(&oSRS))
	//{
	//	copy(inputfile, outputfile);
	//	return;
	//}
	//else
	//{
	//	static char* srcwktbuf = new char[10000];
	//	poLayer->GetSpatialRef()->exportToWkt(&srcwktbuf);
	//	//QFileInfo fileinfo(inputfile.)
	//	std::string srcwktstr = srcwktbuf;
	//	std::ofstream ofs;
	//	ofs.open(srcwktfile);
	//	ofs << srcwktstr;
	//	ofs.close();
	//	//delete[] srcwktbuf;
	//}
	GDALClose(poDS);
	//static char* wktsrc = new char[100000];
	//static char* wktdest = new char[100000];

	//dSRS->exportToWkt(&wktdest);
	//oSRS.exportToWkt(&wktsrc);



	std::stringstream ss;
	if (isWGS84)
	{
		ss << "ogr2ogr -f" << " " << "\"" << "ESRI Shapefile" << "\"" << " " << outputfile << " " << inputfile
			<< " " << "-s_srs" << " " << "EPSG:" << srcEPSG
			<< " " << "-t_srs" << " " << destwktfile << "\n";
	}
	else
	{
		ss << "ogr2ogr -f" << " " << "\"" << "ESRI Shapefile" << "\"" << " " << outputfile << " " << inputfile
			<< " " << "-s_srs" << " " << srcwktfile
			<< " " << "-t_srs" << " " << destwktfile << "\n";

	}
	//std::stringstream ss;
	//ss << "ogr2ogr -f" << " " << "\"" << "ESRI Shapefile" << "\"" << " " << outputfile << " " << inputfile
	//	<< " " << "-s_srs" << " " << "EPSG:" << srcEPSG
	//	<< " " << "-t_srs" << " " << "EPSG:" << destEPSG << "\n";


	//"ogr2ogr -f "ESRI Shapefile" original.shp wgs84.shp -s_srs EPSG:27700 -t_srs EPSG:4326"
	printf("reprojecting: %s.\n", inputfile.data());
	std::string commandline = ss.str().data();
	system(commandline.data());
	QFile::remove(srcwktfile.data());
}
void ShapeFile::reprojectDir(std::string indir, std::string outdir, std::string destwktfile)
{
	if (!QDir(indir.data()).exists())
		return;

	QDir qoutdir(outdir.data());
	if (!qoutdir.exists())
		qoutdir.mkpath(".");
	outdir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();
	std::vector<std::string> files;

	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = fileInfo.absoluteFilePath().toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		std::string output_file = outdir + fileInfo.fileName().toLocal8Bit().data();
		//reprojectWithOGR(input_file, output_file, pythonFile);
		reproject(input_file, output_file, destwktfile);
		printf("%s\n", output_file.data());
	}

}
ShapeFile::~ShapeFile()
{
	close();
}
