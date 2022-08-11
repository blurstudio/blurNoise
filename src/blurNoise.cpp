#include <maya/MFnPlugin.h>
#include <maya/MTypeId.h> 

#include <maya/MMatrix.h>
#include <maya/MMatrixArray.h>
#include <maya/MStringArray.h>
#include <maya/MFloatArray.h>
#include <maya/MDoubleArray.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>

#include <maya/MPxDeformerNode.h> 
#include <maya/MFnMesh.h> 
#include <maya/MItGeometry.h>
#include <maya/MFnPointArrayData.h>
#include <maya/MFnComponentListData.h>
#include <maya/MItMeshVertex.h>
#include <maya/MItSurfaceCV.h>
#include <maya/MFnNurbsSurface.h>
#include <maya/MWeight.h>

#include <vector>

#include "blurNoise.h"
#include "OpenSimplex2F.h"

const MTypeId blurNoise::id(0x00122706);
void* blurNoise::creator() {return new blurNoise();}
MStatus blurNoise::initialize() {return MStatus::kSuccess;}


/*
Get the parameters where each CV has its most influence along the U and V parameters
We will get the normals to move along for these parameters
*/
void blurNoise::getSurfaceCVParams(const MFnNurbsSurface &fnSurf, MDoubleArray &uParams, MDoubleArray &vParams) const {
    auto uCount = fnSurf.numCVsInU();
    auto uForm = fnSurf.formInU();
    auto uDegree = fnSurf.degreeU();
    MDoubleArray uKnots;
    fnSurf.getKnotsInU(uKnots);
    int iStart = 0;
    int iEnd = uCount;
    if (uForm == MFnNurbsSurface::Form::kPeriodic){
        iStart = 1;
        iEnd -= uDegree - 1;
    }
    for (int i = iStart; i < iEnd; ++i) {
        double summer = 0;
        for (int j = 0; j < uDegree; ++j) {
            summer += uKnots[i + j];
        }
        summer /= uDegree;
        uParams.append(summer);
    }

    auto vCount = fnSurf.numCVsInV();
    auto vForm = fnSurf.formInV();
    auto vDegree = fnSurf.degreeV();
    MDoubleArray vKnots;
    fnSurf.getKnotsInV(vKnots);
    iStart = 0;
    iEnd = vCount;
    if (vForm == MFnNurbsSurface::Form::kPeriodic) {
        iStart = 1;
        iEnd -= vDegree - 1;
    }
    for (int i = iStart; i < iEnd; ++i) {
        double summer = 0;
        for (int j = 0; j < vDegree; ++j) {
            summer += vKnots[i + j];
        }
        summer /= vDegree;
        vParams.append(summer);
    }
}

blurNoise::blurNoise() {
    ose = initOpenSimplex();
    osg = newOpenSimplexGradients(ose, 1234);
}
blurNoise::~blurNoise() {
    free(ose);
    free(osg);
}

MStatus blurNoise::deform(
	MDataBlock& dataBlock,
	MItGeometry& geoIter,
	const MMatrix& mat,
	unsigned int multiIndex
){
	MStatus status;

	MDataHandle hEnv = dataBlock.inputValue(envelope);
	float env = hEnv.asFloat();
	double time = 0;

	// Maya automatically copies the input plug to the output plug
	// and then gives you an iterator over that
	// So get the OUTPUT handle for this mutiIndex
	MPlug outPlug(thisMObject(), outputGeom);
	outPlug.selectAncestorLogicalIndex(multiIndex, outputGeom);
	MDataHandle hOutput = dataBlock.outputValue(outPlug);


    auto outType = hOutput.type();
    if (outType == MFnData::kMesh) {
        MObject mesh = hOutput.asMesh();
        MItMeshVertex vertIter(mesh);
        for (; !vertIter.isDone(); vertIter.next()) {
            int idx = vertIter.index();
            float weight =  weightValue(dataBlock, multiIndex, idx);
            if (weight == 0.0) {
                continue;
            }

            MPoint pos = vertIter.position(MSpace::kWorld, &status);
            double offset = noise4_Classic(ose, osg, pos.x, pos.y, pos.z, time);
            MVector norm;
            vertIter.getNormal(norm, MSpace::kObject);
            vertIter.setPosition(vertIter.position() + (norm * weight * offset));
        }
    }
    else if (outType == MFnData::kNurbsSurface) {

        MObject surface = hOutput.asNurbsSurface();
        MItSurfaceCV cvIter(surface);
        MFnNurbsSurface fnSurf(surface);
        MDoubleArray uParams, vParams;
        getSurfaceCVParams(fnSurf, uParams, vParams);
        for (; !cvIter.isDone(); cvIter.nextRow()) {
            for (; !cvIter.isRowDone(); cvIter.next()) {
                int idxU, idxV;
                int idx = cvIter.index();
                cvIter.getIndex(idxU, idxV);
                float weight =  weightValue(dataBlock, multiIndex, idx);
                if (weight == 0.0) {
                    continue;
                }

                MPoint pos = cvIter.position(MSpace::kWorld, &status);
                double offset = noise4_Classic(ose, osg, pos.x, pos.y, pos.z, time);
                double uParam = uParams[idxU];
                double vParam = uParams[idxV];
                MVector norm = fnSurf.normal(uParam, vParam, MSpace::kObject);
                cvIter.setPosition(cvIter.position() + (norm * weight * offset));
            }
        }
        cvIter.updateSurface();

    }
    return MStatus::kSuccess;
}
