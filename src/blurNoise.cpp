#include <maya/MTypeId.h> 

#include <maya/MTime.h> 
#include <maya/MMatrix.h>
#include <maya/MDoubleArray.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MFloatVectorArray.h>
#include <maya/MTransformationMatrix.h>

#include <maya/MItGeometry.h>
#include <maya/MItMeshVertex.h>
#include <maya/MItSurfaceCV.h>

#include <maya/MFnNurbsSurface.h>
#include <maya/MFnMesh.h>

#include <maya/MFnUnitAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnEnumAttribute.h>

#include <vector>

#include "blurNoise.h"


MObject blurNoise::aAmp;
MObject blurNoise::aAmpOffset;
MObject blurNoise::aAmpClampHigh;
MObject blurNoise::aAmpClampLow;
MObject blurNoise::aAmpUseClampHigh;
MObject blurNoise::aAmpUseClampLow;

MObject blurNoise::aTime;

MObject blurNoise::aTranslate;
    MObject blurNoise::aTranslateX;
    MObject blurNoise::aTranslateY;
    MObject blurNoise::aTranslateZ;
MObject blurNoise::aRotate;
    MObject blurNoise::aRotateX;
    MObject blurNoise::aRotateY;
    MObject blurNoise::aRotateZ;
MObject blurNoise::aScale;
    MObject blurNoise::aScaleX;
    MObject blurNoise::aScaleY;
    MObject blurNoise::aScaleZ;
MObject blurNoise::aShear;
    MObject blurNoise::aShearXY;
    MObject blurNoise::aShearXZ;
    MObject blurNoise::aShearYZ;
MObject blurNoise::aRotateOrder;

MObject blurNoise::aOctaveCount;
MObject blurNoise::aOctaveScaleBase;
MObject blurNoise::aOctaveRangeBase;



const MTypeId blurNoise::id(0x00122712);
void* blurNoise::creator() {return new blurNoise();}
MStatus blurNoise::initialize() {
    
    MStatus stat;
    MFnNumericAttribute fnNum;
    MFnUnitAttribute fnUnit;
    MFnEnumAttribute fnEnum;

    aAmp = fnNum.create("amp", "amp", MFnNumericData::kDouble, 1.0);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    stat = addAttribute(aAmp);

    aAmpOffset = fnNum.create("ampOffset", "ampo", MFnNumericData::kDouble, 0.0);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    stat = addAttribute(aAmpOffset);

    aAmpClampHigh = fnNum.create("ampClampHigh", "ach", MFnNumericData::kDouble, 1.0);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    stat = addAttribute(aAmpClampHigh);

    aAmpClampLow = fnNum.create("ampClampLow", "acl", MFnNumericData::kDouble, -1.0);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    stat = addAttribute(aAmpClampLow);

    aAmpUseClampHigh = fnNum.create("ampUseClampHigh", "auh", MFnNumericData::kBoolean, false);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    stat = addAttribute(aAmpUseClampHigh);

    aAmpUseClampLow = fnNum.create("ampUseClampLow", "aul", MFnNumericData::kBoolean, false);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    stat = addAttribute(aAmpUseClampLow);

    aTime = fnUnit.create("time", "time", MFnUnitAttribute::kTime);
    fnUnit.setStorable(true);
    fnUnit.setKeyable(true);
    stat = addAttribute(aTime);


    aTranslateX = fnUnit.create("translateX", "tx", MFnUnitAttribute::kDistance, 0.0);
    fnUnit.setStorable(true);
    fnUnit.setKeyable(true);
    aTranslateY = fnUnit.create("translateY", "ty", MFnUnitAttribute::kDistance, 0.0);
    fnUnit.setStorable(true);
    fnUnit.setKeyable(true);
    aTranslateZ = fnUnit.create("translateZ", "tz", MFnUnitAttribute::kDistance, 0.0);
    fnUnit.setStorable(true);
    fnUnit.setKeyable(true);
    aTranslate = fnNum.create("translate", "t", aTranslateX, aTranslateY, aTranslateZ);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    stat = addAttribute(aTranslate);


    aRotateX = fnUnit.create("rotateX", "rx", MFnUnitAttribute::kAngle, 0.0);
    fnUnit.setStorable(true);
    fnUnit.setKeyable(true);
    aRotateY = fnUnit.create("rotateY", "ry", MFnUnitAttribute::kAngle, 0.0);
    fnUnit.setStorable(true);
    fnUnit.setKeyable(true);
    aRotateZ = fnUnit.create("rotateZ", "rz", MFnUnitAttribute::kAngle, 0.0);
    fnUnit.setStorable(true);
    fnUnit.setKeyable(true);
    aRotate = fnNum.create("rotate", "r", aRotateX, aRotateY, aRotateZ);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    stat = addAttribute(aRotate);


    aScaleX = fnNum.create("scaleX", "sx", MFnNumericData::kDouble, 1.0);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    aScaleY = fnNum.create("scaleY", "sy", MFnNumericData::kDouble, 1.0);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    aScaleZ = fnNum.create("scaleZ", "sz", MFnNumericData::kDouble, 1.0);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    aScale = fnNum.create("scale", "s", aScaleX, aScaleY, aScaleZ);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    stat = addAttribute(aScale);


    aShearXY = fnNum.create("shearXY", "shxy", MFnNumericData::kDouble, 0.0);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    aShearXZ = fnNum.create("shearXZ", "shxz", MFnNumericData::kDouble, 0.0);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    aShearYZ = fnNum.create("shearYZ", "shyz", MFnNumericData::kDouble, 0.0);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    aShear = fnNum.create("shear", "sh", aShearXY, aShearXZ, aShearYZ);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    stat = addAttribute(aShear);


    aRotateOrder = fnEnum.create("rotateOrder", "ro", 0, &stat);
    fnEnum.setKeyable(true);
    fnEnum.addField("xyz", 0);
    fnEnum.addField("yzx", 1);
    fnEnum.addField("zxy", 2);
    fnEnum.addField("xzy", 3);
    fnEnum.addField("yxz", 4);
    fnEnum.addField("zyx", 5);
    stat = addAttribute(aRotateOrder);


    aOctaveCount = fnNum.create("roughnessLoops", "rl", MFnNumericData::kInt, 1);
    fnNum.setMin(1);
    fnNum.setMax(20);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    stat = addAttribute(aOctaveCount);

    aOctaveScaleBase = fnNum.create("roughnessScale", "rs", MFnNumericData::kDouble, 0.25);
    fnNum.setMin(1e-10);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    stat = addAttribute(aOctaveScaleBase);

    aOctaveRangeBase = fnNum.create("roughnessAmplitude", "ra", MFnNumericData::kDouble, 0.25);
    fnNum.setMin(1e-10);
    fnNum.setStorable(true);
    fnNum.setKeyable(true);
    stat = addAttribute(aOctaveRangeBase);

    std::vector<MObject *> iobjs = {
        &aAmp, &aAmpOffset, &aAmpClampHigh, &aAmpClampLow,
        &aAmpUseClampHigh, &aAmpUseClampLow, &aTime,
        &aTranslate, &aTranslateX, &aTranslateY, &aTranslateZ,
        &aRotate, &aRotateX, &aRotateY, &aRotateZ,
        &aScale, &aScaleX, &aScaleY, &aScaleZ,
        &aShear, &aShearXY, &aShearXZ, &aShearYZ,
        &aRotateOrder, &aOctaveCount,
        &aOctaveScaleBase, &aOctaveRangeBase
    };

    std::vector<MObject *> oobjs = {&outputGeom};

    for (auto &ii : iobjs) {
        for (auto &oo : oobjs) {
            attributeAffects(*ii, *oo);
        }
    }

    return MStatus::kSuccess;
}


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
	const MMatrix& wmat,
	unsigned int multiIndex
){
	MStatus status;

	MDataHandle hEnv = dataBlock.inputValue(envelope);
	float env = hEnv.asFloat();
    if (env == 0.0) {
        return MStatus::kSuccess;
    }

    // Get the time input
    MDataHandle hTime = dataBlock.inputValue(aTime);
    double time = hTime.asTime().asUnits(MTime::kSeconds);


    // Get the octave loop data
    MDataHandle hOctaveCount = dataBlock.inputValue(aOctaveCount);
    MDataHandle hOctaveScaleBase = dataBlock.inputValue(aOctaveScaleBase);
    MDataHandle hOctaveRangeBase = dataBlock.inputValue(aOctaveRangeBase);
    int octaveCount = hOctaveCount.asInt();
    double octaveScaleBase = 1.0 / hOctaveScaleBase.asDouble();
    double octaveRangeBase = 1.0 / hOctaveRangeBase.asDouble();


    // Get the amplitude data
    MDataHandle hAmp = dataBlock.inputValue(aAmp);
    MDataHandle hAmpOffset = dataBlock.inputValue(aAmpOffset);
    MDataHandle hAmpClampHigh = dataBlock.inputValue(aAmpClampHigh);
    MDataHandle hAmpClampLow = dataBlock.inputValue(aAmpClampLow);
    MDataHandle hAmpUseClampHigh = dataBlock.inputValue(aAmpUseClampHigh);
    MDataHandle hAmpUseClampLow = dataBlock.inputValue(aAmpUseClampLow);
    double amp = hAmp.asDouble();
    double ampOffset = hAmpOffset.asDouble();
    double ampClampHigh = hAmpClampHigh.asDouble();
    double ampClampLow = hAmpClampLow.asDouble();
    bool ampUseClampHigh = hAmpUseClampHigh.asBool();
    bool ampUseClampLow = hAmpUseClampLow.asBool();


    // Get the SRTSh inputs and compose them into a matrix
    // I should probably set the shear to be hidden
    MDataHandle hTranslate = dataBlock.inputValue(aTranslate);
    MDataHandle hRotate = dataBlock.inputValue(aRotate);
    MDataHandle hScale = dataBlock.inputValue(aScale);
    MDataHandle hShear = dataBlock.inputValue(aShear);
    MDataHandle hRotateOrder = dataBlock.inputValue(aRotateOrder);
    MVector translate = hTranslate.asVector();
    double *rotate = hRotate.asDouble3();
    double *scale = hScale.asDouble3();
    double *shear = hShear.asDouble3();
    short rotateOrder = hRotateOrder.asShort();
    MTransformationMatrix::RotationOrder rotationOrders[6] = {
        MTransformationMatrix::kXYZ, MTransformationMatrix::kYZX,
        MTransformationMatrix::kZXY, MTransformationMatrix::kXZY,
        MTransformationMatrix::kYXZ, MTransformationMatrix::kZYX
    };
    MTransformationMatrix mtm;
    mtm.setScale(scale, MSpace::kWorld);
    mtm.setRotation(rotate, rotationOrders[rotateOrder]);
    mtm.setTranslation(translate, MSpace::kWorld);
    mtm.setShear(shear, MSpace::kWorld);
    MMatrix mat = mtm.asMatrixInverse();


	// Maya automatically copies the input plug to the output plug
	// and then gives you an iterator over that
	// So get the OUTPUT handle for this mutiIndex
	MPlug outPlug(thisMObject(), outputGeom);
	outPlug.selectAncestorLogicalIndex(multiIndex, outputGeom);
	MDataHandle hOutput = dataBlock.outputValue(outPlug);



    auto outType = hOutput.type();
    if (outType == MFnData::kMesh) {

        MObject mesh = hOutput.asMesh();
        MFnMesh meshFn(mesh);
        MFloatVectorArray norms;
        meshFn.getVertexNormals(false, norms);


        MItMeshVertex vertIter(mesh);
        MPointArray newPts;
        newPts.setLength(vertIter.count());
        int pIdx = 0;
        for (; !vertIter.isDone(); vertIter.next()) {
            int idx = vertIter.index();
            float weight =  weightValue(dataBlock, multiIndex, idx);
            if (weight == 0.0) {
                continue;
            }

            MPoint pos = vertIter.position(MSpace::kWorld, &status) * mat;
            double offset = 0.0;
            double scaleBase = 1.0;
            double rangeBase = 1.0;
            for (int oIdx=0; oIdx<octaveCount; ++oIdx){
                //double n = noise4_Classic(ose, osg, pos.x*scaleBase, pos.y*scaleBase, pos.z*scaleBase, time);
                double n = noise4_XYZBeforeW(ose, osg, pos.x*scaleBase, pos.y*scaleBase, pos.z*scaleBase, time);
                offset += n / rangeBase;
                scaleBase *= octaveScaleBase;
                rangeBase *= octaveRangeBase;
            }
            offset = (offset * amp) + ampOffset;
            if (ampUseClampHigh && offset > ampClampHigh)
                offset = ampClampHigh;
            if (ampUseClampLow && offset < ampClampLow)
                offset = ampClampLow;

            //MVector norm;
            //vertIter.getNormal(norm, MSpace::kObject);
            newPts[pIdx++] = vertIter.position() + (norms[idx] * weight * offset * env);
        }
        pIdx = 0;
        vertIter.reset();
        for (; !vertIter.isDone(); vertIter.next()) {
            vertIter.setPosition(newPts[pIdx++]);
        }

    }
    else if (outType == MFnData::kNurbsSurface) {

        MObject surface = hOutput.asNurbsSurface();
        MItSurfaceCV cvIter(surface);
        MFnNurbsSurface fnSurf(surface);
        MDoubleArray uParams, vParams;
        MPointArray newPts;

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



                MPoint pos = cvIter.position(MSpace::kWorld, &status) * mat;
                double offset = 0.0;
                double scaleBase = 1.0;
                double rangeBase = 1.0;
                for (int oIdx=0; oIdx<octaveCount; ++oIdx){
                    double n = noise4_Classic(ose, osg, pos.x*scaleBase, pos.y*scaleBase, pos.z*scaleBase, time);
                    offset += n / rangeBase;
                    scaleBase *= octaveScaleBase;
                    rangeBase *= octaveRangeBase;
                }
                offset = (offset * amp) + ampOffset;
                if (ampUseClampHigh && offset > ampClampHigh)
                    offset = ampClampHigh;
                if (ampUseClampLow && offset < ampClampLow)
                    offset = ampClampLow;



                double uParam = uParams[idxU];
                double vParam = uParams[idxV];
                MVector norm = fnSurf.normal(uParam, vParam, MSpace::kObject);
                newPts.append(cvIter.position() + (norm * weight * offset * env));
            }
        }
        cvIter.reset();
        int pIdx = 0;
        for (; !cvIter.isDone(); cvIter.nextRow()) {
            for (; !cvIter.isRowDone(); cvIter.next()) {
                cvIter.setPosition(newPts[pIdx++]);
            }
        }
        cvIter.updateSurface();
    }
    return MStatus::kSuccess;
}
