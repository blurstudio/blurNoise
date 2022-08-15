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
    float time = (float)hTime.asTime().asUnits(MTime::kSeconds);


    // Get the octave loop data
    MDataHandle hOctaveCount = dataBlock.inputValue(aOctaveCount);
    MDataHandle hOctaveScaleBase = dataBlock.inputValue(aOctaveScaleBase);
    MDataHandle hOctaveRangeBase = dataBlock.inputValue(aOctaveRangeBase);
    int octaveCount = hOctaveCount.asInt();
    double octaveScaleBase = 1.0 / hOctaveScaleBase.asDouble();
    double octaveRangeBase = 1.0 / hOctaveRangeBase.asDouble();
    //if (octaveCount > 1) {fnFractal->SetOctaveCount(octaveCount); }

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


    MPointArray worldPts, localPts, newPts;
    geoIter.allPositions(worldPts, MSpace::kWorld);
    geoIter.allPositions(localPts, MSpace::kObject);
    newPts.setLength(worldPts.length());

    // store the indices of all the geo we're looping over
    std::vector<int> allIdxs(geoIter.count());
    std::vector<float> allWeights(geoIter.count());
    for (; !geoIter.isDone(); geoIter.next()) {
        int pi = geoIter.positionIndex();
        int idx = geoIter.index();
        allIdxs[pi] = idx;
        allWeights[pi] = weightValue(dataBlock, multiIndex, idx);
    }
    geoIter.reset();


    // Get the per-point normals in object space
    MFloatVectorArray allNorms;
    allNorms.setLength(geoIter.count());
    auto outType = hOutput.type();
    if (outType == MFnData::kMesh) {
        MObject mesh = hOutput.asMesh();
        MFnMesh meshFn(mesh);
        MFloatVectorArray tnorms;
        meshFn.getVertexNormals(false, tnorms);
        // match the norms array to the allpts array
        int i = 0;
        for (const int &idx: allIdxs){
            allNorms[i++] = tnorms[idx];
        }
    }
    else if (outType == MFnData::kNurbsSurface) {
        MObject surface = hOutput.asNurbsSurface();
        MFnNurbsSurface fnSurf(surface);
        MDoubleArray uParams, vParams;
        getSurfaceCVParams(fnSurf, uParams, vParams);
        int i = 0;
        // It's possible to do this without the tnorms
        // but this works well enough in practice
        MFloatVectorArray tnorms;
        for (const double &u: uParams){
            for (const double &v: vParams){
                tnorms[i++] = fnSurf.normal(u, v, MSpace::kObject);
            }
        }
        int idx = 0;
        for (const int &idx: allIdxs){
            allNorms[i++] = tnorms[idx];
        }
    }
    else {
        // Only mesh and nurbs supported
        return MStatus::kFailure;
    }

    // Copy the world data from the array-of-arrays mpoint array 
    // into the multiple flat vectors, then run all the noise
    // at once in parallel

    std::vector<float> noiseVals(worldPts.length());
    float scaleBase = 1.0;
    float rangeBase = 1.0;
    for (int octave = 0; octave < octaveCount; ++octave) {

        #pragma omp parallel for
        for (int i=0; i<(int)worldPts.length(); ++i){
            MPoint pos = worldPts[i] * mat;
            float n = (float)noise4_XYZBeforeW(ose, osg, pos.x*scaleBase, pos.y*scaleBase, pos.z*scaleBase, time);
            noiseVals[i] += n / rangeBase;
        }

        scaleBase *= (float)octaveScaleBase;
        rangeBase *= (float)octaveRangeBase;
    }


    // Move the points based off the noise value
    //# pragma omp parallel for
    for (int idx=0; idx<allIdxs.size(); ++idx){
        float weight = allWeights[idx];
        if (weight == 0.0) {
            newPts[idx] = localPts[idx];
            continue;
        }

        double n = ((double)noiseVals[idx] * amp) + ampOffset;
        if (ampUseClampHigh && n > ampClampHigh)
            n = ampClampHigh;
        if (ampUseClampLow && n < ampClampLow)
            n = ampClampLow;

        newPts[idx] = localPts[idx] + (allNorms[idx] * weight * n * env);           
    }

    geoIter.setAllPositions(newPts, MSpace::kObject);
    return MStatus::kSuccess;
}
