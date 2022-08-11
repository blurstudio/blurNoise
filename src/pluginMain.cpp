#include <maya/MFnPlugin.h>
#include <maya/MGlobal.h>
#include "blurNoise.h"

MStatus initializePlugin(MObject obj) {
	MStatus result;
	MFnPlugin plugin(obj, "Blur Studio", "1.0", "Any");
	result = plugin.registerNode(
        DEFORMER_NAME,
        blurNoise::id,
        blurNoise::creator,
        blurNoise::initialize,
        MPxNode::kDeformerNode
    );
	return result;
}

MStatus uninitializePlugin(MObject obj) {
	MStatus result;
	MFnPlugin plugin(obj);
	result = plugin.deregisterNode(blurNoise::id);
	return result;
}

