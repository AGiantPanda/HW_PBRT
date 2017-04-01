
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "camera.h"
#include "paramset.h"
#include "film.h"

struct LenSpec {
	float radius, radius2;
	float oz;
	float n, n2;
	float aperture;
	float aper2;
	bool isConvex;// convex need to choose the outer inter point, concave need to choose the inner one
};

// RealisticCamera Declarations
class RealisticCamera : public Camera {
public:
	// RealisticCamera Public Methods
	RealisticCamera(const AnimatedTransform &cam2world,
						float hither, float yon, float sopen,
						float sclose, float filmdistance, float aperture_diameter, string specfile,
						float filmdiag, Film *film);
	~RealisticCamera();
	float GenerateRay(const CameraSample &sample, Ray *) const;
  
protected:
	// RealisticCamera Protected Data
	Transform CameraToScreen, RasterToCamera;
	Transform ScreenToRaster, RasterToScreen;
	float lensRadius, focalDistance;

private:
	// RealisticCamera Private Methods
	bool ParseCameraSpec(string specfile);

	// RealisticCamera Private Parameters
	float filmdiag;
	float totaldistance;
	LenSpec *lens;
	int lenNum;
};


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);


#endif	// PBRT_CAMERAS_REALISTIC_H