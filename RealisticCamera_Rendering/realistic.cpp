#include "stdafx.h"
#include <fstream>
#include <iostream>
#include "cameras/realistic.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
				 float hither, float yon, 
				 float sopen, float sclose, 
				 float filmdistance, float aperture_diameter, string specfile, 
				 float filmdiag, Film *f)
	: Camera(cam2world, sopen, sclose, f) // pbrt-v2 doesnot specify hither and yon
{
	this->filmdiag = filmdiag;
	focalDistance = filmdistance;
	
	// Read the lens' spec
	// lens = nullptr;
	if (!ParseCameraSpec(specfile))
	{
		std::cout << "Cannot open " << specfile << "! Plz check it." << std::endl;
		exit(1);
	}

	// To create Screen/Camera Space
	float diagDistance = sqrtf(film->xResolution * film->xResolution + 
		film->yResolution*film->yResolution);
	float ratiox = film->xResolution / diagDistance;
	float ratioy = film->yResolution / diagDistance;
	float screen[4];
	screen[0] = float(-filmdiag * ratiox * 0.5);
	screen[1] = -screen[0];
	screen[2] = float(-filmdiag * ratioy * 0.5);
	screen[3] = -screen[2];

	// Compute projective camera transformations
	// camera would be set in (0, 0, 0), the outer len would on the (0, 0, totaldistance)
	// the raster space / image plane would be on the O pos
	CameraToScreen = Scale(1.f, 1.f, 1.f / (yon - hither)) *
		Translate(Vector(0.f, 0.f, 0.f));

	// Compute Projective camera screen transformations
	ScreenToRaster = Scale(float(film->xResolution),
		                   float(film->yResolution), 1.f) *
		Scale(1.f / (screen[1] - screen[0]),
			1.f / (screen[2] - screen[3]), 1.f) *
		Translate(Vector(-screen[0], -screen[3], 0.f));
	RasterToScreen = Scale(-1.f, -1.f, 1.f) * Inverse(ScreenToRaster);
	RasterToCamera = Inverse(CameraToScreen) * RasterToScreen;
}

RealisticCamera::~RealisticCamera()
{
	delete[] this->lens;
}

// Takes a sample position in image space (given by sample.imageX and sampleY) as
// an argument and should return a **random** ray into the scene.
// Also return value is used to calculate the weight of the ray.
// The number of sample gives you exactly that number of samples to a single pixel
float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
	float INFINITE = 3.4e+38;
	Point Pras(sample.imageX, sample.imageY, 0.f);
	Point Pcamera;
	RasterToCamera(Pras, &Pcamera);
	float weight = 0.f;

	// a Monte Carlo Function to generate 2D points, page 667
	float lensU, lensV;
	ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
	lensU *= lensRadius;
	lensV *= lensRadius;
	Vector dir = Normalize(Vector(lensU - Pcamera.x, lensV - Pcamera.y, focalDistance));
	Ray incident(Pcamera, dir, 0.f, INFINITE);

	// count from the inner point
	for (int i = lenNum - 1; i >= 0; i--) {
		// get the intersect point
		Point interP;
		if (lens[i].radius == 0.f)
		{
			// this len is air, just need to check the aperture
			float t = (lens[i].oz - incident.o.z) / incident.d.z;
			interP.x = incident.o.x + t * incident.d.x;
			interP.y = incident.o.y + t * incident.d.y;
			if (interP.x * interP.x + interP.y * interP.y > lens[i].aper2) {
				// The ray is blocked, give it an invalid value and return
				*ray = Ray(Pcamera, Vector(0, 0, 1), INFINITE, 0.f);
				return 0.f;
			}
		}
		else
		{
			// to calc the intersect point, i really hope it is right
			Vector OC = incident.o - Point(0, 0, lens[i].oz);
			float a = incident.d.LengthSquared();
			float b = 2 * Dot(OC, incident.d);
			float c = OC.LengthSquared() -lens[i].radius2;
			float d = b * b - 4 * a * c;
			if (d < 0) {
				*ray = Ray(Pcamera, Vector(0, 0, 1), INFINITE, 0.f);
				return 0.f;
			}
			else
				d = sqrt(d);

			// calc the interPoint and refraction ray with Heckber's method
			Vector Normal;
			Vector I = Normalize(incident.d);
			if (lens[i].isConvex)
			{
				// Convex, get the former point of the intersect
				float t = (-b + d) / a * 0.5;
				interP.x = incident.o.x + t * incident.d.x;
				interP.y = incident.o.y + t * incident.d.y;
				interP.z = incident.o.z + t * incident.d.z;
				Normal = Normalize(Point(0, 0, lens[i].oz) - interP);
			}
			else
			{
				// concave, get the inner point
				float t = (-b - d) / a * 0.5;
				interP.x = incident.o.x + t * incident.d.x;
				interP.y = incident.o.y + t * incident.d.y;
				interP.z = incident.o.z + t * incident.d.z;
				Normal = Normalize(interP - Point(0, 0, lens[i].oz));
			}

			if (interP.x * interP.x + interP.y * interP.y > lens[i].aper2) {
				// The ray is blocked, give it an invalid value and return
				*ray = Ray(Pcamera, Vector(0, 0, 1), INFINITE, 0.f);
				return 0.f;
			}

			float c1 = Dot(-I, Normal);
			float c2 = 1 - lens[i].n2 * (1 - c1 * c1);
			if (c2 < 0)
			{
				*ray = Ray(Pcamera, Vector(0, 0, 1), INFINITE, 0.f);
				return 0.f;
			}
			else
				c2 = sqrt(c2);
			incident.d = Normalize(lens[i].n * I + (lens[i].n * c1 - c2) * Normal);
			incident.o = interP;
		}

	}
	
	// calc the weight, do it when end tracing the lens
	//Vector v1 = Normalize(incident.o - Pcamera);
	Vector v1 = Normalize(Vector(Pcamera.x, Pcamera.y, Pcamera.z - lens[lenNum - 1].oz));
	Vector v2 = Vector(0, 0, -1);
	float cosT = Dot(v1, v2);
	cosT = cosT * cosT * cosT * cosT;
	weight = lens[lenNum - 1].aper2 * M_PI * cosT / (focalDistance * focalDistance);
	Clamp(weight, 0.0, 1.0);
	*ray = incident;
	ray->time = sample.time;
	CameraToWorld(*ray, ray);
	return weight;
}

bool RealisticCamera::ParseCameraSpec(string specfile) {
	totaldistance = focalDistance;
	std::ifstream cameraSpec;
	cameraSpec.open(specfile.c_str());
	if (!cameraSpec.is_open())
		return 0;

	vector<float> rad, axp, dn, ape;
	string line;
	while (std::getline(cameraSpec, line))
	{
		if (line[0] != '#')
		{
			float r, a, n, p;
			sscanf(line.c_str(), "%f\t%f\t%f\t%f", &r, &a, &n, &p);
			rad.push_back(r);
			axp.push_back(a);
			if(r == 0)
				dn.push_back(1);
			else
				dn.push_back(n);
			ape.push_back(p);
		}
	}

	cameraSpec.close();

	lenNum = rad.size();
	lens = new LenSpec[lenNum];
	float total_r = focalDistance;
	
	for (int i = lenNum - 1; i >= 0; i--)
	{
		// len's radius & square(radius)
		lens[i].radius = abs(rad[i]);
		lens[i].radius2 = rad[i] * rad[i];
		// len's totaldistance from camera position (0, 0, 0)
		totaldistance += axp[i];
		// len's z position from (0, 0, 0), 
		// for a convex it should be len.pos - len.radius,
		// for a concave it should be len.pos + len.rad
		lens[i].oz = totaldistance - rad[i];
		// the Eta value between two different materials, the air would be 1
		float right = dn[i];
		float left = (i == 0) ? 1.f : dn[i - 1];
		lens[i].n = right / left;
		lens[i].n2 = lens[i].n * lens[i].n;
		// len.aperture & its square value
		lens[i].aperture = ape[i];
		lens[i].aper2 = ape[i] * ape[i] * 0.25;
		lens[i].isConvex = rad[i] > 0 ? true : false;
	}
	lensRadius = lens[lenNum - 1].aperture * 0.5;

	return 1;
}

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	// Extract common camera parameters from \use{ParamSet}
	float hither = params.FindOneFloat("hither", -1);
	float yon = params.FindOneFloat("yon", -1);
	float shutteropen = params.FindOneFloat("shutteropen", -1);
	float shutterclose = params.FindOneFloat("shutterclose", -1);

	// Realistic camera-specific parameters
	string specfile = params.FindOneString("specfile", "");
	float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film, if found return filmdistance, else return 70.0
 	float fstop = params.FindOneFloat("aperture_diameter", 1.0);	
	float filmdiag = params.FindOneFloat("filmdiag", 35.0);

	Assert(hither != -1 && yon != -1 && shutteropen != -1 &&
		shutterclose != -1 && filmdistance!= -1);
	if (specfile == "") {
	    Severe( "No lens spec file supplied!\n" );
	}
	return new RealisticCamera(cam2world, hither, yon,
				   shutteropen, shutterclose, filmdistance, fstop, 
				   specfile, filmdiag, film);
}
