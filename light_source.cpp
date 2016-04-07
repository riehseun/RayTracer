/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

void PointLight::shade( Ray3D& ray ) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  
	/*
	Colour ra = ray.intersection.mat->ambient;
	Colour rd = ray.intersection.mat->diffuse;
	Colour rs = ray.intersection.mat->specular;
	Colour Ia = _col_ambient;
	Colour Id = _col_diffuse;
	Colour Is = _col_specular;

	// direction of light source (incident direction)
	Vector3D s = _pos - ray.intersection.point;
	s.normalize();
	// outward unit surface normal
	Vector3D n = ray.intersection.normal;
	n.normalize();
	// emittant direction
	Vector3D r = -s + 2*(n.dot(s)) * n;
	r.normalize();
	// unit vector in camera direction
	double magnitude = sqrt(ray.dir[0] * ray.dir[0] + ray.dir[1] * ray.dir[1] + ray.dir[2] * ray.dir[2]); 
	Vector3D b(-ray.dir[0]/magnitude, -ray.dir[1]/magnitude, -ray.dir[2]/magnitude);
	b.normalize();

	/****************** Scene signiture ************************/
	//ray.col = ray.intersection.mat->diffuse;
	/****************** Ambient + diffuse  ************************/
	//ray.col = ra * Ia + fmax(0, s.dot(n)) * rd * Id;
	/****************** Phong shading  ************************/
	//ray.col = ra * Ia + fmax(0, s.dot(n)) * rd * Id + pow(fmax(0, r.dot(b)), ray.intersection.mat->specular_exp) * rs * Is; 

	//ray.col.clamp();


	// Normalized vectors needed for phong shading
    Vector3D N = ray.intersection.normal; // normal
    Vector3D L = _pos - ray.intersection.point; // light source direction
    Vector3D V = -ray.dir; // reflection
    Vector3D R = 2.*(L.dot(N)*N)-L; // Perfect Mirror Specular Reflection

    N.normalize();
    L.normalize();
    V.normalize();
    R.normalize();

    //intensity due to ambient light
    Colour Ia = (ray.intersection.mat->ambient) * _col_ambient;

    //intensity due to diffuse light
    Colour Id = (ray.intersection.mat->diffuse) * (std::max(0.0, N.dot(L))*_col_diffuse);

    //intensity due to specular light
    Colour Is = (ray.intersection.mat->specular) * (std::max(0.0, pow(V.dot(R),(*ray.intersection.mat).specular_exp))*_col_specular);

    //Phong shading
    Ia.clamp();
    Id.clamp();
    Is.clamp();

    ray.col = ray.col + Ia + Id + Is;
}

