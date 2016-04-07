/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		Implementations of functions in raytracer.h, 
		and the main function which specifies the 
		scene to be rendered.	

***********************************************************/


#include "raytracer.h"
#include "bmp_io.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

Raytracer::Raytracer() : _lightSource(NULL) {
	_root = new SceneDagNode();
}

Raytracer::~Raytracer() {
	delete _root;
}

SceneDagNode* Raytracer::addObject( SceneDagNode* parent, 
		SceneObject* obj, Material* mat ) {
	SceneDagNode* node = new SceneDagNode( obj, mat );
	node->parent = parent;
	node->next = NULL;
	node->child = NULL;
	
	// Add the object to the parent's child list, this means
	// whatever transformation applied to the parent will also
	// be applied to the child.
	if (parent->child == NULL) {
		parent->child = node;
	}
	else {
		parent = parent->child;
		while (parent->next != NULL) {
			parent = parent->next;
		}
		parent->next = node;
	}
	
	return node;;
}

LightListNode* Raytracer::addLightSource( LightSource* light ) {
	LightListNode* tmp = _lightSource;
	_lightSource = new LightListNode( light, tmp );
	return _lightSource;
}

void Raytracer::addAreaLightSource(int orient, double l, double w, double density, Point3D origin, Colour col)
{
	/*
	col[0] = col[0] / ((density+1)*(density+1));
	col[1] = col[1] / ((density+1)*(density+1));
	col[2] = col[2] / ((density+1)*(density+1));
	*/
	if( orient == 1 )
	{// Area light in the yz plane
		for (double i = origin[1]-l/2; i <= origin[1]+l/2; i+=l/density)
		{
			for (double j= origin[2]-w/2; j<=origin[2]+w/2; j+=w/density)
			{
				addLightSource( new PointLight(Point3D(origin[0],i,j), col));
			}
		}
	}
	if( orient == 2 )
	{// Area light in the xz plane
		for (double i = origin[0]-l/2; i <= origin[0]+l/2; i+=l/density)
		{
			for (double j= origin[2]-w/2; j<=origin[2]+w/2; j+=w/density)
			{
				addLightSource( new PointLight(Point3D(i,origin[1],j), col));
			}
		}
	}
	if( orient == 3 )
	{// Area light in the xy plane
		for (double i = origin[0]-l/2; i <= origin[0]+l/2; i+=l/density)
		{
			for (double j= origin[1]-w/2; j<=origin[1]+w/2; j+=w/density)
			{
				addLightSource( new PointLight(Point3D(i,j,origin[2]), col));
			}
		}
	}
}

void Raytracer::rotate( SceneDagNode* node, char axis, double angle ) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;
	
	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
			break;
		}
		if (i == 0) {
		    node->trans = node->trans*rotation; 	
			angle = -angle;
		} 
		else {
			node->invtrans = rotation*node->invtrans; 
		}	
	}
}

void Raytracer::translate( SceneDagNode* node, Vector3D trans ) {
	Matrix4x4 translation;
	
	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	node->trans = node->trans*translation; 	
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	node->invtrans = translation*node->invtrans; 
}

void Raytracer::scale( SceneDagNode* node, Point3D origin, double factor[3] ) {
	Matrix4x4 scale;
	
	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	node->trans = node->trans*scale; 	
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	node->invtrans = scale*node->invtrans; 
}

Matrix4x4 Raytracer::initInvViewMatrix( Point3D eye, Vector3D view, 
		Vector3D up ) {
	Matrix4x4 mat; 
	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);

	mat[0][0] = w[0];
	mat[1][0] = w[1];
	mat[2][0] = w[2];
	mat[0][1] = up[0];
	mat[1][1] = up[1];
	mat[2][1] = up[2];
	mat[0][2] = -view[0];
	mat[1][2] = -view[1];
	mat[2][2] = -view[2];
	mat[0][3] = eye[0];
	mat[1][3] = eye[1];
	mat[2][3] = eye[2];

	return mat; 
}

void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray ) {
	SceneDagNode *childPtr;

	// Applies transformation of the current node to the global
	// transformation matrices.
	_modelToWorld = _modelToWorld*node->trans;
	_worldToModel = node->invtrans*_worldToModel; 
	if (node->obj) {
		// Perform intersection.
		if (node->obj->intersect(ray, _worldToModel, _modelToWorld)) {
			ray.intersection.mat = node->mat;
		}
	}
	// Traverse the children.
	childPtr = node->child;
	while (childPtr != NULL) {
		traverseScene(childPtr, ray);
		childPtr = childPtr->next;
	}

	// Removes transformation of the current node from the global
	// transformation matrices.
	_worldToModel = node->trans*_worldToModel;
	_modelToWorld = _modelToWorld*node->invtrans;
}

void Raytracer::computeShading( Ray3D& ray ) {
	LightListNode* curLight = _lightSource;
	for (;;) {
		if (curLight == NULL) break;

		/****************** shadow ************************/
		/*
		bool inShadow = false;
		Ray3D shadowRay;
		Point3D lightPos = (curLight->light)->get_position();
		shadowRay.origin = lightPos;
		shadowRay.dir = ray.intersection.point - lightPos;
		traverseScene(_root, shadowRay); 
		if (shadowRay.intersection.t_value > 0.0001 && shadowRay.intersection.t_value < 0.999) {
			inShadow = true;
		}

		if (!inShadow) {
			curLight->light->shade(ray);
			curLight = curLight->next;
		}
		else {
			curLight = curLight->next;
		}
		*/

		/****************** No shadow ************************/
		
		curLight->light->shade(ray);
		curLight = curLight->next;
		
	}
}

void Raytracer::initPixelBuffer() {
	int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
	_rbuffer = new unsigned char[numbytes];
	_gbuffer = new unsigned char[numbytes];
	_bbuffer = new unsigned char[numbytes];
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			_rbuffer[i*_scrWidth+j] = 0;
			_gbuffer[i*_scrWidth+j] = 0;
			_bbuffer[i*_scrWidth+j] = 0;
		}
	}
}

void Raytracer::flushPixelBuffer( char *file_name ) {
	bmp_write( file_name, _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer );
	delete _rbuffer;
	delete _gbuffer;
	delete _bbuffer;
}

// Calculate the direction of the reflected ray
Ray3D reflectRay(Ray3D& ray){
	ray.dir.normalize();
	ray.intersection.normal.normalize();
	Vector3D direction = ray.dir - 2 * (ray.dir.dot(ray.intersection.normal)) * ray.intersection.normal;
	return Ray3D (ray.intersection.point, direction);
}

Colour Raytracer::shadeRay( Ray3D& ray, int depth ) {
	Colour col(0.0, 0.0, 0.0); 
	traverseScene(_root, ray); 
	
	// Don't bother shading if the ray didn't hit 
	// anything.
	if (!ray.intersection.none) {
		computeShading(ray);

		// You'll want to call shadeRay recursively (with a different ray, 
	    // of course) here to implement reflection/refraction effects.  
		if (depth < _maxDepth) {

			/****************** No glossy reflection or refraction ************************/
			
			col = ray.col;
			col.clamp();
			

			/****************** Glossy reflection ************************/
			/*
	        double reflectionCoef = 0;

	        // Calculate reflection ray
	        Vector3D reflection_N = ray.intersection.normal;
	        Vector3D reflection_D = ray.dir;
	        Vector3D reflectionDir = reflection_D - 2 * reflection_D.dot(reflection_N) * reflection_N;
	        reflectionDir.normalize();
	        Point3D reflectionOrigin = ray.intersection.point + 0.01*reflectionDir;
	        Ray3D reflectionRay = Ray3D(reflectionOrigin, reflectionDir);

	        shadeRay(reflectionRay, depth+1);

	        if (reflectionRay.intersection.t_value > 10 || reflectionRay.intersection.t_value <= 0) {
	            col = ray.col;
	        }
	        else {
		        reflectionCoef = fabs(1.0/reflectionRay.intersection.t_value);
		        reflectionCoef = fmax(0, fmin(reflectionCoef, 0.8));
		        col = ray.col + reflectionCoef * reflectionRay.col;
	        }
			col.clamp();
			*/

	        /****************** Refraction ************************/
			/*
	        double refractionCoef = 0;
	        double index = 0.8;

	        Vector3D refraction_N = ray.intersection.normal;
	        Vector3D refraction_D = ray.dir;	
	        double k = 1 - index * index * (1 - refraction_D.dot(refraction_N) * refraction_D.dot(refraction_N));
	        Vector3D refractionDir(0, 0, 0);
	        if (k >= 0) {
	        	refractionDir = index * refraction_D - (index * refraction_D.dot(refraction_N) + sqrt(k)) * refraction_N;
	        }
			refractionDir.normalize();
			Point3D refractionOrigin = ray.intersection.point + 0.01*refractionDir;
			Ray3D refractionRay = Ray3D(refractionOrigin, refractionDir);

			shadeRay(refractionRay, depth+1);

	        if (refractionRay.intersection.t_value > 10 || refractionRay.intersection.t_value <= 0) {
	            col = ray.col;
	        }
	        else {
				refractionCoef = fabs(1/refractionRay.intersection.t_value);
				refractionCoef = fmax(0, fmin(refractionCoef, 0.8));
				col = ray.col + refractionCoef * refractionRay.col;
			}
	        col.clamp();
	        */
	    }
	}
	return col; 
}	

void Raytracer::render( int width, int height, Point3D eye, Vector3D view, 
		Vector3D up, double fov, char* fileName ) {
	Matrix4x4 viewToWorld;
	_scrWidth = width;
	_scrHeight = height;
	double factor = (double(height)/2)/tan(fov*M_PI/360.0);

	_maxDepth = 20;
	
	initPixelBuffer();
	viewToWorld = initInvViewMatrix(eye, view, up);

	
	/****************** Construct 4 rays for each pixel. (Anti-alising) ************************/ 
	/*
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			for(float k = i; k < i + 1.0f; k += 0.5f){
				for(float l = j; l < j + 1.0f; l += 0.5f){
					// Sets up ray origin and direction in view space,
					// image plane is at z = -1.
					Point3D origin(0, 0, 0);
					Point3D imagePlane;
					imagePlane[0] = (-double(width)/2 + 0.5 + l)/factor;
					imagePlane[1] = (-double(height)/2 + 0.5 + k)/factor;
					imagePlane[2] = -1;

					// TODO: Convert ray to world space and call
					// shadeRay(ray) to generate pixel colour.
					float incre = 0.25f;
					Ray3D ray;
					ray.origin = viewToWorld*origin;
					ray.dir = viewToWorld*(imagePlane-origin);

					Colour col = shadeRay(ray);
					_rbuffer[i*width+j] += int(col[0]*255*incre);
					_gbuffer[i*width+j] += int(col[1]*255*incre);
					_bbuffer[i*width+j] += int(col[2]*255*incre);
				}
			}
		}
	}
	*/


	/****************** Construct 1 ray for each pixel. (No anti-alising) ************************/ 
	/*
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			// Sets up ray origin and direction in view space,
			// image plane is at z = -1.
			Point3D origin(0, 0, 0);
			Point3D imagePlane;
			imagePlane[0] = (-double(width)/2 + 0.5 + j)/factor;
			imagePlane[1] = (-double(height)/2 + 0.5 + i)/factor;
			imagePlane[2] = -1;

			// TODO: Convert ray to world space and call
			// shadeRay(ray) to generate pixel colour.
			Ray3D ray;
			viewToWorld = initInvViewMatrix(eye, view, up);
			ray.origin = viewToWorld*origin;
			ray.dir = viewToWorld*(imagePlane-origin);

			Colour col = shadeRay(ray, 0);
			_rbuffer[i*width+j] += int(col[0]*255);
			_gbuffer[i*width+j] += int(col[1]*255);
			_bbuffer[i*width+j] += int(col[2]*255);
		}
	}
	*/

	/****************** Depth of field ************************/ 
	for (int e=0; e<8; e++) {
		// Shake eyes to create blur effect
		eye[0] = eye[0] + 0.02*cos(0 + e * (2*M_PI/8.0));
		eye[1] = eye[1] + 0.02*sin(0 + e * (2*M_PI/8.0));

		for (int i = 0; i < _scrHeight; i++) {
			for (int j = 0; j < _scrWidth; j++) {
				// Sets up ray origin and direction in view space,
				// image plane is at z = -1.
				Point3D origin(0, 0, 0);
				Point3D imagePlane;
				imagePlane[0] = (-double(width)/2 + 0.5 + j)/factor;
				imagePlane[1] = (-double(height)/2 + 0.5 + i)/factor;
				imagePlane[2] = -1;

				// TODO: Convert ray to world space and call
				// shadeRay(ray) to generate pixel colour.
				Ray3D ray;
				
				ray.origin = viewToWorld*imagePlane;
				// ray.dir = viewToWorld*(imagePlane-origin);
				ray.dir = ray.origin - eye;
				ray.dir.normalize();

				Colour col = shadeRay(ray, 0);
				_rbuffer[i*width+j] += int(col[0]*255*0.125);
				_gbuffer[i*width+j] += int(col[1]*255*0.125);
				_bbuffer[i*width+j] += int(col[2]*255*0.125);
			}
		}
	}
	flushPixelBuffer(fileName);
}

int main(int argc, char* argv[])
{	
	
	// Build your scene and setup your camera here, by calling 
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it from two different view points, DO NOT
	// change this if you're just implementing part one of the 
	// assignment.  
	Raytracer raytracer;
	int width = 320; 
	int height = 240; 

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}

	// Camera parameters.
	Point3D eye(0, 0, 1);
	Vector3D view(0, 0, -1);
	Vector3D up(0, 1, 0);
	double fov = 60;

	// Defines a material for shading.
	Material gold( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648), 
			Colour(0.628281, 0.555802, 0.366065), 
			51.2 );
	Material jade( Colour(0, 0, 0), Colour(0.54, 0.89, 0.63), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8 );

	/****************** An area light source. ************************/ 
	//raytracer.addAreaLightSource(3, 2,2,8, Point3D(0,0,5), Colour(0.01, 0.01, 0.01) );
	
	/****************** A point light source. ************************/ 
	raytracer.addLightSource( new PointLight(Point3D(0, 0, 5), Colour(0.9, 0.9, 0.9) ) );

	/*
	// Add a unit square into the scene with material mat.
	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &gold );
	SceneDagNode* plane = raytracer.addObject( new UnitSquare(), &jade );
	SceneDagNode* cylinder = raytracer.addObject(new UnitCylinder(), &gold);

	// Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 2.0, 1.0 };
	double factor2[3] = { 6.0, 6.0, 6.0 };
	double factor3[3] = { 1.5, 1.5, 3.0 };
	*/

	SceneDagNode* sphere1 = raytracer.addObject( new UnitSphere(), &gold );
	SceneDagNode* sphere2 = raytracer.addObject( new UnitSphere(), &gold );
	SceneDagNode* sphere3 = raytracer.addObject( new UnitSphere(), &gold );
	SceneDagNode* sphere4 = raytracer.addObject( new UnitSphere(), &jade );
	double factor4[3] = { 1.0, 1.0, 1.0 };
	double factor5[3] = { 1.5, 1.5, 1.5 };
	double factor6[3] = { 2.0, 2.0, 2.0 };
	double factor7[3] = { 6.0, 6.0, 6.0 };
	
	raytracer.translate(sphere4, Vector3D(0, 0, -10));	
	raytracer.rotate(sphere4, 'x', -45); 
	raytracer.rotate(sphere4, 'z', 45); 
	raytracer.scale(sphere4, Point3D(0, 0, 0), factor7);

	raytracer.translate(sphere3, Vector3D(0, 2, -6));	
	raytracer.rotate(sphere3, 'x', -45); 
	raytracer.rotate(sphere3, 'z', 45); 
	raytracer.scale(sphere3, Point3D(0, 0, 0), factor6);

	raytracer.translate(sphere2, Vector3D(0, 0, -3));	
	raytracer.rotate(sphere2, 'x', -45); 
	raytracer.rotate(sphere2, 'z', 45); 
	raytracer.scale(sphere2, Point3D(0, 0, 0), factor5);

	raytracer.translate(sphere1, Vector3D(0, -1, 0));	
	raytracer.rotate(sphere1, 'x', -45); 
	raytracer.rotate(sphere1, 'z', 45); 
	raytracer.scale(sphere1, Point3D(0, 0, 0), factor4);

	/*
	raytracer.translate(sphere, Vector3D(0, 0, -5));	
	raytracer.rotate(sphere, 'x', -45); 
	raytracer.rotate(sphere, 'z', 45); 
	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);

	raytracer.translate(plane, Vector3D(0, 0, -7));	
	raytracer.rotate(plane, 'z', 45); 
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);

	raytracer.translate(cylinder, Vector3D(-4, 0, -5));
	raytracer.scale(cylinder, Point3D(0,0,0), factor3);
	*/

	// Render the scene, feel free to make the image smaller for
	// testing purposes.	
	raytracer.render(width, height, eye, view, up, fov, "view1.bmp");
	
	//raytracer.translate(cylinder, Vector3D(-1, 0, 0));
	// Render it from a different point of view.
	Point3D eye2(4, 2, 1);
	Vector3D view2(-4, -2, -6);
	raytracer.render(width, height, eye2, view2, up, fov, "view2.bmp");
	
	return 0;
}

