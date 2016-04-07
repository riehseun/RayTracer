/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"
using namespace std;

bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.

	// transform the ray into object space  
	Ray3D r;
	r.origin = worldToModel * ray.origin;
	r.dir = worldToModel * ray.dir;

	// normal vector and a point on the plane
	Vector3D normal(0.0,0.0,1.0);
	Point3D pointOnPlane(0.5,0.5,0.0);

	double d;
	if (r.dir.dot(normal) != 0) {
		d = ((pointOnPlane - r.origin).dot(normal)) / r.dir.dot(normal);
	}
	else {
		d = 0;
	}
	double x = r.origin[0] + (r.dir[0] * d);
	double y = r.origin[1] + (r.dir[1] * d);
	double z = 0;
	Point3D intersection(x,y,z);
	
	if (x <= 0.5 && x >= -0.5 && y <= 0.5 && y >= -0.5) {
		// if line and plane are parallel
		if (d == 0) {
			return false;
		}
		else {
			if (ray.intersection.none || d < ray.intersection.t_value) {
				ray.intersection.point = modelToWorld * intersection;
				ray.intersection.normal = worldToModel.transpose()*normal;
				ray.intersection.normal.normalize();
				ray.intersection.none = false;
				ray.intersection.t_value = d;
				return true;	
			}
			else {
				return false;
			}
		}
	}  
	else {
		return false;	
	}
}

bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
	
	// transform the ray into object space  
	Ray3D r;
	r.origin = worldToModel * ray.origin;
	r.dir = worldToModel * ray.dir;

	//sphere parameters
	Point3D centerOfSphere(0,0,0); 
	double radius = 1;
	Vector3D rayToSphere(r.origin-centerOfSphere);

	/*
	// calculate unit vector for the ray direction
	double magnitude = sqrt(r.dir[0] * r.dir[0] + r.dir[1] * r.dir[1] + r.dir[2] * r.dir[2]);  
	Vector3D unitDir(r.dir[0]/magnitude, r.dir[1]/magnitude, r.dir[2]/magnitude); 

	//compute the number inside the square root
	double squareRoot = (unitDir.dot(rayToSphere)) * unitDir.dot(rayToSphere) - rayToSphere.dot(rayToSphere) + radius * radius;
	if (squareRoot < 0) {
		return false;
	}
	*/

	double A = r.dir.dot(r.dir);
    double B = 2*(r.dir.dot(r.origin-centerOfSphere)); 
    double C = (r.origin-centerOfSphere).dot(r.origin-centerOfSphere) - radius * radius;
    double dis = B*B - 4*A*C;
    if ( dis < 0 ) {
    	return false;
    }
	else {
		//double sol1 = -1 * unitDir.dot(rayToSphere) + sqrt(squareRoot);
		//double sol2 = -1 * unitDir.dot(rayToSphere) - sqrt(squareRoot);
		double sol1 = (-B + sqrt(dis))/(2 * A);
    	double sol2 = (-B - sqrt(dis))/(2 * A);
		double d = fmin(sol1,sol2);
		Point3D intersection(r.origin[0] + r.dir[0] * d, r.origin[1] + r.dir[1] * d, r.origin[2] + r.dir[2] * d);
		Vector3D normal(r.origin[0] + r.dir[0] * d, r.origin[1] + r.dir[1] * d, r.origin[2] + r.dir[2] * d);
		if (ray.intersection.none || d < ray.intersection.t_value) {
			ray.intersection.point = modelToWorld * intersection;
			ray.intersection.normal = worldToModel.transpose()*normal;
			ray.intersection.normal.normalize();
			ray.intersection.none = false;
			ray.intersection.t_value = d;
			return true;
		}
		else {
			return false;
		}
	}
}

bool UnitCylinder::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// intersection for unit cylinder with height = 1 and top and bottom base circles are at
	// (0,0,0.5), (0,0,-0.5), radius = 1
	// idea, mathematically intersect cylinder of inifinite height with two planes
	// call these the caps at z =-0.5, 0.5
	Ray3D ro;
	ro.origin = worldToModel * ray.origin;
    ro.dir = worldToModel * ray.dir;
	Point3D O(0,0,0);
	Point3D intersectionPoint;
	Vector3D normal_1;
	Vector3D normal_2;
	Vector3D normal;

	double lambdaStar;
	double lambdaStar_1;
	double lambdaStar_2;
	double lambda_1;
	double lambda_2;

	//Use quadratic formula to find intersection
	double A = ro.dir[0]*ro.dir[0] + ro.dir[1]*ro.dir[1];
	double B = (ro.origin[0]*ro.dir[0] + ro.origin[1] * ro.dir[1]);
	double C = ro.origin[0]*ro.origin[0] + ro.origin[1]*ro.origin[1] - 1;

	//Find discriminant
	double d = B*B-A*C;

	//If the discriminant is negative there is no intersection
	//Else, get the lambda for the side of the cylinder
	if (d<0) {return false;}
    else
    {
        // calculate solutions, and take the minimum (closest) non-negative
        // number
        lambda_1 = -B/A + sqrt(d) / A;
        lambda_2 = -B/A - sqrt(d) / A;
        if (lambda_1 < 0 && lambda_2 < 0) {return false;}
        else if (lambda_1 < 0) {lambdaStar_2 = lambda_2;}
        else if (lambda_2 < 0) {lambdaStar_2 = lambda_1;}
        else {lambdaStar_2 = fmin(lambda_1,lambda_2);}
    }


	//See if ray intersetions either of the caps
	//and take the minimum position value (may intersect both caps)
	if (ro.dir[2] != 0)
    {
        lambda_1 = (-0.5-ro.origin[2])/ro.dir[2];
        lambda_2 = (0.5-ro.origin[2])/ro.dir[2];
        if (lambda_1 < lambda_2)
        {
            lambdaStar_1 = lambda_1;
            Point3D normal_temp(0,0,-1);
            normal_1 = normal_temp - O;
            normal_1.normalize();
        }
        else
        {
            lambdaStar_1 = lambda_2;
            Point3D normal_temp(0,0,1);
            normal_1 = normal_temp - O;
            normal_1.normalize();
        }
    }

	intersectionPoint = ro.origin + lambdaStar_1 * ro.dir;
	if (lambdaStar_1* lambdaStar_1 < 0.0000001) {return false; }
	
	//Use the first lambda to check if it intersects with the cap, top or bottom
	if (intersectionPoint[0]*intersectionPoint[0] + intersectionPoint[1] * intersectionPoint[1] <= 1)
	{
		if (!ray.intersection.none && lambdaStar > ray.intersection.t_value)
		{
			return false;
		}
		else
		{
			ray.intersection.point = intersectionPoint;
			ray.intersection.normal = normal_1;
			ray.intersection.t_value = lambdaStar_1;
			ray.intersection.none = false;
			return true;
		}
	}
    else
    {
        //if not intersected with the caps, use the second lamdba to check
         //if intersects with the side
        intersectionPoint = ro.origin + lambdaStar_2 * ro.dir;
        if (lambdaStar_2 * lambdaStar_2 < 0.00000001)
            return false;

        normal_2[0] = intersectionPoint[0];
        normal_2[1] = intersectionPoint[1];
        normal_2[2] = 0;
        normal_2.normalize();



        if (intersectionPoint[2] < 0.5 && intersectionPoint[2] > -0.5)
        {
            if (!ray.intersection.none > ray.intersection.t_value)
                return false;

            ray.intersection.point = modelToWorld * intersectionPoint;
            Point3D normalll;
            normalll[0] = intersectionPoint[0];
            normalll[1] = intersectionPoint[1];
            normalll[2] = 0;
            ray.intersection.normal = modelToWorld * (normalll - O);
            ray.intersection.t_value = lambdaStar_2;
            ray.intersection.none = false;
            return true;
        }
        else
            return false;
    }
}
