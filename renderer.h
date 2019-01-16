#pragma once

#include "global.h"

class Renderer
{
public:
	struct {
		int x = 800;
		int y = 600;
	} window_size;

	SDL_Renderer *r;
	SDL_Window *w;

	Renderer();

	~Renderer();

	void clearScreen();

	struct Line {
		int x1, x2, y1, y2 = 0;
	};
	void drawLine(Line line, uint8_t red, uint8_t green, uint8_t blue);

	/* 3d stuff begins here*/



	vec3d mat_multiplyVec(mat4x4 &m, vec3d i);
	vec3d vec_add(vec3d &a, vec3d &b);
	vec3d vec_sub(vec3d &a, vec3d &b);
	vec3d vec_mul(vec3d &a, float b);
	vec3d vec_div(vec3d &a, float b);
	float vec_dotProduct(vec3d &a, vec3d &b);
	float vec_lenght(vec3d &a);
	vec3d vec_normalize(vec3d &a);
	vec3d vec_crossProduct(vec3d &a, vec3d &b);
	vec3d vec_intersectPlane(vec3d &plane_p, vec3d &plane_n, vec3d &lineStart, vec3d &lineEnd, float &t);
	int	  tri_clipAgPlane(vec3d plane_p, vec3d plane_n, triangle &in_tri, triangle &out_tri1, triangle &out_tri2);

	mat4x4 makeIdentity();

	mat4x4 makeRotX(float angle);
	mat4x4 makeRotY(float angle);
	mat4x4 makeRotZ(float angle);

	mat4x4 makeTrans(float x, float y, float z);
	mat4x4 makeProj(float fovDeg, float aRatio, float near, float far);
	mat4x4 mulMatrix(mat4x4 &m1, mat4x4 &m2);
	mat4x4 matPointat(vec3d &pos, vec3d &target, vec3d &up);
	mat4x4 matQuickInv(mat4x4 &m);

	void drawTri(triangle tri);
};

