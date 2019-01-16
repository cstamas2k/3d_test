#pragma once
#include <string>
#include <iostream>
#include <SDL.h>
#include <fstream>
#include <strstream>
#include <algorithm>
#include <vector>

struct color {
	uint8_t r, g, b = 255;
};

struct vec2d {
	float u, v = 0.0f;
	float w = 1;
};

struct vec3d {
	float x, y, z = 0.0f;
	float w = 1;
};

struct triangle {
	vec3d p[3];
	vec2d t[3];
	color c;
};

struct mesh {
	std::vector<triangle> tris;

	bool loadFromOBJ(std::string path, bool hasTexture = false);
};

struct mat4x4 {
	float m[4][4] = { 0 };
};