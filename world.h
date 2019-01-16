#pragma once

#include "renderer.h"
#include "global.h"

class World {
	mesh cubeMesh;
public:

	World();
	std::vector<mesh> v_Cubes;
	mesh makeCube();

	void generateWorld(int x, int y, int z);
};