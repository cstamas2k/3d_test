#include "global.h"
#include "renderer.h"
#include "world.h"
#include <list>

int main(int argc, char* argv[]) {

	Renderer engine;
	SDL_Event sdl_event;


	bool run = true;

	engine.clearScreen();

	vec3d cam_loc = { 0, 0, 0 };
	vec3d vec_Up = { 0, 1, 0 };
	vec3d vec_lookAt = { 0, 1, 0 };
	float f_yaw = 1.0f;

	World worldHandler;
	mesh meshCube = worldHandler.makeCube();

	std::vector<mesh> v_world;

	vec3d startingLoc = { 0, 0, 0 };

		for (int y = 0; y < 10; y++) {
			mesh genCube;
			for (auto tris : meshCube.tris) {
				
				tris.p[0].x += y*10;  tris.p[1].x += y*10; tris.p[2].x += y*10;
				tris.p[0].y += y*10;  tris.p[1].y += y*10; tris.p[2].y += y*10;
				tris.p[0].z += y*10;  tris.p[1].z += y*10; tris.p[2].z += y*10;


				genCube.tris.push_back(tris);
			}
			v_world.push_back(genCube);
		}


	mat4x4 matProj;
	matProj = engine.makeProj(90.0f, 800 / 600, 0.1f, 1000.0f);
	float *depthBuff = new float[800 * 600];
	while (run) {

		//main loop start
		SDL_PollEvent(&sdl_event);
		if (sdl_event.type == SDL_KEYDOWN && sdl_event.key.keysym.sym == SDLK_ESCAPE) run = false;
		

		mat4x4 matTrans;
		matTrans = engine.makeTrans(0.0f, 0.0f, 5.0f);

		mat4x4 matWorld;
		matWorld = engine.makeIdentity();	// Form World Matrix
		matWorld = engine.mulMatrix(matWorld, matTrans); // Transform by translation

		// Create "Point At" Matrix for camera
		vec3d vUp = { 0,1,0 };
		vec3d vTarget = { 0,0,1 };
		mat4x4 matCameraRot = engine.makeRotY(f_yaw);
		vec_lookAt = engine.mat_multiplyVec(matCameraRot, vTarget);
		vTarget = engine.vec_add(cam_loc, vec_lookAt);
		mat4x4 matCamera = engine.matPointat(cam_loc, vTarget, vUp);

		// Make view matrix from camera
		mat4x4 mat_view = engine.matQuickInv(matCamera);


		std::vector<triangle> triToRaster;

		vec3d vec_fwd = engine.vec_mul(vec_lookAt, 2.0f);
		if (sdl_event.type == SDL_KEYDOWN) {
			if (sdl_event.key.keysym.sym == SDLK_w) cam_loc = engine.vec_add(cam_loc, vec_fwd);
			if (sdl_event.key.keysym.sym == SDLK_s) cam_loc = engine.vec_sub(cam_loc, vec_fwd);
			if (sdl_event.key.keysym.sym == SDLK_a) f_yaw -= 0.2f;
			if (sdl_event.key.keysym.sym == SDLK_d) f_yaw += 0.2f;
			if (sdl_event.key.keysym.sym == SDLK_UP) cam_loc.y += 2.0f;
			if (sdl_event.key.keysym.sym == SDLK_DOWN) cam_loc.y -= 2.0f;
		}

		for (auto meshes : v_world) {

			//3d math stuff
			for (auto tri : meshCube.tris) {
				triangle triProjected, triTransformed, triViewed;

				triTransformed.p[0] = engine.mat_multiplyVec(matWorld, tri.p[0]);
				triTransformed.p[1] = engine.mat_multiplyVec(matWorld, tri.p[1]);
				triTransformed.p[2] = engine.mat_multiplyVec(matWorld, tri.p[2]);
				triTransformed.t[0] = tri.t[0];
				triTransformed.t[1] = tri.t[1];
				triTransformed.t[2] = tri.t[2];

				vec3d normal, line1, line2;
				line1 = engine.vec_sub(triTransformed.p[1], triTransformed.p[0]);
				line2 = engine.vec_sub(triTransformed.p[2], triTransformed.p[0]);

				normal = engine.vec_crossProduct(line1, line2);
				normal = engine.vec_normalize(normal);

				vec3d cam_ray = engine.vec_sub(triTransformed.p[0], cam_loc);

				if (engine.vec_dotProduct(normal, cam_ray) < 0.0f) {
					vec3d light_coming_from = { 0.0f, 1.0f, -1.0f };
					light_coming_from = engine.vec_normalize(light_coming_from);

					//float dp = std::max(0.1f, engine.vec_dotProduct(light_coming_from, normal));

					//triTransformed.c.r = dp * 255;

					triViewed.p[0] = engine.mat_multiplyVec(mat_view, triTransformed.p[0]);
					triViewed.p[1] = engine.mat_multiplyVec(mat_view, triTransformed.p[1]);
					triViewed.p[2] = engine.mat_multiplyVec(mat_view, triTransformed.p[2]);
					triViewed.t[0] = triTransformed.t[0];
					triViewed.t[1] = triTransformed.t[1];
					triViewed.t[2] = triTransformed.t[2];
					triViewed.c = triTransformed.c;

					int nClippedTriangles = 0;
					triangle clipped[2];
					nClippedTriangles = engine.tri_clipAgPlane({ 0.0f, 0.0f, 0.1f }, { 0.0f, 0.0f, 1.0f }, triViewed, clipped[0], clipped[1]);

					for (int n = 0; n < nClippedTriangles; n++)
					{
						// Project triangles from 3D --> 2D
						triProjected.p[0] = engine.mat_multiplyVec(matProj, clipped[n].p[0]);
						triProjected.p[1] = engine.mat_multiplyVec(matProj, clipped[n].p[1]);
						triProjected.p[2] = engine.mat_multiplyVec(matProj, clipped[n].p[2]);
						triProjected.c = clipped[n].c;
						triProjected.t[0] = clipped[n].t[0];
						triProjected.t[1] = clipped[n].t[1];
						triProjected.t[2] = clipped[n].t[2];


						triProjected.t[0].u = triProjected.t[0].u / triProjected.p[0].w;
						triProjected.t[1].u = triProjected.t[1].u / triProjected.p[1].w;
						triProjected.t[2].u = triProjected.t[2].u / triProjected.p[2].w;

						triProjected.t[0].v = triProjected.t[0].v / triProjected.p[0].w;
						triProjected.t[1].v = triProjected.t[1].v / triProjected.p[1].w;
						triProjected.t[2].v = triProjected.t[2].v / triProjected.p[2].w;

						triProjected.t[0].w = 1.0f / triProjected.p[0].w;
						triProjected.t[1].w = 1.0f / triProjected.p[1].w;
						triProjected.t[2].w = 1.0f / triProjected.p[2].w;


						// Scale into view, we moved the normalising into cartesian space
						// out of the matrix.vector function from the previous videos, so
						// do this manually
						triProjected.p[0] = engine.vec_div(triProjected.p[0], triProjected.p[0].w);
						triProjected.p[1] = engine.vec_div(triProjected.p[1], triProjected.p[1].w);
						triProjected.p[2] = engine.vec_div(triProjected.p[2], triProjected.p[2].w);

						// X/Y are inverted so put them back
						triProjected.p[0].x *= -1.0f;
						triProjected.p[1].x *= -1.0f;
						triProjected.p[2].x *= -1.0f;
						triProjected.p[0].y *= -1.0f;
						triProjected.p[1].y *= -1.0f;
						triProjected.p[2].y *= -1.0f;

						// Offset verts into visible normalised space
						vec3d vOffsetView = { 1,1,0 };
						triProjected.p[0] = engine.vec_add(triProjected.p[0], vOffsetView);
						triProjected.p[1] = engine.vec_add(triProjected.p[1], vOffsetView);
						triProjected.p[2] = engine.vec_add(triProjected.p[2], vOffsetView);
						triProjected.p[0].x *= 0.5f * (float)800;
						triProjected.p[0].y *= 0.5f * (float)600;
						triProjected.p[1].x *= 0.5f * (float)800;
						triProjected.p[1].y *= 0.5f * (float)600;
						triProjected.p[2].x *= 0.5f * (float)800;
						triProjected.p[2].y *= 0.5f * (float)600;

						// Store triangle for sorting
						triToRaster.push_back(triProjected);
					}

				}
			}

		}

		//end of 3d math stuff
		//actual rendering

		engine.clearScreen();
		for (int i = 0; i < 800 * 600; i++)
			depthBuff[i] = 0.0f;

		for (auto &tri_Rast : triToRaster)
		{
			// Clip triangles against all four screen edges, this could yield
			// a bunch of triangles, so create a queue that we traverse to 
			//  ensure we only test new triangles generated against planes
			triangle clipped[2];
			std::list<triangle> listTriangles;

			// Add initial triangle
			listTriangles.push_back(tri_Rast);
			int nNewTriangles = 1;

			for (int p = 0; p < 4; p++)
			{
				int nTrisToAdd = 0;
				while (nNewTriangles > 0)
				{
					// Take triangle from front of queue
					triangle test = listTriangles.front();
					listTriangles.pop_front();
					nNewTriangles--;

					// Clip it against a plane. We only need to test each 
					// subsequent plane, against subsequent new triangles
					// as all triangles after a plane clip are guaranteed
					// to lie on the inside of the plane. I like how this
					// comment is almost completely and utterly justified
					switch (p)
					{
					case 0:	nTrisToAdd = engine.tri_clipAgPlane({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					case 1:	nTrisToAdd = engine.tri_clipAgPlane({ 0.0f, (float)600 - 1, 0.0f }, { 0.0f, -1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					case 2:	nTrisToAdd = engine.tri_clipAgPlane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					case 3:	nTrisToAdd = engine.tri_clipAgPlane({ (float)800 - 1, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					}

					// Clipping may yield a variable number of triangles, so
					// add these new ones to the back of the queue for subsequent
					// clipping against next planes
					for (int w = 0; w < nTrisToAdd; w++)
						listTriangles.push_back(clipped[w]);
				}
				nNewTriangles = listTriangles.size();
			}


			// Draw the transformed, viewed, clipped, projected, sorted, clipped triangles
			for (auto &t : listTriangles)
			{
				engine.drawTri(t);
			}
		}


		//main loop end
		SDL_Delay(16);
		SDL_RenderPresent(engine.r);
		//std::cout << cam_loc.x << "|" << cam_loc.y << "|" << cam_loc.z << std::endl;
	}

	return 0;
}