#pragma once
#include <iostream>

namespace Graphics {

	struct Point {
		int x, y = 0;
	};

	struct Color {
		uint8_t red, green, blue = 0;
	};

	struct Rect {
		Point p[4];
		Color c;
	};


}