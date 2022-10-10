#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <CanvasPoint.h>
#include <TextureMap.h>
#include <Colour.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>


#define WIDTH 320
#define HEIGHT 240

// ---------------------- Week 1 ---------------------- //
void draw(DrawingWindow &window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = rand() % 256;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

// ---------------------- Week 2 ---------------------- //
std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValue) {
    std::vector<float> result;
    float difference = (to - from) / float(numberOfValue - 1);
    for(int i=0; i<numberOfValue; i++) result.push_back(from + (float)i*difference);
    return result;
}

void greyscaleInterpolation(DrawingWindow &window) {
    window.clearPixels();
    std::vector<float> interpolation = interpolateSingleFloats(255, 0, int(window.width));
    for (size_t x = 0; x < window.width; x++) {
        float red = interpolation[x];
        float green = interpolation[x];
        float blue = interpolation[x];
        for (size_t y = 0; y < window.height; y++) {
            uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
            window.setPixelColour(x, y, colour);
        }
    }
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValue) {
    std::vector<glm::vec3> result;
    glm::vec3 difference = (to - from) * (1 / (float(numberOfValue) - 1));
    for(int i=0; i<numberOfValue; i++) {
        result.push_back(from + difference*float(i));
    }
    return result;
}

void twoDimensionalColourInterpolation(DrawingWindow &window) {
    window.clearPixels();
    glm::vec3 topLeft(255, 0, 0);        // red
    glm::vec3 topRight(0, 0, 255);       // blue
    glm::vec3 bottomRight(0, 255, 0);    // green
    glm::vec3 bottomLeft(255, 255, 0);   // yellow
    std::vector<glm::vec3> left = interpolateThreeElementValues(topLeft, bottomLeft, int(window.height));
    std::vector<glm::vec3> right = interpolateThreeElementValues(topRight, bottomRight, int(window.width));
    for(size_t y=0; y<window.height; y++) {
        std::vector<glm::vec3> row = interpolateThreeElementValues(left[y], right[y], int(window.width));
        for(size_t x=0; x<window.width; x++) {
            float red = row[x].x;
            float green = row[x].y;
            float blue = row[x].z;
            uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
            window.setPixelColour(x, y, colour);
        }
    }
}

// ---------------------- Week 3 ---------------------- //
void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, const Colour& colour) {
    // Calculate step size
    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float numberOfSteps = fmax(abs(xDiff), abs(yDiff));
    float xStepSize = xDiff/numberOfSteps;
    float yStepSize = yDiff/numberOfSteps;
    // Calculate colour
    int red = colour.red;
    int green = colour.green;
    int blue = colour.blue;
    uint32_t c = (255 << 24) + (red << 16) + (green << 8) + (blue);
    // Set colour for pixel
    for (int i = 0.0; float(i) < numberOfSteps; i++) {
        float x = from.x + (xStepSize * float(i));
        float y = from.y + (yStepSize * float(i));
        window.setPixelColour(int(x), int(y), c);
    }
}

void drawStrokedTriangle(DrawingWindow &window, CanvasTriangle triangle, const Colour& colour) {
    drawLine(window, triangle.v0(), triangle.v1(), colour);
    drawLine(window, triangle.v1(), triangle.v2(), colour);
    drawLine(window, triangle.v2(), triangle.v0(), colour);
}

// judge triangle is top-to-bottom or bottom-to-top
int sign(float d){
    if (d<0) return -1; // bottom-to-top
    else if (d == 0) return 0;
    else return 1; // top-to-bottom
}

void triangleRasteriser(DrawingWindow &window, CanvasPoint v1, CanvasPoint v2, CanvasPoint v3, const Colour& colour) {
    float dy =abs(v2.y-v1.y);
    float x = v1.x;
    float y = v1.y;
    float x2 = v2.x;
    float x3 = v3.x;
    float y3 = v3.y;
    int s = sign(y - y3);
    CanvasPoint newV1=v1;
    CanvasPoint newV2=v1;
    for (int i = 1; float(i)<=dy; i++)  {
        float f = float(i) / dy;
        newV1.x = x + (x2 - x) * f; // Calculate the start point
        newV1.y -= float(s); // move y with 1 step each time
        newV2.x = x + (x3 - x) * f; // Calculate the end point
        newV2.y -= float(s); // move y with 1 step each time
        drawLine(window, newV1, newV2, colour);
    }

}

void drawFilledTriangle(DrawingWindow &window, CanvasTriangle triangle, const Colour& colour) {
    drawStrokedTriangle(window, triangle, Colour(255, 255, 255));
    // sort point
    CanvasPoint topPoint;
    CanvasPoint bottomPoint;
    CanvasPoint givenPoint;
    CanvasPoint extraPoint;
    if(triangle.v0().y <= triangle.v1().y && triangle.v0().y <= triangle.v2().y) {
        topPoint = triangle.v0();
        if(triangle.v1().y >= triangle.v2().y) {
            bottomPoint = triangle.v1();
            givenPoint = triangle.v2();
        } else {
            bottomPoint = triangle.v2();
            givenPoint = triangle.v1();
        }
    } else if(triangle.v1().y <= triangle.v0().y && triangle.v1().y <= triangle.v2().y) {
        topPoint = triangle.v1();
        if(triangle.v0().y >= triangle.v2().y) {
            bottomPoint = triangle.v0();
            givenPoint = triangle.v2();
        } else {
            bottomPoint = triangle.v2();
            givenPoint = triangle.v0();
        }
    } else if(triangle.v2().y <= triangle.v0().y && triangle.v2().y <= triangle.v1().y) {
        topPoint = triangle.v2();
        if(triangle.v1().y >= triangle.v0().y) {
            bottomPoint = triangle.v1();
            givenPoint = triangle.v0();
        } else {
            bottomPoint = triangle.v0();
            givenPoint = triangle.v1();
        }
    }
    // separate triangle
    double similarTriangleRadio = (bottomPoint.x - topPoint.x)/(bottomPoint.y - topPoint.y);
    float extraPointX;
    // when radio is negative, the extra point is behind on the topPoint (-); when it is positive, the extraPoint is in the front of topPoint (+)
    if(similarTriangleRadio >= 0) {
        extraPointX = float(topPoint.x + ((givenPoint.y - topPoint.y) * abs(similarTriangleRadio)));
    } else {
        extraPointX = float(topPoint.x - ((givenPoint.y - topPoint.y) * abs(similarTriangleRadio)));
    }
    extraPoint = CanvasPoint(extraPointX, givenPoint.y);

    // test separate line
    // std::cout << topPoint << std::endl;
    // std::cout << bottomPoint << extraPointX << std::endl;
    // std::cout << givenPoint << extraPointX << std::endl;
    // std::cout << "Extra Point X " << extraPointX << std::endl;
    // window.setPixelColour(extraPoint.x, extraPoint.y, (255 << 24) + (255 << 16) + (0 << 8) + (0));
    // drawLine(window, CanvasPoint(extraPoint.x, extraPoint.y), CanvasPoint(givenPoint.x, givenPoint.y), colour);

    // filled the first triangle
    triangleRasteriser(window, topPoint, givenPoint, extraPoint, colour);
    // filled the second triangle
    triangleRasteriser(window, bottomPoint, givenPoint, extraPoint, colour);
}

// Read texture map from image file
TextureMap getTextureMap(const std::string& image) {
    TextureMap textureMap = TextureMap(image);
    std::cout << "width " << textureMap.width << " height " << textureMap.height << std::endl;
    return textureMap;
}

// Match the point in texture map
uint32_t textureMapper(CanvasTriangle triangle, CanvasPoint point, TextureMap textureMap){
    // Calculate barycentric coordinates
    float alpha = (-(point.x-triangle.v1().x)*(triangle.v2().y-triangle.v1().y)+(point.y-triangle.v1().y)*(triangle.v2().x-triangle.v1().x))/(-(triangle.v0().x-triangle.v1().x)*(triangle.v2().y-triangle.v1().y)+(triangle.v0().y-triangle.v1().y)*(triangle.v2().x-triangle.v1().x));
    float beta = (-(point.x-triangle.v2().x)*(triangle.v0().y-triangle.v2().y)+(point.y-triangle.v2().y)*(triangle.v0().x-triangle.v2().x))/(-(triangle.v1().x-triangle.v2().x)*(triangle.v0().y-triangle.v2().y)+(triangle.v1().y-triangle.v2().y)*(triangle.v0().x-triangle.v2().x));
    float gamma = 1 - alpha - beta;
    // Get texture point
    CanvasPoint texturePoint((alpha * triangle.v0().texturePoint.x + beta * triangle.v1().texturePoint.x + gamma * triangle.v2().texturePoint.x), (alpha * triangle.v0().texturePoint.y + beta * triangle.v1().texturePoint.y + gamma * triangle.v2().texturePoint.y));
    // Locate the texture point position in texture map pixel list
    int index = int(texturePoint.y) * int(textureMap.width) + int(texturePoint.x);
    std::cout << "Index of texture is:" << index << " corresponds to texture point x:" << texturePoint.x << " and y: " << texturePoint.y << "\n" << std::endl;
    uint32_t colour = textureMap.pixels[index - 1];
    return colour;
}

void textureTriangle(DrawingWindow &window, const TextureMap& textureMap, CanvasTriangle triangle, CanvasPoint v1, CanvasPoint v2, CanvasPoint v3) {
    float dy =abs(v2.y-v1.y);
    float x = v1.x;
    float y = v1.y;
    float x2 = v2.x;
    float x3 = v3.x;
    float y3 = v3.y;
    int s = sign(y - y3);
    CanvasPoint newV1=v1;
    CanvasPoint newV2=v1;
    for (int i = 1; float(i)<=dy; i++)  {
        float f = float(i) / dy;
        newV1.x = x + (x2 - x) * f; // Calculate the start point
        newV1.y -= float(s); // move y with 1 step each time
        newV2.x = x + (x3 - x) * f; // Calculate the end point
        newV2.y -= float(s); // move y with 1 step each time
        // Calculate step size
        float numberOfSteps = newV2.x - newV1.x;
        // Set colour for pixel on one line
        for (int j = 0.0; float(j) < numberOfSteps; j++) {
            float pixel_x = newV1.x + float(j);
            float pixel_y = newV1.y;
            uint32_t colour = textureMapper(triangle, CanvasPoint(pixel_x, pixel_y), textureMap);
            window.setPixelColour(int(pixel_x), int(pixel_y), colour);
        }
    }
}

void drawTextureTriangle(DrawingWindow &window, const TextureMap& textureMap, CanvasTriangle triangle) {
    // Draw the triangle
    drawStrokedTriangle(window, triangle, Colour(255,255,255));
    // sort point
    CanvasPoint topPoint;
    CanvasPoint bottomPoint;
    CanvasPoint givenPoint;
    CanvasPoint extraPoint;
    if(triangle.v0().y <= triangle.v1().y && triangle.v0().y <= triangle.v2().y) {
        topPoint = triangle.v0();
        if(triangle.v1().y >= triangle.v2().y) {
            bottomPoint = triangle.v1();
            givenPoint = triangle.v2();
        } else {
            bottomPoint = triangle.v2();
            givenPoint = triangle.v1();
        }
    } else if(triangle.v1().y <= triangle.v0().y && triangle.v1().y <= triangle.v2().y) {
        topPoint = triangle.v1();
        if(triangle.v0().y >= triangle.v2().y) {
            bottomPoint = triangle.v0();
            givenPoint = triangle.v2();
        } else {
            bottomPoint = triangle.v2();
            givenPoint = triangle.v0();
        }
    } else if(triangle.v2().y <= triangle.v0().y && triangle.v2().y <= triangle.v1().y) {
        topPoint = triangle.v2();
        if(triangle.v1().y >= triangle.v0().y) {
            bottomPoint = triangle.v1();
            givenPoint = triangle.v0();
        } else {
            bottomPoint = triangle.v0();
            givenPoint = triangle.v1();
        }
    }
    // separate triangle
    double similarTriangleRadio = (bottomPoint.x - topPoint.x)/(bottomPoint.y - topPoint.y);
    float extraPointX;
    // when radio is negative, the extra point is behind on the topPoint (-); when it is positive, the extraPoint is in the front of topPoint (+)
    if(similarTriangleRadio >= 0) {
        extraPointX = float(topPoint.x + ((givenPoint.y - topPoint.y) * abs(similarTriangleRadio)));
    } else {
        extraPointX = float(topPoint.x - ((givenPoint.y - topPoint.y) * abs(similarTriangleRadio)));
    }
    extraPoint = CanvasPoint(extraPointX, givenPoint.y);

    // filled the first triangle
    textureTriangle(window, textureMap, triangle, topPoint, givenPoint, extraPoint);
    // filled the second triangle
    textureTriangle(window, textureMap, triangle, bottomPoint, givenPoint, extraPoint);
}

// ---------------------- Show in window ---------------------- //
void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
        else if (event.key.keysym.sym == SDLK_u){
            CanvasPoint point1,point2,point3;
            CanvasTriangle triangle;
            point1 = CanvasPoint(float(rand()%window.width), float(rand()%window.height));
            point2 = CanvasPoint(float(rand()%window.width), float(rand()%window.height));
            point3 = CanvasPoint(float(rand()%window.width), float(rand()%window.height));
            triangle = CanvasTriangle(point1, point2, point3);
            drawStrokedTriangle(window, triangle, Colour(rand()%256, rand()%256, rand()%256));
        }
        else if (event.key.keysym.sym == SDLK_f){
            CanvasPoint point1,point2,point3;
            CanvasTriangle triangle;
            point1 = CanvasPoint(float(rand()%window.width), float(rand()%window.height));
            point2 = CanvasPoint(float(rand()%window.width), float(rand()%window.height));
            point3 = CanvasPoint(float(rand()%window.width), float(rand()%window.height));
            triangle = CanvasTriangle(point1, point2, point3);
            drawFilledTriangle(window, triangle, Colour(rand()%256, rand()%256, rand()%256));
        }
        else if(event.key.keysym.sym == SDLK_o){
            CanvasPoint point1,point2,point3;
            point1 = CanvasPoint(160, 10);
            point2 = CanvasPoint(300, 230);
            point3 = CanvasPoint(10, 150);
            point1.texturePoint.x = 195; point1.texturePoint.y = 5;
            point2.texturePoint.x = 395; point2.texturePoint.y = 380;
            point3.texturePoint.x = 65; point3.texturePoint.y = 330;
            drawTextureTriangle(window, getTextureMap("texture.ppm"), CanvasTriangle(point1, point2, point3));
        }
        else if (event.key.keysym.sym == SDLK_c) window.clearPixels();
    } else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
     // Single Element Numerical Interpolation
//     std::vector<float> result;
//     result = interpolateSingleFloats(2.2, 8.5, 7);
//     for(float i : result) std::cout << i << " ";
//     std::cout << std::endl;

    // Three Element Numerical Interpolation
//     std::vector<glm::vec3> result;
//     glm::vec3 from(1, 4, 9.2);
//     glm::vec3 to(4, 1, 9.8);
//     result = interpolateThreeElementValues(from, to, 4);
//     for(size_t i=0; i<result.size(); i++) std::cout << result[i].x << " " << result[i].y << " " << result[i].z << " " << std::endl;make

	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
        // week 1
//        draw(window);
        // week 2
//        greyscaleInterpolation(window);
//        twoDimensionalColourInterpolation(window);
        // week 3
//        drawLine(window, CanvasPoint(window.width/2, window.height/2), CanvasPoint(window.width-1, window.height-1), Colour(255, 255, 255));
//        drawLine(window, CanvasPoint(window.width/2, 0), CanvasPoint(window.width/2, window.height/2), Colour(255, 0, 0));
//        drawLine(window, CanvasPoint(0, window.height-1), CanvasPoint(window.width-1, 0), Colour(0, 0, 255));

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
