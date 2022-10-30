#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <CanvasPoint.h>
#include <ModelTriangle.h>
#include <TextureMap.h>
#include <Colour.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <map>
#include <glm/glm.hpp>
#include <RayTriangleIntersection.h>


#define WIDTH 320
#define HEIGHT 240
float depthBuffer[HEIGHT][WIDTH];

void clearDepthBuffer(){
    for(int y = 0; y < HEIGHT; y++){
        for(int x = 0; x < WIDTH; x++){
            depthBuffer[y][x] = INT32_MIN;
        }
    }

}

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
    float depthDiff = to.depth - from.depth;
    float numberOfSteps = fmax(abs(xDiff), abs(yDiff));
    float xStepSize = xDiff/numberOfSteps;
    float yStepSize = yDiff/numberOfSteps;
    float depthStepSize = depthDiff/numberOfSteps;

    // Calculate colour
    int red = colour.red;
    int green = colour.green;
    int blue = colour.blue;
    uint32_t c = (255 << 24) + (red << 16) + (green << 8) + (blue);
    // Set colour for pixel
    for (int i = 0.0; float(i) <= numberOfSteps; i++) {
        float x = from.x + (xStepSize * float(i));
        float y = from.y + (yStepSize * float(i));
        float z = from.depth + (depthStepSize * float(i));
        if (int(round(y)) < HEIGHT && int(round(y)) >= 0 && int(round(x)) < WIDTH && int(round(x)) >= 0) {
            if (depthBuffer[int(round(y))][int(round(x))] <= z) {
                depthBuffer[int(round(y))][int(round(x))] = z;
                window.setPixelColour(int(round(x)), int(round(y)), c);
            }
        }
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
    float dy = abs(v2.y-v1.y);
    float x = v1.x;
    float y = v1.y;
    float z1Diff = v2.depth - v1.depth;
    float x2 = v2.x;
    float z2Diff = v3.depth - v1.depth;
    float x3 = v3.x;
    float y3 = v3.y;
    int s = sign(y - y3);
    CanvasPoint newV1=v1;
    CanvasPoint newV2=v1;
    for (int i = 1; float(i)<=dy; i++)  {
        float f = float(i) / dy;
        newV1.x = x + (x2 - x) * f; // Calculate the start point
        newV1.y -= float(s); // move y with 1 step each time
        newV1.depth = v1.depth + z1Diff * f;
        newV2.x = x + (x3 - x) * f; // Calculate the end point
        newV2.y -= float(s); // move y with 1 step each time
        newV2.depth = v1.depth + z2Diff * f;
        // v1.depth = - f * (v2.depth - v1.depth)

        // f = v1.depth/(v1.depth-v2.depth)
        drawLine(window, newV1, newV2, colour);
    }

}

void drawFilledTriangle(DrawingWindow &window, CanvasTriangle triangle, const Colour& colour) {
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
    float extraPointY;
    float extraPointDepth;
    // when radio is negative, the extra point is behind on the topPoint (-); when it is positive, the extraPoint is in the front of topPoint (+)
    if(similarTriangleRadio >= 0) {
        extraPointX = float(topPoint.x + ((givenPoint.y - topPoint.y) * abs(similarTriangleRadio)));
    } else {
        extraPointX = float(topPoint.x - ((givenPoint.y - topPoint.y) * abs(similarTriangleRadio)));
    }
    extraPointY = givenPoint.y;
    // calculate extraPoint depth
//    float numberOfSteps = bottomPoint.y - topPoint.y;
//    std::vector<float> depthInterpolation = interpolateSingleFloats(topPoint.depth, bottomPoint.depth, int(numberOfSteps));
//    float index = extraPointY - topPoint.y;
//    extraPointDepth = depthInterpolation[round(index)];
    float alpha = (-(extraPointX-triangle.v1().x)*(triangle.v2().y-triangle.v1().y)+(extraPointY-triangle.v1().y)*(triangle.v2().x-triangle.v1().x))/(-(triangle.v0().x-triangle.v1().x)*(triangle.v2().y-triangle.v1().y)+(triangle.v0().y-triangle.v1().y)*(triangle.v2().x-triangle.v1().x));
    float beta = (-(extraPointX-triangle.v2().x)*(triangle.v0().y-triangle.v2().y)+(extraPointY-triangle.v2().y)*(triangle.v0().x-triangle.v2().x))/(-(triangle.v1().x-triangle.v2().x)*(triangle.v0().y-triangle.v2().y)+(triangle.v1().y-triangle.v2().y)*(triangle.v0().x-triangle.v2().x));
    float gamma = 1 - alpha - beta;
    extraPointDepth = alpha * triangle.v0().depth + beta * triangle.v1().depth + gamma * triangle.v2().depth;

    extraPoint = CanvasPoint(extraPointX, extraPointY, extraPointDepth);

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

// ---------------------- Week 4 ---------------------- //
// Convert mtl colours to packed integer colours
Colour mtlConverter(double r, double g, double b) {
    Colour packedColour = Colour(int(r*255), int(g*255), int(b*255));
    return packedColour;
}

// Read colour from mtl file
std::map<std::string, Colour> mtlReader(const std::string& file) {
    std::map<std::string, Colour> colours;
    // Create a text string
    std::string myText;
    std::string name;
    // Read from the text file
    std::ifstream File(file);
    // Use a while loop together with the getline() function to read the file line by line
    while(getline(File, myText)) {
        std::vector<std::string> text = split(myText, ' ');
        if(text[0] == "newmtl") {
            name = text[1];
        } else if (text[0] == "Kd") {
            colours[name] = mtlConverter(std::stod(text[1]), std::stod(text[2]), std::stod(text[3]));
        }
    }
    // Close the file
    File.close();
    return colours;
}

// Read data from obj file
std::vector<ModelTriangle> objReader(const std::string& objFile, const std::string& mtlFile, float scalingFactor) {
    std::vector<glm::vec3> vertex;
    std::vector<std::vector<std::string>> facets;
    std::vector<ModelTriangle> modelTriangles;
    // Create a text string, which is used to output the text objFile
    std::string myText;
    std::string colourName;
    // Read from the text objFile
    std::ifstream File(objFile);
    // Use a while loop together with the getline() function to read the objFile line by line
    while(getline(File, myText)) {
        std::vector<std::string> text = split(myText, ' ');
        // get colour
        if(text[0] == "usemtl") {
            colourName = text[1];
        // get vertex coordinates
        } else if(text[0] == "v") {
            glm::vec3 v = glm::vec3(std::stod(text[1]), std::stod(text[2]), std::stod(text[3]));
            vertex.push_back(v);
        // get triangle facets
        } else if(text[0] == "f") {
            std::vector<std::string> f {text[1], text[2], text[3], colourName};
            facets.push_back(f);
        }
    }
    // Close the File
    File.close();
    // Get mtl file
    std::map<std::string, Colour> colourMap = mtlReader(mtlFile);
    // Create model triangle list
    for(std::vector<std::string> i : facets) {
        glm::vec3 v0 = vertex[std::stoi(i[0])-1];
        glm::vec3 v1 = vertex[std::stoi(i[1])-1];
        glm::vec3 v2 = vertex[std::stoi(i[2])-1];
        Colour colour = colourMap[i[3]];
        ModelTriangle triangle = ModelTriangle(v0*scalingFactor, v1*scalingFactor, v2*scalingFactor, colour);
        modelTriangles.push_back(triangle);
    }

    return modelTriangles;
}

glm::vec3 cameraCoordinateSystemConverter(glm::vec3 cameraPosition, glm::vec3 vertexPosition) {
    glm::vec3 converted = {vertexPosition.x - cameraPosition.x, vertexPosition.y - cameraPosition.y,  vertexPosition.z - cameraPosition.z};
    return converted;
}

CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::mat3 cameraOrientation, glm::vec3 vertexPosition, float focalLength, float scalingFactor) {
    CanvasPoint twoDPosition;
    glm::vec3 convertedPosition = cameraCoordinateSystemConverter(cameraPosition, vertexPosition);
    glm::vec3 adjustedPosition = convertedPosition * cameraOrientation;
    twoDPosition.x = round(focalLength * (adjustedPosition.x / adjustedPosition.z) * (-scalingFactor) + (float(WIDTH)/2));
    twoDPosition.y = round(focalLength * (adjustedPosition.y / adjustedPosition.z) * scalingFactor + (float(HEIGHT)/2));
    twoDPosition.depth = abs(1/adjustedPosition.z);
    return twoDPosition;
}

CanvasTriangle getCanvasTriangle(ModelTriangle modelTriangle, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, float focalLength, float scalingFactor) {
    CanvasPoint v0 = getCanvasIntersectionPoint(cameraPosition, cameraOrientation, modelTriangle.vertices[0], focalLength, scalingFactor);
    CanvasPoint v1 = getCanvasIntersectionPoint(cameraPosition, cameraOrientation, modelTriangle.vertices[1], focalLength, scalingFactor);
    CanvasPoint v2 = getCanvasIntersectionPoint(cameraPosition, cameraOrientation, modelTriangle.vertices[2], focalLength, scalingFactor);
    CanvasTriangle canvasTriangle = CanvasTriangle(v0, v1, v2);
    return canvasTriangle;
}

void pointCloudRender(DrawingWindow &window, const std::vector<ModelTriangle>& modelTriangles, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, float focalLength, float scalingFactor, const Colour& colour) {
    // Calculate colour
    int red = colour.red;
    int green = colour.green;
    int blue = colour.blue;
    uint32_t c = (255 << 24) + (red << 16) + (green << 8) + (blue);
    // Draw point
    for(const ModelTriangle& modelTriangle : modelTriangles) {
        CanvasTriangle canvasTriangle = getCanvasTriangle(modelTriangle, cameraPosition,  cameraOrientation, focalLength, scalingFactor);
        window.setPixelColour(int(round(canvasTriangle.v0().x)), int(round(canvasTriangle.v0().y)), c);
        window.setPixelColour(int(round(canvasTriangle.v1().x)), int(round(canvasTriangle.v1().y)), c);
        window.setPixelColour(int(round(canvasTriangle.v2().x)), int(round(canvasTriangle.v2().y)), c);
    }
}

void wireframeRender(DrawingWindow &window, const std::vector<ModelTriangle>& modelTriangles, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, float focalLength, float scalingFactor, const Colour& colour) {
    // Draw triangle
    for(const ModelTriangle& modelTriangle : modelTriangles) {
        CanvasTriangle canvasTriangle = getCanvasTriangle(modelTriangle, cameraPosition, cameraOrientation, focalLength, scalingFactor);
        drawStrokedTriangle(window, canvasTriangle, colour);
    }
}

void rasterisedRender(DrawingWindow &window, const std::vector<ModelTriangle>& modelTriangles, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, float focalLength, float scalingFactor) {
    // Rasterise triangle
    for(const ModelTriangle& modelTriangle : modelTriangles) {
        CanvasTriangle canvasTriangle = getCanvasTriangle(modelTriangle, cameraPosition, cameraOrientation, focalLength, scalingFactor);
        drawFilledTriangle(window, canvasTriangle, modelTriangle.colour);
    }
}

// ---------------------- Week 5 ---------------------- //
glm::vec3 moveCamera(SDL_Event event, glm::vec3 cameraPosition) {
    glm::vec3 newCameraPosition = cameraPosition;
    float theta = 0.1;
    float movement = 0.1;
    if (event.type == SDL_KEYDOWN) {
        // Translate
        if (event.key.keysym.sym == SDLK_LEFT) {
            newCameraPosition.x += movement;
        } else if (event.key.keysym.sym == SDLK_RIGHT) {
            newCameraPosition.x -= movement;
        } else if (event.key.keysym.sym == SDLK_UP) {
            newCameraPosition.y += movement;
        } else if (event.key.keysym.sym == SDLK_DOWN) {
            newCameraPosition.y -= movement;
        } else if (event.key.keysym.sym == SDLK_f) {
            newCameraPosition.z -= movement;
        } else if (event.key.keysym.sym == SDLK_b) {
            newCameraPosition.z += movement;
        // Rotate
        } else if (event.key.keysym.sym == SDLK_d) {
            glm::mat3 rotation = glm::mat3(cos(theta), 0, -sin(theta), 0, 1, 0, sin(theta), 0, cos(theta));
            newCameraPosition = rotation * cameraPosition;
        } else if (event.key.keysym.sym == SDLK_a) {
            glm::mat3 rotation = glm::mat3(cos(theta), 0, sin(theta), 0, 1, 0, -sin(theta), 0, cos(theta));
            newCameraPosition = rotation * cameraPosition;
        } else if (event.key.keysym.sym == SDLK_w) {
            glm::mat3 rotation = glm::mat3(1, 0, 0, 0, cos(theta), -sin(theta), 0, sin(theta), cos(theta));
            newCameraPosition = rotation * cameraPosition;
        } else if (event.key.keysym.sym == SDLK_s) {
            glm::mat3 rotation = glm::mat3(1, 0, 0, 0, cos(theta), sin(theta), 0, -sin(theta), cos(theta));
            newCameraPosition = rotation * cameraPosition;
        }
    }
    return newCameraPosition;
}

glm::mat3 rotateCameraOrientation(SDL_Event event, glm::mat3 cameraOrientation) {
    float theta = 0.1;
    glm::mat3 newCameraOrientation;
    if (event.key.keysym.sym == SDLK_d) {
        glm::mat3 rotation = glm::mat3(cos(theta), 0, -sin(theta), 0, 1, 0, sin(theta), 0, cos(theta));
        newCameraOrientation = rotation * cameraOrientation;
    } else if (event.key.keysym.sym == SDLK_a) {
        glm::mat3 rotation = glm::mat3(cos(theta), 0, sin(theta), 0, 1, 0, -sin(theta), 0, cos(theta));
        newCameraOrientation = rotation * cameraOrientation;
    } else if (event.key.keysym.sym == SDLK_w) {
        glm::mat3 rotation = glm::mat3(1, 0, 0, 0, cos(theta), -sin(theta), 0, sin(theta), cos(theta));
        newCameraOrientation = rotation * cameraOrientation;
    } else if (event.key.keysym.sym == SDLK_s) {
        glm::mat3 rotation = glm::mat3(1, 0, 0, 0, cos(theta), sin(theta), 0, -sin(theta), cos(theta));
        newCameraOrientation = rotation * cameraOrientation;
    }
    return newCameraOrientation;
}

glm::vec3 orbit(glm::vec3 cameraPosition, float theta) {
    glm::mat3 rotation = glm::mat3(cos(theta), 0, sin(theta), 0, 1, 0, -sin(theta), 0, cos(theta));
    glm::vec3 newCameraPosition = rotation * cameraPosition;
    return newCameraPosition;
}

glm::mat3 lookAt(glm::vec3 cameraPosition) {
    glm::vec3 vertical (0, 1, 0);
    glm::vec3 center (0, 0, 0);
    glm::vec3 forward = glm::normalize(cameraPosition - center);
    glm::vec3 right = glm::cross(vertical, forward);
    glm::vec3 up = glm::cross(forward, right);
    glm::mat3 newCameraOrientation = glm::mat3(right.x, right.y, right.z, up.x, up.y, up.z, forward.x, forward.y, forward.z);
    return newCameraOrientation;
}

// ---------------------- Week 6 ---------------------- //
glm::vec3 get3DPoint(CanvasPoint point, glm::vec3 cameraPosition, float focalLength, float scalingFactor) {
    float u = point.x;
    float v = point.y;
    float z = cameraPosition.z - focalLength;
    float x = ((u - (float(WIDTH)/2)) / (focalLength * (scalingFactor))) * z;
    float y = ((v - (float(HEIGHT)/2)) / (focalLength * (-scalingFactor))) * z;
    glm::vec3 threeDPoint (x, y, z);
    return threeDPoint;
}

RayTriangleIntersection getClosestIntersection(glm::vec3 source, glm::vec3 rayDirection, const std::vector<ModelTriangle>& modelTriangles) {
    RayTriangleIntersection closestIntersection;
    float closestT = FLT_MAX;
    for(int i = 0; i < int(modelTriangles.size()); i++) {
        glm::vec3 e0 = modelTriangles[i].vertices[1] - modelTriangles[i].vertices[0];
        glm::vec3 e1 = modelTriangles[i].vertices[2] - modelTriangles[i].vertices[0];
        glm::vec3 SPVector = source - modelTriangles[i].vertices[0];
        glm::mat3 DEMatrix(-rayDirection, e0, e1);
        glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
        float u = possibleSolution.y;
        float v = possibleSolution.z;
        if(possibleSolution.x < closestT && possibleSolution.x > 0 && u >= 0 && u <= 1 && v >= 0 && v <= 1 && (u + v) <= 1) {
            closestT = possibleSolution.x;
            glm::vec3 intersectionPoint = source + rayDirection * closestT;
            closestIntersection = RayTriangleIntersection(intersectionPoint, closestT, modelTriangles[i], i);
        }
    }
    return closestIntersection;
}

void rayTrace(DrawingWindow &window, const std::vector<ModelTriangle>& modelTriangles, glm::vec3 cameraPosition, glm::vec3 lightPosition, float focalLength, float scalingFactor) {
    int red;
    int green;
    int blue;
    for(int y = 0; y < HEIGHT; y++) {
        for(int x = 0; x < WIDTH; x++) {
            CanvasPoint point = CanvasPoint(float(x), float(y));
            glm::vec3 threeDPoint = get3DPoint(point, cameraPosition, focalLength, scalingFactor);
            glm::vec3 rayDirection = threeDPoint - cameraPosition;
            RayTriangleIntersection closestIntersection = getClosestIntersection(cameraPosition, rayDirection, modelTriangles);
            glm::vec3 lightDirection = closestIntersection.intersectionPoint - lightPosition;
            RayTriangleIntersection lightIntersection = getClosestIntersection(lightPosition, lightDirection, modelTriangles);
            if (closestIntersection.triangleIndex != lightIntersection.triangleIndex) {
                red = 0;
                green = 0;
                blue = 0;
            } else {
                Colour colour = closestIntersection.intersectedTriangle.colour;
                red = colour.red;
                green = colour.green;
                blue = colour.blue;
            }

            uint32_t c = (255 << 24) + (red << 16) + (green << 8) + (blue);
            window.setPixelColour(x, y, c);
        }
    }
}

// ---------------------- Box model ---------------------- //
glm::vec3 cameraPosition (0.0, 0.0, 4.0);
glm::vec3 lightPosition (0.0, 0.8, 0.5);
glm::mat3 cameraOrientation (1, 0, 0,0,1, 0 ,0, 0, 1);
std::vector<ModelTriangle> modelTriangles = objReader("cornell-box.obj", "cornell-box.mtl", 0.35);
float focalLength = 2.0;
bool orbitYes = false;
bool rasterisedRenderYes = false;

// ---------------------- Show in window ---------------------- //
void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) {
            std::cout << "LEFT" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
        } else if (event.key.keysym.sym == SDLK_RIGHT) {
            std::cout << "RIGHT" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
        } else if (event.key.keysym.sym == SDLK_UP) {
            std::cout << "UP" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
        } else if (event.key.keysym.sym == SDLK_DOWN) {
            std::cout << "DOWN" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
        } else if (event.key.keysym.sym == SDLK_f) {
            std::cout << "FORWARD" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
        } else if (event.key.keysym.sym == SDLK_b) {
            std::cout << "BACKWARD" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
        } else if (event.key.keysym.sym == SDLK_d) {
            std::cout << "Rotate about X-axis Clockwise" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
            cameraOrientation = rotateCameraOrientation(event, cameraOrientation);
            //cameraOrientation = lookAt(cameraPosition);
        } else if (event.key.keysym.sym == SDLK_a) {
            std::cout << "Rotate about X-axis Anticlockwise" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
            cameraOrientation = rotateCameraOrientation(event, cameraOrientation);
            //cameraOrientation = lookAt(cameraPosition);
        } else if (event.key.keysym.sym == SDLK_w) {
            std::cout << "Rotate about Y-axis Clockwise" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
            cameraOrientation = rotateCameraOrientation(event, cameraOrientation);
            //cameraOrientation = lookAt(cameraPosition);
        } else if (event.key.keysym.sym == SDLK_s) {
            std::cout << "Rotate about Y-axis Anticlockwise" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
            cameraOrientation = rotateCameraOrientation(event, cameraOrientation);
            //cameraOrientation = lookAt(cameraPosition);
        } else if (event.key.keysym.sym == SDLK_j) {
            clearDepthBuffer();
            window.clearPixels();
            rayTrace(window, modelTriangles, cameraPosition, lightPosition, focalLength, float(HEIGHT)*2/3);
        } else if (event.key.keysym.sym == SDLK_u) {
            CanvasPoint point1,point2,point3;
            CanvasTriangle triangle;
            point1 = CanvasPoint(float(rand()%window.width), float(rand()%window.height));
            point2 = CanvasPoint(float(rand()%window.width), float(rand()%window.height));
            point3 = CanvasPoint(float(rand()%window.width), float(rand()%window.height));
            triangle = CanvasTriangle(point1, point2, point3);
            drawStrokedTriangle(window, triangle, Colour(rand()%256, rand()%256, rand()%256));
        } else if (event.key.keysym.sym == SDLK_t) {
            CanvasPoint point1,point2,point3;
            CanvasTriangle triangle;
            point1 = CanvasPoint(float(rand()%window.width), float(rand()%window.height));
            point2 = CanvasPoint(float(rand()%window.width), float(rand()%window.height));
            point3 = CanvasPoint(float(rand()%window.width), float(rand()%window.height));
            triangle = CanvasTriangle(point1, point2, point3);
            drawStrokedTriangle(window, triangle, Colour(255, 255, 255));
            drawFilledTriangle(window, triangle, Colour(rand()%256, rand()%256, rand()%256));
        } else if(event.key.keysym.sym == SDLK_o) {
            CanvasPoint point1, point2, point3;
            point1 = CanvasPoint(160, 10);
            point2 = CanvasPoint(300, 230);
            point3 = CanvasPoint(10, 150);
            point1.texturePoint.x = 195; point1.texturePoint.y = 5;
            point2.texturePoint.x = 395; point2.texturePoint.y = 380;
            point3.texturePoint.x = 65; point3.texturePoint.y = 330;
            drawTextureTriangle(window, getTextureMap("texture.ppm"), CanvasTriangle(point1, point2, point3));
        } else if(event.key.keysym.sym == SDLK_p) {
            pointCloudRender(window, modelTriangles, cameraPosition, cameraOrientation, focalLength, float(HEIGHT)*2/3, Colour(255, 255, 255));
        } else if(event.key.keysym.sym == SDLK_l) {
            wireframeRender(window, modelTriangles, cameraPosition, cameraOrientation, focalLength, float(HEIGHT)*2/3, Colour(255, 255, 255));
        } else if(event.key.keysym.sym == SDLK_r) {
            if(rasterisedRenderYes) {
                // reset
                clearDepthBuffer();
                window.clearPixels();
                cameraPosition = glm::vec3(0.0, 0.0, 4.0);
                cameraOrientation  = glm::mat3(1, 0, 0,0,1, 0 ,0, 0, 1);
                focalLength = 2.0;
                orbitYes = false;
                rasterisedRenderYes = false;
            } else {
                clearDepthBuffer();
                window.clearPixels();
                rasterisedRenderYes = true;
            }
            //rasterisedRender(window, modelTriangles, cameraPosition, cameraOrientation, focalLength, float(HEIGHT)*2/3);
        } else if(event.key.keysym.sym == SDLK_c) {
            // reset the matrix
            clearDepthBuffer();
            window.clearPixels();
        } else if(event.key.keysym.sym == SDLK_q) {
            if(orbitYes) {
                orbitYes = false;
            } else {
                orbitYes = true;
            }
        } else if(event.key.keysym.sym == SDLK_SPACE) {
            orbitYes = false;
        }
    } else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
    clearDepthBuffer();
    // week 1
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

    // week 4
    // Print modelTriangle list
//    std::vector<ModelTriangle> modelTriangles = objReader("cornell-box.obj", "cornell-box.mtl", 1);
//    for(const ModelTriangle& t : modelTriangles) {std::cout << t << t.colour << std::endl;}

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
        if (rasterisedRenderYes) {
            if (orbitYes) {
                window.clearPixels();
                clearDepthBuffer();
                cameraOrientation = lookAt(cameraPosition);
                cameraPosition = orbit(cameraPosition, 0.005);

            }
            rasterisedRender(window, modelTriangles, cameraPosition, cameraOrientation, focalLength, float(HEIGHT)*2/3);
        }
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
