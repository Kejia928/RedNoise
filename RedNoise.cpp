#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <CanvasPoint.h>
#include <ModelTriangle.h>
#include <TextureMap.h>
#include <Colour.h>
#include <Utils.h>
#include <utility>
#include <vector>
#include <map>
#include <glm/glm.hpp>
#include <RayTriangleIntersection.h>


#define WIDTH 320
#define HEIGHT 240
#define PI 3.1415926f

enum shading {
    FLAT, GOURAND, PHONG
};

enum shadow {
    SOFT, HARD, NONE
};

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

uint32_t modelTextureMapper(float alpha, float beta, float gamma, ModelTriangle triangle, TextureMap textureMap){
    // Get texture point
    CanvasPoint texturePoint(round(alpha * triangle.texturePoints[0].x * float(textureMap.width) + beta * triangle.texturePoints[1].x * float(textureMap.width) + gamma * triangle.texturePoints[2].x * float(textureMap.width)), round(alpha * triangle.texturePoints[0].y * float(textureMap.height) + beta * triangle.texturePoints[1].y * float(textureMap.height) + gamma * triangle.texturePoints[2].y * float(textureMap.height)));
    // Locate the texture point position in texture map pixel list
    int index = int(texturePoint.y) * int(textureMap.width) + int(texturePoint.x);
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
glm::vec3 calculateSphereCentre(const glm::vec3 p1, const glm::vec3 p2, const glm::vec3 p3, const glm::vec3 p4) {
    std::cout << p1.x << " " << p1.y << " " << p1.z << std::endl;
    std::cout << p2.x << " " << p2.y << " " << p2.z << std::endl;
    std::cout << p3.x << " " << p3.y << " " << p3.z << std::endl;
    std::cout << p4.x << " " << p4.y << " " << p4.z << std::endl;

    glm::vec3 centre(0, 0, 0);
    float a = p1[0] - p2[0], b = p1[1] - p2[1], c = p1[2] - p2[2];
    float a1 = p3[0] - p4[0], b1 = p3[1] - p4[1], c1 = p3[2] - p3[2];
    float a2 = p2[0] - p3[0], b2 = p2[1] - p3[1], c2 = p2[2] - p3[2];
    float A = p1[0] * p1[0] - p2[0] * p2[0];
    float B = p1[1] * p1[1] - p2[1] * p2[1];
    float C = p1[2] * p1[2] - p2[2] * p2[2];
    float A1 = p3[0] * p3[0] - p4[0] * p4[0];
    float B1 = p3[1] * p3[1] - p4[1] * p4[1];
    float C1 = p3[2] * p3[2] - p4[2] * p4[2];
    float A2 = p2[0] * p2[0] - p3[0] * p3[0];
    float B2 = p2[1] * p2[1] - p3[1] * p3[1];
    float C2 = p2[2] * p2[2] - p3[2] * p3[2];
    float P = (A + B + C) / 2;
    float Q = (A1 + B1 + C1) / 2;
    float R = (A2 + B2 + C2) / 2;

    float D = a * b1 * c2 + a2 * b * c1 + c * a1 * b2 - (a2 * b1 * c + a1 * b * c2 + a * b2 * c1);
    float Dx = P * b1 * c2 + b * c1 * R + c * Q * b2 - (c * b1 * R + P * c1 * b2 + Q * b * c2);
    float Dy = a * Q * c2 + P * c1 * a2 + c * a1 * R - (c * Q * a2 + a * c1 * R + c2 * P * a1);
    float Dz = a * b1 * R + b * Q * a2 + P * a1 * b2 - (a2 * b1 * P + a * Q * b2 + R * b * a1);

    if (D == 0) {
        return centre;
    } else {
        centre.x = Dx / D;
        centre.y = Dy / D;
        centre.z = Dz / D;
        return centre;
    }
}

float calculateRadius(glm::vec3 p1, glm::vec3 centre){
    return sqrt((p1[0]-centre[0])*(p1[0]-centre[0]) +
                (p1[1]-centre[1])*(p1[1]-centre[1]) +
                (p1[2]-centre[2])*(p1[2]-centre[2]));
}

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
    std::vector<TexturePoint> textureVertex;
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
        } else if(text[0] == "vt"){
            std::cout << text[3] <<std::endl;
            TexturePoint vt = TexturePoint(std::stof(text[1]), std::stof(text[2]));
            textureVertex.push_back(vt);
        } else if(text[0] == "f") {
            if(colourName == "Cobbles") {
                std::vector<std::string> faceVertex;
                std::vector<std::string> faceTextureVertex;
                for(int i = 1; i < 4; i++) {
                    std::vector<std::string> t = split(text[i], '/');
                    faceVertex.push_back(t[0]);
                    faceTextureVertex.push_back(t[1]);
                }
                std::vector<std::string> f {faceVertex[0], faceVertex[1], faceVertex[2], colourName, faceTextureVertex[0], faceTextureVertex[1], faceTextureVertex[2]};
                facets.push_back(f);
            } else {
                std::vector<std::string> f {text[1], text[2], text[3], colourName};
                facets.push_back(f);
            }

        }
    }
    // Close the File
    File.close();
    // Get mtl file
    std::map<std::string, Colour> colourMap;
    if(mtlFile != "null") {
        colourMap = mtlReader(mtlFile);
    }

    // Create model triangle list
    for(std::vector<std::string> i : facets) {
        ModelTriangle triangle;
        glm::vec3 v0 = vertex[std::stoi(i[0])-1];
        glm::vec3 v1 = vertex[std::stoi(i[1])-1];
        glm::vec3 v2 = vertex[std::stoi(i[2])-1];
        triangle.vertices = {v0*scalingFactor, v1*scalingFactor, v2*scalingFactor};
        if(mtlFile != "null") {
            triangle.colour = colourMap[i[3]];
            triangle.colour.name = i[3];
            if(triangle.colour.name == "Cobbles") {
                TexturePoint vt0 = textureVertex[std::stoi(i[4])-1];
                TexturePoint vt1 = textureVertex[std::stoi(i[5])-1];
                TexturePoint vt2 = textureVertex[std::stoi(i[6])-1];
                triangle.texturePoints = {vt0, vt1, vt2};
            }
        } else {
            triangle.colour = Colour(255, 165, 0);
            triangle.colour.name = "Orange";
        }

        // normal
        triangle.normal = glm::normalize(glm::cross((triangle.vertices[2] - triangle.vertices[0]), (triangle.vertices[1] - triangle.vertices[0])));
        modelTriangles.push_back(triangle);
    }

    return modelTriangles;
}

glm::vec3 sphereCenter (0, 0, 0);
TextureMap Lunar = TextureMap("Venus.ppm");
TextureMap Logo = TextureMap("texture_cw.ppm");
float sphereRadius = 0;

std::vector<ModelTriangle> objReader2(const std::vector<std::string>& objFiles, const std::vector<std::string>& mtlFiles, float scalingFactor) {
    std::vector<ModelTriangle> modelTriangles;
    std::vector<glm::vec3> sphereVertex;
    for(const std::string& objFile : objFiles) {
        for(const std::string& mtlFile : mtlFiles) {
            std::vector<glm::vec3> vertex;
            std::vector<std::vector<std::string>> facets;
            std::vector<TexturePoint> textureVertex;
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
                } else if(text[0] == "vt"){
                    std::cout << text[3] <<std::endl;
                    TexturePoint vt = TexturePoint(std::stof(text[1]), std::stof(text[2]));
                    textureVertex.push_back(vt);
                } else if(text[0] == "f") {
                    if(colourName == "Cobbles" || objFile == "logo.obj") {
                        std::vector<std::string> faceVertex;
                        std::vector<std::string> faceTextureVertex;
                        for(int i = 1; i < 4; i++) {
                            std::vector<std::string> t = split(text[i], '/');
                            faceVertex.push_back(t[0]);
                            faceTextureVertex.push_back(t[1]);
                        }
                        std::vector<std::string> f {faceVertex[0], faceVertex[1], faceVertex[2], colourName, faceTextureVertex[0], faceTextureVertex[1], faceTextureVertex[2]};
                        facets.push_back(f);
                    } else {
                        std::vector<std::string> f {text[1], text[2], text[3], colourName};
                        facets.push_back(f);
                    }

                }
            }
            // Close the File
            File.close();
            // Get mtl file
            std::map<std::string, Colour> colourMap;
            if(mtlFile != "null") {
                colourMap = mtlReader(mtlFile);
            }

            // Create model triangle list
            for(std::vector<std::string> i : facets) {
                ModelTriangle triangle;
                glm::vec3 v0 = vertex[std::stoi(i[0])-1];
                glm::vec3 v1 = vertex[std::stoi(i[1])-1];
                glm::vec3 v2 = vertex[std::stoi(i[2])-1];
                if(mtlFile != "null") {
                    triangle.colour = colourMap[i[3]];
                    triangle.colour.name = i[3];
                    if(triangle.colour.name == "Cobbles" || objFile == "logo.obj") {
                        TexturePoint vt0 = textureVertex[std::stoi(i[4])-1];
                        TexturePoint vt1 = textureVertex[std::stoi(i[5])-1];
                        TexturePoint vt2 = textureVertex[std::stoi(i[6])-1];
                        triangle.texturePoints = {vt0, vt1, vt2};
                    }
                }
                if(objFile == "sphere.obj") {
                    glm::vec3 movement {0.5, 1.134, 0};
                    //glm::vec3 movement {-1, 1.6, 0};
                    triangle.colour = Colour(255, 165, 0);
                    triangle.colour.name = "Orange";
                    v0 = v0*scalingFactor - movement;
                    v1 = v1*scalingFactor - movement;
                    v2 = v2*scalingFactor - movement;
                    sphereVertex.push_back(v0);
                    sphereVertex.push_back(v1);
                    sphereVertex.push_back(v2);
                } else if(objFile == "logo.obj") {
                    glm::vec3 movement {-0.2, 0, 1};
                    triangle.colour = Colour(128, 64, 0);
                    triangle.colour.name = "Brown";
                    v0 = (v0/250.0f)*scalingFactor - movement;
                    v1 = (v1/250.0f)*scalingFactor - movement;
                    v2 = (v2/250.0f)*scalingFactor - movement;
                } else {
                    v0 = v0*scalingFactor;
                    v1 = v1*scalingFactor;
                    v2 = v2*scalingFactor;
                }

                triangle.vertices = {v0, v1, v2};
                // normal
                triangle.normal = glm::normalize(glm::cross((triangle.vertices[2] - triangle.vertices[0]), (triangle.vertices[1] - triangle.vertices[0])));
                modelTriangles.push_back(triangle);
            }
            break;
        }
    }
    // calculate sphere sphereCenter and radius
    sphereCenter = calculateSphereCentre(sphereVertex[3], sphereVertex[10], sphereVertex[20], sphereVertex[30]);
    sphereRadius = calculateRadius(sphereVertex[0], sphereCenter);
    std::cout << sphereCenter.x << " " << sphereCenter.y << " " << sphereCenter.z << std::endl;
    std::cout << sphereRadius << std::endl;
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
    //glm::mat3 rotation = glm::mat3(1, 0, 0, 0, cos(theta), sin(theta), 0, -sin(theta), cos(theta));
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
//std::vector<glm::vec3> multiLightPosition(glm::vec3 sphereCenter) {
//    std::vector<glm::vec3> lightPositions;
//    float i = -0.2f;
//    for(;;) {
//        if(i > 0.2f) { break;}
//        float j = -0.2f;
//        for(;;) {
//            if(j > 0.2f) { break;}
//            lightPositions.emplace_back(sphereCenter.x+i, sphereCenter.y, sphereCenter.z+j);
//            j = j + 0.02f;
//        }
//        i = i + 0.02f;
//    }
//    return lightPositions;
//}
//
std::vector<glm::vec3> multiLightPosition(glm::vec3 lightCenterPoint) {
    int radius = 10;
    float spread = 0.05;
    std::vector<glm::vec3> lightPoints;
    for (int i = -radius; i <= radius; i++) {
        for (int j = -radius; j <= radius; j++) {
            if (std::abs(j) + std::abs(i) <= radius)
                lightPoints.push_back(glm::vec3(lightCenterPoint.x + float(i) * spread, lightCenterPoint.y, lightCenterPoint.z + float(j) * spread));
        }
    }
    return lightPoints;
}

glm::vec3 get3DPoint(CanvasPoint point, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, float focalLength, float scalingFactor) {
    float u = point.x;
    float v = point.y;
    // The focal length is the distance between the window and camera
    float z = cameraPosition.z - focalLength;
    float x = ((u - (float(WIDTH)/2)) / ((scalingFactor)));
    float y = ((v - (float(HEIGHT)/2)) / ((-scalingFactor)));
    glm::vec3 threeDPoint (x, y, z);
    // Using cameraOrientation to correct the 3D point
//    threeDPoint = cameraOrientation * threeDPoint;
    threeDPoint = {threeDPoint.x + cameraPosition.x, threeDPoint.y + cameraPosition.y, threeDPoint.z};
    // std::cout << threeDPoint.x << " " << threeDPoint.y << " " << threeDPoint.z << std::endl;
    return threeDPoint;
}

RayTriangleIntersection getClosestIntersection(glm::vec3 source, glm::vec3 rayDirection, const std::vector<ModelTriangle>& modelTriangles) {
    RayTriangleIntersection closestIntersection;
    // Set the initial distance is as large as possible
    closestIntersection.distanceFromCamera = FLT_MAX;
    for(int i = 0; i < int(modelTriangles.size()); i++) {
        glm::vec3 e0 = modelTriangles[i].vertices[1] - modelTriangles[i].vertices[0];
        glm::vec3 e1 = modelTriangles[i].vertices[2] - modelTriangles[i].vertices[0];
        glm::vec3 SPVector = source - modelTriangles[i].vertices[0];
        glm::mat3 DEMatrix(-rayDirection, e0, e1);
        glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
        float u = possibleSolution.y;
        float v = possibleSolution.z;
        if(possibleSolution.x <= closestIntersection.distanceFromCamera && possibleSolution.x >= 0.000001 && u >= 0.000001 && u <= 1 && v >= 0.000001 && v <= 1 && (u + v) <= 1) {
            closestIntersection.distanceFromCamera = possibleSolution.x;
            glm::vec3 intersectionPoint = source + rayDirection * closestIntersection.distanceFromCamera;
            closestIntersection = RayTriangleIntersection(intersectionPoint, closestIntersection.distanceFromCamera, modelTriangles[i], i);
        }
    }
    return closestIntersection;
}

// ---------------------- Week 7 ---------------------- //
glm::vec3 moveLight(SDL_Event event, glm::vec3 lightPosition) {
    glm::vec3 newLightPosition = lightPosition;
    float movement = 0.1;
    if (event.type == SDL_KEYDOWN) {
        // Translate
        if (event.key.keysym.sym == SDLK_j) {
            newLightPosition.x += movement;
        } else if (event.key.keysym.sym == SDLK_l) {
            newLightPosition.x -= movement;
        } else if (event.key.keysym.sym == SDLK_i) {
            newLightPosition.y += movement;
        } else if (event.key.keysym.sym == SDLK_k) {
            newLightPosition.y -= movement;
        } else if (event.key.keysym.sym == SDLK_m) {
            newLightPosition.z -= movement;
        } else if (event.key.keysym.sym == SDLK_n) {
            newLightPosition.z += movement;
        }
    }
    return newLightPosition;
}

// calculate the vertex normal vector by taking the average neighbor face normal vector
glm::vec3 vertexNormalCalculator(glm::vec3 vertex, const std::vector<ModelTriangle>& sphereModel) {
    glm::vec3 vertexNormal;
    float faceNumber = 0.0f;
    for (ModelTriangle triangle : sphereModel) {
        if (triangle.vertices[0] == vertex || triangle.vertices[1] == vertex || triangle.vertices[2] == vertex) {
            faceNumber++;
            vertexNormal += triangle.normal;
        }
    }
    // The vertex normal vector is equal to the average neighbor face normal vector
    vertexNormal = vertexNormal * (1 / faceNumber);
    return glm::normalize(vertexNormal);
}

glm::vec3 getRefractionDirection(float air, float glass, glm::vec3 incidenceRay, glm::vec3 normal) {
    glm::vec3 direction (0, 0, 0);
    auto incidenceAngle = glm::dot(incidenceRay, normal);
    float eta = 0.0f;
    float k = 0.0f;
    if(incidenceAngle > 0) {
        eta = air / glass;
        normal = -normal;
    } else {
        eta = glass / air;
    }
    k = 1 - eta * eta * (1 - abs(incidenceAngle) * abs(incidenceAngle));
    if(k < 0){
        return direction;
    }
    direction = glm::normalize(incidenceRay*eta + normal * (eta * abs(incidenceAngle) - sqrt(k)));
    return direction;
}

glm::vec4 refractionDirection(glm::vec3 incident, glm::vec3 normal, float refractiveIndex) {
    float incidentAngle = dot(incident, normal);
    float v1 = 1, v2 = refractiveIndex;
    float c = v1 / v2;
    if (incidentAngle > 0){
        c = 1 / c;
        normal = -normal;
    }

    float cosTheta = abs(incidentAngle);

    float k = 1 - c * c * (1 - cosTheta * cosTheta);
    if (k < 0) return {0,0,0,0};
    glm::vec3 refracted = glm::normalize(incident * c + normal * (c * cosTheta - sqrtf(k)));

    return {refracted[0], refracted[1], refracted[2], sign(incidentAngle)};
}

float fresnel(glm::vec3 incident, glm::vec3 normal, float refractiveIndex) {
    float cosi = glm::dot(incident, normal);
    float etai = 1;
    float etat = refractiveIndex;

    if (cosi > 0){
        std::swap(etai, etat);
    }
    float eta = etai / etat;
    float sint = eta * sqrtf(std::max(0.f, 1 - cosi * cosi));

    if (sint >= 1){
        return 1;
    } else {
        float cost = sqrtf(std::max(0.f, 1 - sint * sint));
        cosi = fabsf(cosi);
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        return (Rs * Rs + Rp * Rp) * 0.5f;
    }
}

Colour castRay(glm::vec3 cameraPosition, glm::vec3 lightPosition, const std::vector<glm::vec3>& lightPositions, glm::vec3 rayDirection, const std::vector<ModelTriangle>& modelTriangles, float lightPower, float ambient, const TextureMap& textureMap, int index, shadow shadowType, shading shadeType) {
    if(index >= 5) {
        return {255, 255, 255};
    }
    // Get the closest intersect model with ray
    RayTriangleIntersection closestIntersection = getClosestIntersection(cameraPosition, rayDirection, modelTriangles);
    glm::vec3 point = closestIntersection.intersectionPoint;
    glm::vec3 v0 = closestIntersection.intersectedTriangle.vertices[0];
    glm::vec3 v1 = closestIntersection.intersectedTriangle.vertices[1];
    glm::vec3 v2 = closestIntersection.intersectedTriangle.vertices[2];
    // Shoot a ray from lightPosition to intersected point
    glm::vec3 lightDirection = glm::normalize(point - lightPosition);
    // Get the closest intersect model with light
    RayTriangleIntersection lightIntersection = getClosestIntersection(lightPosition, lightDirection, modelTriangles);
    glm::vec3 view = glm::normalize(point - cameraPosition);

    // Draw colour
    Colour colour = closestIntersection.intersectedTriangle.colour;

    if(closestIntersection.intersectedTriangle.colour.name == "Magenta") {
        glm::vec3 reflectionRay = glm::normalize(rayDirection - (2.0f * closestIntersection.intersectedTriangle.normal * glm::dot(rayDirection, closestIntersection.intersectedTriangle.normal)));
        return castRay(closestIntersection.intersectionPoint, lightPosition, lightPositions, reflectionRay, modelTriangles, lightPower, ambient, textureMap, index+1,shadowType, shadeType);
    } else if(closestIntersection.intersectedTriangle.colour.name == "Cobbles") {
        float alpha = ((-(point.x - v1.x) * (v2.y - v1.y) + (point.y - v1.y) * (v2.x - v1.x)) /
                       (-(v0.x - v1.x) * (v2.y - v1.y) + (v0.y - v1.y) * (v2.x - v1.x))) / abs(v0.z);
        float beta = ((-(point.x - v2.x) * (v0.y - v2.y) + (point.y - v2.y) * (v0.x - v2.x)) /
                      (-(v1.x - v2.x) * (v0.y - v2.y) + (v1.y - v2.y) * (v0.x - v2.x))) / abs(v1.z);
        float gamma = (1 - alpha - beta) / abs(v2.z);
        uint32_t c = modelTextureMapper(alpha, beta, gamma, closestIntersection.intersectedTriangle, textureMap);
        colour.blue = int((c) & 0xff);
        colour.green = int((c >> 8) & 0xff);
        colour.red = int((c >> 16) & 0xff);
    } else if(closestIntersection.intersectedTriangle.colour.name == "Red") {
        //glm::vec3 refractDirection = getRefractionDirection(1.0f, 1.3f, rayDirection, closestIntersection.intersectedTriangle.normal);
        glm::vec4 refractionD = refractionDirection(rayDirection, closestIntersection.intersectedTriangle.normal, 1.3f);
        glm::vec3 refractDirection (refractionD[0], refractionD[1], refractionD[2]);
        glm::vec3 reflectionRay = glm::normalize(rayDirection - (2.0f * closestIntersection.intersectedTriangle.normal * glm::dot(rayDirection, closestIntersection.intersectedTriangle.normal)));
        Colour reflectionColour = castRay(closestIntersection.intersectionPoint, lightPosition, lightPositions, reflectionRay, modelTriangles, lightPower, ambient, textureMap, index+1, shadowType, shadeType);
        glm::vec3 updatePoint = closestIntersection.intersectionPoint;
        if(refractionD[3] == -1){
            updatePoint = point - (closestIntersection.intersectedTriangle.normal * float(0.0001));
        } else {
            updatePoint = point + (closestIntersection.intersectedTriangle.normal * float(0.0001));
        }
        Colour refractionColour = castRay(updatePoint, lightPosition, lightPositions, refractDirection, modelTriangles, lightPower, ambient, textureMap, index+1, shadowType, shadeType);
        float reflectiveConstant = fresnel(rayDirection, closestIntersection.intersectedTriangle.normal, 1.3f);
        float refractiveConstant = 1 - reflectiveConstant;

        Colour mixColour;
        mixColour.red = int(reflectiveConstant * float(reflectionColour.red) + (refractiveConstant * float(refractionColour.red)));
        mixColour.green = int(reflectiveConstant * float(reflectionColour.green) + (refractiveConstant * float(refractionColour.green)));
        mixColour.blue = int(reflectiveConstant * float(reflectionColour.blue) + (refractiveConstant * float(refractionColour.blue)));
        return mixColour;
    } else if(closestIntersection.intersectedTriangle.colour.name == "Orange") {
        glm::vec3 spherePoint = (point - sphereCenter) / sphereRadius;
        float phi = atan2(spherePoint.z, spherePoint.x);
        float theta = asin(spherePoint.y);
        float x = 1.0f - ((phi+PI) / (2*PI));
        float y = (theta + PI/2) / PI;
        int i = std::max(std::min(int((x * float(Lunar.width))-1), int(Lunar.width-1)), 0);
        int j = std::max(std::min((int((1-y) * float(Lunar.height) - 0.001f)), int(Lunar.height-1)), 0);
        uint32_t c = Lunar.pixels[i + j * int(Lunar.width)];
        colour.blue = int((c) & 0xff);
        colour.green = int((c >> 8) & 0xff);
        colour.red = int((c >> 16) & 0xff);
    } else if(closestIntersection.intersectedTriangle.colour.name == "Brown") {
        // 2D Texture
        float alpha = ((-(point.x - v1.x) * (v2.y - v1.y) + (point.y - v1.y) * (v2.x - v1.x)) /
                       (-(v0.x - v1.x) * (v2.y - v1.y) + (v0.y - v1.y) * (v2.x - v1.x))) / abs(v0.z);
        float beta = ((-(point.x - v2.x) * (v0.y - v2.y) + (point.y - v2.y) * (v0.x - v2.x)) /
                      (-(v1.x - v2.x) * (v0.y - v2.y) + (v1.y - v2.y) * (v0.x - v2.x))) / abs(v1.z);
        float gamma = (1 - alpha - beta) / abs(v2.z);
        uint32_t c = modelTextureMapper(alpha, beta, gamma, closestIntersection.intersectedTriangle, Logo);
        colour.blue = int((c) & 0xff);
        colour.green = int((c >> 8) & 0xff);
        colour.red = int((c >> 16) & 0xff);
    }

    // Calculate brightness
    float brightness = 0;
    if(shadeType == GOURAND || shadeType == PHONG) {
        // Calculate Barycentric coordinates
        float v0X = v0.x;
        float v0Y = v0.y;
        float v0Z = v0.z;
        float v1X = v1.x;
        float v1Y = v1.y;
        float v1Z = v1.z;
        float v2X = v2.x;
        float v2Y = v2.y;
        float v2Z = v2.z;
        float pointX = point.x;
        float pointY = point.y;
        float pointZ = point.z;

        if (v0X == v1X && v1X == v2X) {
            v0X = v0Z;
            v1X = v1Z;
            v2X = v2Z;
            pointX = pointZ;
        } else if(v0Y == v1Y && v1Y == v2Y) {
            v0Y = v0Z;
            v1Y = v1Z;
            v2Y = v2Z;
            pointX = pointZ;
        }
        float alpha = (-(pointX-v1X)*(v2Y-v1Y)+(pointY-v1Y)*(v2X-v1X))/(-(v0X-v1X)*(v2Y-v1Y)+(v0Y-v1Y)*(v2X-v1X));
        float beta = (-(pointX-v2X)*(v0Y-v2Y)+(pointY-v2Y)*(v0X-v2X))/(-(v1X-v2X)*(v0Y-v2Y)+(v1Y-v2Y)*(v0X-v2X));
        float gamma = 1 - alpha - beta;

        // Calculate the vertex normal
        glm::vec3 vn0 = vertexNormalCalculator(v0, modelTriangles);
        glm::vec3 vn1 = vertexNormalCalculator(v1, modelTriangles);
        glm::vec3 vn2 = vertexNormalCalculator(v2, modelTriangles);

        if(shadeType == GOURAND){
            glm::vec3 lightDirection0 = glm::normalize(v0 - lightPosition);
            glm::vec3 lightDirection1 = glm::normalize(v1 - lightPosition);
            glm::vec3 lightDirection2 = glm::normalize(v2 - lightPosition);

            // Proximity Lighting
            glm::vec3 v0ToLight = lightPosition - v0;
            glm::vec3 v1ToLight = lightPosition - v1;
            glm::vec3 v2ToLight = lightPosition - v2;
            float distance0 = glm::length(v0ToLight);
            float distance1 = glm::length(v1ToLight);
            float distance2 = glm::length(v2ToLight);
            float pL0 = lightPower / (4 * PI * distance0 * distance0);
            float pL1 = lightPower / (4 * PI * distance1 * distance1);
            float pL2 = lightPower / (4 * PI * distance2 * distance2);
            // Angle of Incidence Lighting
            float incidenceAngle0 = std::max(0.0f, glm::dot(lightDirection0, vn0));
            float incidenceAngle1 = std::max(0.0f, glm::dot(lightDirection1, vn1));
            float incidenceAngle2 = std::max(0.0f, glm::dot(lightDirection2, vn2));
            pL0 = pL0 * incidenceAngle0;
            pL1 = pL1 * incidenceAngle1;
            pL2 = pL2 * incidenceAngle2;

            // Specular Lighting
            float glossy = 256;
            glm::vec3 view0 = glm::normalize(v0 - cameraPosition);
            glm::vec3 view1 = glm::normalize(v1 - cameraPosition);
            glm::vec3 view2 = glm::normalize(v2 - cameraPosition);
            glm::vec3 reflectionVector0 = lightDirection0 - (2.0f * vn0 * glm::dot(vn0, lightDirection0));
            glm::vec3 reflectionVector1 = lightDirection1 - (2.0f * vn1 * glm::dot(vn1, lightDirection1));
            glm::vec3 reflectionVector2 = lightDirection2 - (2.0f * vn2 * glm::dot(vn2, lightDirection2));
            float sL0 = glm::pow(glm::dot(reflectionVector0, view0), glossy);
            float sL1 = glm::pow(glm::dot(reflectionVector1, view1), glossy);
            float sL2 = glm::pow(glm::dot(reflectionVector2, view2), glossy);

            // Calculate pixel brightness
            // Add all type of light
            float brightness0 = pL0 + sL0 + ambient;
            float brightness1 = pL1 + sL1 + ambient;
            float brightness2 = pL2 + sL2 + ambient;
            brightness = alpha*brightness0 + beta*brightness1 + gamma*brightness2;
        } else if(shadeType == PHONG){
            // Interpolation for normal vector
            glm::vec3 normal = alpha*vn0 + beta*vn1 + gamma*vn2;

            // Proximity Lighting
            glm::vec3 pointToLight = lightPosition - point;
            float distance = glm::length(pointToLight);
            float pL = lightPower / (4 * PI * distance * distance);
            // Angle of Incidence Lighting
            float incidenceAngle = std::max(0.0f, glm::dot(lightDirection, normal));
            pL = pL * incidenceAngle;

            // Specular Lighting
            float glossy = 256;
            glm::vec3 reflectionVector = lightDirection - (2.0f * normal * glm::dot(normal, lightDirection));
            float sL = glm::pow(glm::dot(reflectionVector, view), glossy);

            brightness = pL + sL + ambient;
        }
    } else {
        // Proximity Lighting
        glm::vec3 pointToLight = lightPosition - closestIntersection.intersectionPoint;
        float distance = glm::length(pointToLight);
        float pL = lightPower / (4 * PI * distance * distance);
        // Angle of Incidence Lighting
        // If the angle is less than 0, the PL will be zero, which means there do not have proximity light in that point
        float incidenceAngle = std::max(0.0f, glm::dot(lightDirection, closestIntersection.intersectedTriangle.normal));
        pL = pL * incidenceAngle;

        // Specular Lighting
        float glossy = 256;
        glm::vec3 reflectionVector = lightDirection - (2.0f * closestIntersection.intersectedTriangle.normal * glm::dot(closestIntersection.intersectedTriangle.normal,lightDirection));
        float sL = glm::pow(glm::dot(reflectionVector, view), glossy);

        // Calculate pixel brightness
        // Add all type of light
        brightness = pL + sL + ambient;
    }

    // Using the brightness to multiply each RGB channel
    // The colour can not greater than 255
    float red = std::min((float(colour.red) * brightness), 255.0f);
    float green = std::min((float(colour.green) * brightness), 255.0f);
    float blue = std::min((float(colour.blue) * brightness), 255.0f);
    if(shadowType==HARD) {
        if(closestIntersection.triangleIndex != lightIntersection.triangleIndex && closestIntersection.intersectedTriangle.colour.name != "Yellow" && closestIntersection.intersectedTriangle.colour.name != "Orange") {
            if(lightIntersection.intersectedTriangle.colour.name == "Red") {
                ambient = 0.7;
            }

            red = red * ambient;
            green = green * ambient;
            blue = blue * ambient;


        }
    } else if (shadowType==SOFT){
        // Soft shadow
        if (closestIntersection.triangleIndex != lightIntersection.triangleIndex && closestIntersection.intersectedTriangle.colour.name != "Yellow" && closestIntersection.intersectedTriangle.colour.name != "Orange") {
            float number = 0.0f;
            for(glm::vec3 i : lightPositions) {
                glm::vec3 direction = glm::normalize(closestIntersection.intersectionPoint - i);
                RayTriangleIntersection intersection = getClosestIntersection(i, direction, modelTriangles);
                if(intersection.triangleIndex == closestIntersection.triangleIndex || intersection.intersectedTriangle.colour.name == "Red") {
                    number = number + 1.0f;
                }
            }
            auto shadowValue = glm::clamp<float>((number/float(lightPositions.size()) + ambient), 0, 1);
//            auto shadowValue = glm::clamp<float>(1.3f/ambient *(number/float(lightPositions.size()) ), 0, 1);

            if(lightIntersection.intersectedTriangle.colour.name == "Red") {
                shadowValue = 0.8;
            }

            red = red * shadowValue;
            green = green  * shadowValue;
            blue = blue * shadowValue;

        }
    } else {
        return {int(red), int(green), int(blue)};
    }
    return {int(red), int(green), int(blue)};
}

void rayTrace(DrawingWindow &window, const std::vector<ModelTriangle>& modelTriangles, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, glm::vec3 lightPosition, float focalLength, float scalingFactor, float lightPower, float ambient, const TextureMap& textureMap, shadow shadowType, shading shadeType) {
    // std::vector<glm::vec3> lightPositions = multiLightPosition(lightPosition);
    std::vector<glm::vec3> lightPositions = multiLightPosition(lightPosition);
    for(int y = 0; y < HEIGHT; y++) {
        for(int x = 0; x < WIDTH; x++) {
            // Convert pixel to 3D coordinate
            glm::vec3 threeDPoint = get3DPoint(CanvasPoint(float(x), float(y)), cameraPosition, cameraOrientation, focalLength, scalingFactor);
            // Shoot a ray from cameraPosition to pixel world coordinate
            glm::vec3 rayDirection = glm::normalize(threeDPoint - cameraPosition);
            rayDirection = rayDirection * glm::inverse(cameraOrientation);
            // Get the closest intersect model with ray
            RayTriangleIntersection closestIntersection = getClosestIntersection(cameraPosition, rayDirection, modelTriangles);
            // Set the point on the outside of box is black
            if (closestIntersection.intersectionPoint == glm::vec3(0, 0, 0)) {
                uint32_t c = (255 << 24) + (0 << 16) + (0 << 8) + (0);
                window.setPixelColour(x, y, c);
                continue;
            }
            Colour colour = castRay(cameraPosition, lightPosition, lightPositions, rayDirection, modelTriangles, lightPower, ambient, textureMap, 1, shadowType, shadeType);
            int red = colour.red;
            int blue = colour.blue;
            int green = colour.green;
            uint32_t c = (255 << 24) + (red << 16) + (green << 8) + (blue);
            window.setPixelColour(x, y, c);
        }
    }
}

void drawSphereWithFlatShading(DrawingWindow &window, const std::vector<ModelTriangle>& sphereModel, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, glm::vec3 lightPosition, float focalLength, float scalingFactor, float lightPower, float ambient) {
    std::cout << "Light: " << lightPosition[0] << " " << lightPosition[1] << " " << lightPosition[2] << " " << std::endl;
    std::cout << "Camera: " << cameraPosition[0] << " " << cameraPosition[1] << " " << cameraPosition[2] << " " << std::endl;
    for(int y = 0; y < HEIGHT; y++) {
        for(int x = 0; x < WIDTH; x++) {
            // Convert pixel to 3D coordinate
            CanvasPoint point = CanvasPoint(float(x), float(y));
            glm::vec3 threeDPoint = get3DPoint(point, cameraPosition, cameraOrientation, focalLength, scalingFactor);
            // Calculate the ray direction
            glm::vec3 rayDirection = glm::normalize(threeDPoint - cameraPosition);
            rayDirection = rayDirection * glm::inverse(cameraOrientation);
            // Get the closest intersect model with ray
            RayTriangleIntersection closestIntersection = getClosestIntersection(cameraPosition, rayDirection, sphereModel);
            // Get the light direction
            glm::vec3 lightDirection = glm::normalize(closestIntersection.intersectionPoint - lightPosition);
            // Get the closest intersect model with light
            RayTriangleIntersection lightIntersection = getClosestIntersection(lightPosition, lightDirection, sphereModel);

            // Proximity Lighting
            glm::vec3 pointToLight = lightPosition - closestIntersection.intersectionPoint;
            float distance = glm::length(pointToLight);
            float pL = lightPower / (4 * PI * distance * distance);
            // Angle of Incidence Lighting
            // If the angle is less than 0, the PL will be zero, which means there do not have proximity light in that point
            auto incidenceAngle = glm::clamp<float>(glm::dot(-lightDirection, closestIntersection.intersectedTriangle.normal), 0.0, 1.0);
            pL = pL * incidenceAngle;

            // Specular Lighting
            float glossy = 256;
            glm::vec3 view = glm::normalize(closestIntersection.intersectionPoint - cameraPosition);
            glm::vec3 reflectionVector = lightDirection - (2.0f * closestIntersection.intersectedTriangle.normal * glm::dot(closestIntersection.intersectedTriangle.normal, lightDirection));
            float sL = std::fabs(glm::pow(glm::dot(reflectionVector, view), glossy));

            // Calculate pixel brightness
            // Add all type of light
            point.brightness = pL + sL + ambient;

            // Draw colour
            Colour colour = closestIntersection.intersectedTriangle.colour;
            // Using the brightness to multiply each RGB channel
            // The colour can not greater than 255
            float red = std::min((float(colour.red) * point.brightness), 255.0f);
            float green = std::min((float(colour.green) * point.brightness), 255.0f);
            float blue = std::min((float(colour.blue) * point.brightness), 255.0f);

            // Set the point not in the ball is black
            if(closestIntersection.intersectionPoint == glm::vec3(0, 0, 0)) {
                red = 0;
                green = 0;
                blue = 0;
            } else {
                std::cout << closestIntersection.intersectedTriangle.normal[0] << " " << closestIntersection.intersectedTriangle.normal[1] << " " << closestIntersection.intersectedTriangle.normal[2] << std::endl;
            }

            uint32_t c = (255 << 24) + (int(round(red)) << 16) + (int(round(green)) << 8) + (int(round(blue)));
            window.setPixelColour(x, y, c);
        }
    }
}

void drawSphereWithGourandShading(DrawingWindow &window, const std::vector<ModelTriangle>& sphereModel, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, glm::vec3 lightPosition, float focalLength, float scalingFactor, float lightPower, float ambient) {
    std::cout << "Light: " << lightPosition[0] << " " << lightPosition[1] << " " << lightPosition[2] << " " << std::endl;
    std::cout << "Camera: " << cameraPosition[0] << " " << cameraPosition[1] << " " << cameraPosition[2] << " " << std::endl;
    for(int y = 0; y < HEIGHT; y++) {
        for(int x = 0; x < WIDTH; x++) {
            // Convert pixel to 3D coordinate
            CanvasPoint canvasPoint = CanvasPoint(float(x), float(y));
            glm::vec3 threeDPoint = get3DPoint(canvasPoint, cameraPosition, cameraOrientation, focalLength, scalingFactor);
            // Calculate the ray direction
            glm::vec3 rayDirection = glm::normalize(threeDPoint - cameraPosition);
            rayDirection = rayDirection * glm::inverse(cameraOrientation);
            // Get the closest intersect model with ray
            RayTriangleIntersection closestIntersection = getClosestIntersection(cameraPosition, rayDirection, sphereModel);
            // Set the point not in the ball is black
            // If the ray not intersect with model, the intersection point will be (0, 0, 0)
            if(closestIntersection.intersectionPoint == glm::vec3(0, 0, 0)) {
                uint32_t c = (255 << 24) + (0 << 16) + (0 << 8) + (0);
                window.setPixelColour(x, y, c);
                continue;
            }
            // Get the light direction
            glm::vec3 lightDirection = glm::normalize(closestIntersection.intersectionPoint - lightPosition);
            // Get the light intersection with each vertex
            glm::vec3 v0 = closestIntersection.intersectedTriangle.vertices[0];
            glm::vec3 v1 = closestIntersection.intersectedTriangle.vertices[1];
            glm::vec3 v2 = closestIntersection.intersectedTriangle.vertices[2];
            glm::vec3 lightDirection0 = glm::normalize(v0 - lightPosition);
            glm::vec3 lightDirection1 = glm::normalize(v1 - lightPosition);
            glm::vec3 lightDirection2 = glm::normalize(v2 - lightPosition);

            // Get the closest intersect model with light
            RayTriangleIntersection lightIntersection = getClosestIntersection(lightPosition, lightDirection, sphereModel);

            glm::vec3 point = closestIntersection.intersectionPoint;
            // Calculate the vertex normal
            glm::vec3 vn0 = vertexNormalCalculator(v0, sphereModel);
            glm::vec3 vn1 = vertexNormalCalculator(v1, sphereModel);
            glm::vec3 vn2 = vertexNormalCalculator(v2, sphereModel);

            // Proximity Lighting
            glm::vec3 v0ToLight = lightPosition - v0;
            glm::vec3 v1ToLight = lightPosition - v1;
            glm::vec3 v2ToLight = lightPosition - v2;
            float distance0 = glm::length(v0ToLight);
            float distance1 = glm::length(v1ToLight);
            float distance2 = glm::length(v2ToLight);
            float pL0 = lightPower / (4 * PI * distance0 * distance0);
            float pL1 = lightPower / (4 * PI * distance1 * distance1);
            float pL2 = lightPower / (4 * PI * distance2 * distance2);
            // Angle of Incidence Lighting
            float incidenceAngle0 = std::max(0.0f, glm::dot(-lightDirection0, vn0));
            float incidenceAngle1 = std::max(0.0f, glm::dot(-lightDirection1, vn1));
            float incidenceAngle2 = std::max(0.0f, glm::dot(-lightDirection2, vn2));
            pL0 = pL0 * incidenceAngle0;
            pL1 = pL1 * incidenceAngle1;
            pL2 = pL2 * incidenceAngle2;

            // Specular Lighting
            float glossy = 256;
            glm::vec3 view0 = glm::normalize(v0 - cameraPosition);
            glm::vec3 view1 = glm::normalize(v1 - cameraPosition);
            glm::vec3 view2 = glm::normalize(v2 - cameraPosition);
            glm::vec3 reflectionVector0 = lightDirection0 - (2.0f * vn0 * glm::dot(vn0, lightDirection0));
            glm::vec3 reflectionVector1 = lightDirection1 - (2.0f * vn1 * glm::dot(vn1, lightDirection1));
            glm::vec3 reflectionVector2 = lightDirection2 - (2.0f * vn2 * glm::dot(vn2, lightDirection2));
            float sL0 = glm::pow(glm::dot(reflectionVector0, view0), glossy);
            float sL1 = glm::pow(glm::dot(reflectionVector1, view1), glossy);
            float sL2 = glm::pow(glm::dot(reflectionVector2, view2), glossy);

            // Calculate pixel brightness
            // Add all type of light
            float brightness0 = pL0 + sL0 + ambient;
            float brightness1 = pL1 + sL1 + ambient;
            float brightness2 = pL2 + sL2 + ambient;


            // Calculate Barycentric coordinates
            float alpha = (-(point.x-v1.x)*(v2.y-v1.y)+(point.y-v1.y)*(v2.x-v1.x))/(-(v0.x-v1.x)*(v2.y-v1.y)+(v0.y-v1.y)*(v2.x-v1.x));
            float beta = (-(point.x-v2.x)*(v0.y-v2.y)+(point.y-v2.y)*(v0.x-v2.x))/(-(v1.x-v2.x)*(v0.y-v2.y)+(v1.y-v2.y)*(v0.x-v2.x));
            float gamma = 1 - alpha - beta;
            // Interpolation
            canvasPoint.brightness = alpha*brightness0 + beta*brightness1 + gamma*brightness2;

            // Draw colour
            Colour colour = closestIntersection.intersectedTriangle.colour;
            float red = std::min((float(colour.red) * canvasPoint.brightness), 255.0f);
            float green = std::min((float(colour.green) * canvasPoint.brightness), 255.0f);
            float blue = std::min((float(colour.blue) * canvasPoint.brightness), 255.0f);

            uint32_t c = (255 << 24) + (int(round(red)) << 16) + (int(round(green)) << 8) + (int(round(blue)));
            window.setPixelColour(x, y, c);
        }
    }
}

void drawSphereWithPhoneShading(DrawingWindow &window, const std::vector<ModelTriangle>& sphereModel, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, glm::vec3 lightPosition, float focalLength, float scalingFactor, float lightPower, float ambient) {
    std::cout << "Light: " << lightPosition[0] << " " << lightPosition[1] << " " << lightPosition[2] << " " << std::endl;
    std::cout << "Camera: " << cameraPosition[0] << " " << cameraPosition[1] << " " << cameraPosition[2] << " " << std::endl;
    for(int y = 0; y < HEIGHT; y++) {
        for(int x = 0; x < WIDTH; x++) {
            // Convert pixel to 3D coordinate
            CanvasPoint canvasPoint = CanvasPoint(float(x), float(y));
            glm::vec3 threeDPoint = get3DPoint(canvasPoint, cameraPosition, cameraOrientation, focalLength, scalingFactor);
            // Calculate the ray direction
            glm::vec3 rayDirection = glm::normalize(threeDPoint - cameraPosition);
            rayDirection = rayDirection * glm::inverse(cameraOrientation);
            // Get the closest intersect model with ray
            RayTriangleIntersection closestIntersection = getClosestIntersection(cameraPosition, rayDirection, sphereModel);
            // Get the light direction
            glm::vec3 lightDirection = glm::normalize(closestIntersection.intersectionPoint - lightPosition);
            // Get the closest intersect model with light
            RayTriangleIntersection lightIntersection = getClosestIntersection(lightPosition, lightDirection, sphereModel);

            glm::vec3 point = closestIntersection.intersectionPoint;
            // Calculate the vertex normal
            glm::vec3 v0 = closestIntersection.intersectedTriangle.vertices[0];
            glm::vec3 v1 = closestIntersection.intersectedTriangle.vertices[1];
            glm::vec3 v2 = closestIntersection.intersectedTriangle.vertices[2];
            glm::vec3 vn0 = vertexNormalCalculator(v0, sphereModel);
            glm::vec3 vn1 = vertexNormalCalculator(v1, sphereModel);
            glm::vec3 vn2 = vertexNormalCalculator(v2, sphereModel);

            // Calculate Barycentric coordinates
            float alpha = (-(point.x-v1.x)*(v2.y-v1.y)+(point.y-v1.y)*(v2.x-v1.x))/(-(v0.x-v1.x)*(v2.y-v1.y)+(v0.y-v1.y)*(v2.x-v1.x));
            float beta = (-(point.x-v2.x)*(v0.y-v2.y)+(point.y-v2.y)*(v0.x-v2.x))/(-(v1.x-v2.x)*(v0.y-v2.y)+(v1.y-v2.y)*(v0.x-v2.x));
            float gamma = 1 - alpha - beta;

            // Interpolation for normal vector
            glm::vec3 normal = alpha*vn0 + beta*vn1 + gamma*vn2;

            // Proximity Lighting
            glm::vec3 pointToLight = lightPosition - closestIntersection.intersectionPoint;
            float distance = glm::length(pointToLight);
            float pL = lightPower / (4 * PI * distance * distance);
            // Angle of Incidence Lighting
            float incidenceAngle = std::max(0.0f, glm::dot(lightDirection, normal));
            pL = pL * incidenceAngle;

            // Specular Lighting
            float glossy = 256;
            glm::vec3 view = glm::normalize(closestIntersection.intersectionPoint - cameraPosition);
            glm::vec3 reflectionVector = lightDirection - (2.0f * normal * glm::dot(normal, lightDirection));
            float sL = glm::pow(glm::dot(reflectionVector, view), glossy);

            // Calculate pixel brightness
            // Add all type of light
            canvasPoint.brightness = pL + sL + ambient;

            // Draw colour
            Colour colour = closestIntersection.intersectedTriangle.colour;
            float red = std::min((float(colour.red) * canvasPoint.brightness), 255.0f);
            float green = std::min((float(colour.green) * canvasPoint.brightness), 255.0f);
            float blue = std::min((float(colour.blue) * canvasPoint.brightness), 255.0f);

            // Set the point not in the ball is black
            if(closestIntersection.intersectionPoint == glm::vec3(0, 0, 0)) {
                red = 0;
                green = 0;
                blue = 0;
            }
            uint32_t c = (255 << 24) + (int(round(red)) << 16) + (int(round(green)) << 8) + (int(round(blue)));
            window.setPixelColour(x, y, c);
        }
    }
}

// ---------------------- Setting ---------------------- //
glm::vec3 cameraPosition (0.0, 0.0, 4.0);
glm::vec3 lightPosition (0.0, 0.5, 0.0);
glm::mat3 cameraOrientation (1, 0, 0, 0, 1, 0, 0, 0, 1);
float lightPower = 4.0f;
float ambient = 0.4;
float focalLength = 2.0;
bool orbitYes = false;
bool rasterisedRenderYes = false;
int scene = 0;
shadow shadowType = HARD;
shading shadeType = FLAT;

// ---------------------- Box Model ---------------------- //
std::vector<std::string> objFileName {"textured-cornell-box.obj", "sphere.obj", "logo.obj"};
std::vector<std::string> mtlFileName {"textured-cornell-box.mtl", "null", "null"};
std::vector<ModelTriangle> modelTriangles = objReader2(objFileName, mtlFileName, 0.35);
TextureMap textureMap = getTextureMap("chessboard.ppm");

// ---------------------- Sphere Model ---------------------- //
std::vector<ModelTriangle> sphereModel = objReader("sphere.obj", "null", 0.35);

// ---------------------- Show in window ---------------------- //
void handleEvent(SDL_Event event, DrawingWindow &window) {
	if(event.type == SDL_KEYDOWN) {
		if(event.key.keysym.sym == SDLK_LEFT) {
            std::cout << "LEFT" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_RIGHT) {
            std::cout << "RIGHT" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_UP) {
            std::cout << "UP" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_DOWN) {
            std::cout << "DOWN" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_f) {
            std::cout << "FORWARD" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_b) {
            std::cout << "BACKWARD" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_d) {
            std::cout << "Rotate about X-axis Clockwise" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
            cameraOrientation = rotateCameraOrientation(event, cameraOrientation);
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_a) {
            std::cout << "Rotate about X-axis Anticlockwise" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
            cameraOrientation = rotateCameraOrientation(event, cameraOrientation);
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_w) {
            std::cout << "Rotate about Y-axis Clockwise" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
            cameraOrientation = rotateCameraOrientation(event, cameraOrientation);
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_s) {
            std::cout << "Rotate about Y-axis Anticlockwise" << std::endl;
            window.clearPixels();
            clearDepthBuffer();
            cameraPosition = moveCamera(event, cameraPosition);
            cameraOrientation = rotateCameraOrientation(event, cameraOrientation);
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_1) {
            clearDepthBuffer();
            window.clearPixels();
            scene = 1;
            shadowType = NONE;
            shadeType = FLAT;
            rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                     float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
        } else if(event.key.keysym.sym == SDLK_2) {
            clearDepthBuffer();
            window.clearPixels();
            scene = 2;
            shadowType = SOFT;
            rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                     float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
        } else if(event.key.keysym.sym == SDLK_3) {
            clearDepthBuffer();
            window.clearPixels();
            scene = 3;
            shadowType = HARD;
            shadeType = FLAT;
            rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                     float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
        } else if(event.key.keysym.sym == SDLK_4) {
            clearDepthBuffer();
            window.clearPixels();
            cameraPosition = glm::vec3(0.0, 0.7, 5.2);
            lightPosition = glm::vec3(-0.3, 1.0, 0.8);
            cameraOrientation = glm::mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
            lightPower = 4.0f;
            ambient = 0.0;
            focalLength = 2.0;
            scene = 4;
            drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
        } else if(event.key.keysym.sym == SDLK_5) {
            clearDepthBuffer();
            window.clearPixels();
            cameraPosition = glm::vec3(0.0, 0.7, 5.2);
            lightPosition = glm::vec3(-0.3, 1.0, 0.8);
            cameraOrientation = glm::mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
            lightPower = 4.0f;
            ambient = 0.4;
            focalLength = 2.0;
            scene = 5;
            drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
        } else if(event.key.keysym.sym == SDLK_6) {
            clearDepthBuffer();
            window.clearPixels();
            cameraPosition = glm::vec3(0.0, 0.7, 5.2);
            lightPosition = glm::vec3(-0.3, 1.0, 0.8);
            cameraOrientation = glm::mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
            lightPower = 4.0f;
            ambient = 0.4;
            focalLength = 2.0;
            scene = 6;
            drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
        } else if(event.key.keysym.sym == SDLK_x) {
            clearDepthBuffer();
            window.clearPixels();
            cameraOrientation = glm::mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
            lightPower = 4.0f;
            ambient = 0.4;
            focalLength = 2.0;
            scene = 7;
            shadeType = GOURAND;
            rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                     float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
        } else if(event.key.keysym.sym == SDLK_z) {
            clearDepthBuffer();
            window.clearPixels();
            cameraOrientation = glm::mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
            lightPower = 4.0f;
            ambient = 0.4;
            focalLength = 2.0;
            scene = 8;
            shadeType = PHONG;
            rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                     float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
        } else if(event.key.keysym.sym == SDLK_i) {
            std::cout << "Light UP" << std::endl;
            window.clearPixels();
            lightPosition = moveLight(event, lightPosition);
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_k) {
            std::cout << "Light DOWN" << std::endl;
            window.clearPixels();
            lightPosition = moveLight(event, lightPosition);
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_j) {
            std::cout << "Light LEFT" << std::endl;
            window.clearPixels();
            lightPosition = moveLight(event, lightPosition);
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_l) {
            std::cout << "Light RIGHT" << std::endl;
            window.clearPixels();
            lightPosition = moveLight(event, lightPosition);
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_m) {
            std::cout << "Light FORWARD" << std::endl;
            window.clearPixels();
            lightPosition = moveLight(event, lightPosition);
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_n) {
            std::cout << "Light BACKWARD" << std::endl;
            window.clearPixels();
            lightPosition = moveLight(event, lightPosition);
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_9) {
            std::cout << "Plus power" << std::endl;
            window.clearPixels();
            lightPower = lightPower + 1.0f;
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_0) {
            std::cout << "Minus power" << std::endl;
            window.clearPixels();
            lightPower = lightPower - 1.0f;
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_7) {
            std::cout << "Plus ambient power" << std::endl;
            window.clearPixels();
            ambient = ambient + 0.1f;
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_8) {
            std::cout << "Minus ambient power" << std::endl;
            window.clearPixels();
            ambient = ambient - 0.1f;
            if (scene == 4){
                drawSphereWithFlatShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition,focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 5){
                drawSphereWithGourandShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else if (scene == 6){
                drawSphereWithPhoneShading(window, sphereModel, cameraPosition, cameraOrientation, lightPosition, focalLength, float(HEIGHT) * 2 / 3, lightPower, ambient);
            } else {
                rayTrace(window, modelTriangles, cameraPosition, cameraOrientation, lightPosition, focalLength,
                         float(HEIGHT) * 2 / 3, lightPower, ambient, textureMap, shadowType, shadeType);
            }
        } else if(event.key.keysym.sym == SDLK_u) {
            CanvasPoint point1, point2, point3;
            CanvasTriangle triangle;
            point1 = CanvasPoint(float(rand()%window.width), float(rand()%window.height));
            point2 = CanvasPoint(float(rand()%window.width), float(rand()%window.height));
            point3 = CanvasPoint(float(rand()%window.width), float(rand()%window.height));
            triangle = CanvasTriangle(point1, point2, point3);
            drawStrokedTriangle(window, triangle, Colour(rand()%256, rand()%256, rand()%256));
        } else if(event.key.keysym.sym == SDLK_t) {
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
            drawTextureTriangle(window, textureMap, CanvasTriangle(point1, point2, point3));
        } else if(event.key.keysym.sym == SDLK_p) {
            pointCloudRender(window, modelTriangles, cameraPosition, cameraOrientation, focalLength, float(HEIGHT)*2/3, Colour(255, 255, 255));
        } else if(event.key.keysym.sym == SDLK_g) {
            clearDepthBuffer();
            window.clearPixels();
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
            cameraPosition = glm::vec3(0.0, 0.0, 4.0);
            lightPosition = glm::vec3(0.0, 0.6, 0.3);
            cameraOrientation = glm::mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
            lightPower = 4.0f;
            ambient = 0.4;
            focalLength = 2.0;
            orbitYes = false;
            rasterisedRenderYes = false;
            scene = 0;
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
