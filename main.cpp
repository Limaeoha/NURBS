
#include <GLFW/glfw3.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <iostream>
#include <math.h>
#include <chrono>
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include <vector>

#define TWO_PI 6.28318530718
#define CURVE_ORDER_INCREMENT 1
#define PARAMETER_INCREMENT_INCREMENT 0.01
#define GEOMETRIC_FASHION_U_INCREMENT 0.01
#define RADIUS_RATIO 0.03
#define RADIUS_INCREMENT 0.1
#define GLM_FORCE_RADIANS

using namespace std::chrono;

GLFWwindow *window, *surfaceWindow;

int curveOrder;
double parameterIncrement;
enum ChangeMode {
	CHANGE_MODE_CURVE_ORDER,
	CHANGE_MODE_PARAMETER_INCREMENT,
} changeMode;

std::vector<std::pair<glm::vec2, float>> points;
int selectedPoint;
bool movingAPoint;

std::vector<double> knots;

bool displayGeomatricFashion;
double currentUForGeometricFashion;

//surface of revolution params
glm::vec2 surfacePoint1, surfacePoint2;
bool creatingSurface;
bool point1Created;

float rotX, rotZ;

bool performingRotation;
double prevMouseX, prevMouseY;

void renderPoint(glm::vec2 pos, double radius, int pointNum) {
	if (pointNum == selectedPoint) {
		glColor3f(0.3, 0.5, 0.3);
	} else {
		glColor3f(0.7, 0.2, 0.2);
	}

	glBegin(GL_POLYGON);
	for(double i = 0; i < TWO_PI; i += TWO_PI / 15) {
	 	glVertex3f(pos.x + cos(i) * radius * RADIUS_RATIO, pos.y + sin(i) * radius * RADIUS_RATIO, 0.0);
	}
	glEnd();
}

void renderPoints() {
	for (unsigned int i = 0; i < points.size(); i++) {
		renderPoint(points[i].first, points[i].second, i);
	}
}

void renderLineToMouse(glm::vec2 surfacePoint1) {
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);

	int w, h;
	glfwGetFramebufferSize(window, &w, &h);

	xpos = 2 * xpos / w - 1;
	ypos = 1 - 2 * ypos / h;

	glBegin(GL_LINES);
		glVertex3f(surfacePoint1.x, surfacePoint1.y, 0);
		glVertex3f(xpos, ypos, 0);
	glEnd();
}

void drawAxis() {
	glColor3f(0, 0, 1);
	glBegin(GL_LINES);
		glVertex3f(-1, 0, 0);
		glVertex3f(1, 0, 0);
	glEnd();

	glBegin(GL_LINES);
		glVertex3f(0, -1, 0);
		glVertex3f(0, 1, 0);
	glEnd();
}

float calculateDenominator(double u, double delta) {
	float *c = new float[curveOrder];

	for (int i=0; i<=curveOrder-1; i++) {
		c[i] = points[delta-i].second;
	}
	for (int r=curveOrder; r>=2; r--) {
		int i = delta;
		for (int s=0; s<=r-2; s++) {
			float omega = (u - knots[i])/(float)(knots[i+r-1] - knots[i]);
			c[s] = omega*c[s] + (1-omega) * c[s+1];
			i--;
		}
	}

	float toReturn = c[0];
	delete[] c;
	return toReturn;
}

void renderGeometricFashion() {
	double u = currentUForGeometricFashion;
	double delta = 0;
	for (int i=0; i<knots.size(); i++) {
		if (u >= knots[i] && u < knots[i+1]) {
			delta = i;
			break;
		}
	}
	glm::vec2 *c = new glm::vec2[curveOrder];
	float *d = new float[curveOrder];			//denominator

	for (int i=0; i<=curveOrder-1; i++) {
		c[i] = points[delta - i].first * points[delta-i].second;
		d[i] = points[delta - i].second;
	}
	for (int r=curveOrder; r>=2; r--) {
		glColor3f(1.0, 1.0, 0.0);
		for (int i=0; i<=r-2; i++) {
			glBegin(GL_LINES);
			glVertex3f(c[i].x / d[i], c[i].y / d[i], 0.0f);
			glVertex3f(c[i+1].x / d[i+1], c[i+1].y / d[i+1], 0.0f);
			glEnd();
		}
		int i = delta;
		for (int s=0; s<=r-2; s++) {
			float omega = (u - knots[i])/(float)(knots[i+r-1] - knots[i]);
			c[s] = omega*c[s] + (1-omega) * c[s+1];
			d[s] = omega*d[s] + (1-omega) * d[s+1];
			i--;
		}
	}

	renderPoint(c[0] / d[0], .5, -2);

	delete[] c;
}

glm::vec2 getCurvePoint(double u, double &nextKnot, unsigned int &delta) {
	glm::vec2 *c = new glm::vec2[curveOrder];

	if (u >= nextKnot) {
		delta++;
		nextKnot = knots[delta+1];
	}
	for (int i=0; i<=curveOrder-1; i++) {
		c[i] = points[delta - i].first * points[delta-i].second;
	}
	for (int r=curveOrder; r>=2; r--) {
		int i = delta;
		for (int s=0; s<=r-2; s++) {
			float omega = (u - knots[i])/(float)(knots[i+r-1] - knots[i]);
			c[s] = omega*c[s] + (1-omega) * c[s+1];
			i--;
		}
	}
	c[0] = c[0] / calculateDenominator(u, delta);

	glm::vec2 toReturn = c[0];

	delete[] c;

	return toReturn;
}

void renderCurve() {
	if (points.size() < 2) {
		return;
	}

	glColor3f(0.7, 0.2, 0.2);
	glBegin(GL_LINE_STRIP);
	unsigned int delta = curveOrder-1;
	double nextKnot = knots[curveOrder];

	for (double u=knots[delta]; u <= knots[points.size()]; u+=parameterIncrement/points.size()) {
		glm::vec2 point = getCurvePoint(u, nextKnot, delta);
		glVertex3f(point.x, point.y, 0.0f);
	}

	glEnd();
}

void renderSurfaceOfRevolutionLine() {
	glColor3f(0, 0.85, 0);

	glBegin(GL_LINES);
	glVertex3f(surfacePoint1.x, surfacePoint1.y, 0.0);
	glVertex3f(surfacePoint2.x, surfacePoint2.y, 0.0);
	glEnd();
}

void render () {
	glEnable (GL_DEPTH_TEST);
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//Functions for changing transformation matrix
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity ();
	//glRotatef (rotationFactor, 0.0f, 0.0f, 1.0f);
	glScalef (1, 1, 1);

	//Functions for changing projection matrix
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	glOrtho (-1, 1, -1, 1, -1, 1);

	drawAxis();
	renderPoints();
	renderCurve();

	if (point1Created) {
		renderPoint(surfacePoint1, 1, -2);
		renderLineToMouse(surfacePoint1);
	}

	if (displayGeomatricFashion) {
		renderGeometricFashion();
	}

	if (surfaceWindow != NULL) {
		renderSurfaceOfRevolutionLine();
	}
}

void createKnots() {
	knots.clear();
	for (int i=0; i<curveOrder; i++) {
		knots.push_back(0);
	}

	double inc = 1.0 / (points.size() - curveOrder + 1);
	double total = 0;
	for (unsigned int i=curveOrder; i<points.size(); i++) {
		total += inc;
		knots.push_back(total);
	}

	for (unsigned int i=points.size(); i<curveOrder + points.size(); i++) {
		knots.push_back(1);
	}
}

void renderSurface() {
	glEnable (GL_DEPTH_TEST);
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity ();
	glRotatef(rotX, 1, 0, 0);
	glRotatef(rotZ, 0, 0, 1);
	glScalef (0.5, 0.5, 0.5);

	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	glFrustum (-1, 1, -1, 1, -1, 1);

	std::vector<glm::vec2> curvePoints;
	unsigned int delta = curveOrder-1;
	double nextKnot = knots[curveOrder];

	glm::vec2 diff = (surfacePoint2 + surfacePoint1)/2.0f;
	for (double u=knots[delta]; u <= knots[points.size()]; u+=parameterIncrement/points.size()) {
		glm::vec2 point = getCurvePoint(u, nextKnot, delta);
		point = point - diff;
		curvePoints.push_back(point);
	}


	float vInc = 0.1f;
	for (double v = 0; v < TWO_PI; v+=vInc) {
		for (unsigned int i=0; i<curvePoints.size() - 1; i++) {
			glBegin(GL_POLYGON);
			glm::mat4 rotationMatrix = glm::rotate(glm::mat4(1.0f), (float)v, glm::vec3(surfacePoint2, 0.0f) - glm::vec3(surfacePoint1, 0.0f));
			glm::vec3 newPoint = rotationMatrix * glm::vec4(curvePoints[i], 0.0f, 1.0f);
			glColor3f((float)i/curvePoints.size(), v / TWO_PI, (float)i/curvePoints.size() * v / TWO_PI);
			glVertex3f(newPoint.x, newPoint.y, newPoint.z);


			newPoint = rotationMatrix * glm::vec4(curvePoints[i+1], 0.0f, 1.0f);
			glColor3f((float)i/curvePoints.size(), v / TWO_PI, (float)i/curvePoints.size() * v / TWO_PI);
			glVertex3f(newPoint.x, newPoint.y, newPoint.z);

			rotationMatrix = glm::rotate(glm::mat4(1.0f), (float)v + vInc, glm::vec3(surfacePoint2, 0.0f) - glm::vec3(surfacePoint1, 0.0f));
			newPoint = rotationMatrix * glm::vec4(curvePoints[i+1], 0.0f, 1.0f);
			glColor3f((float)i/curvePoints.size(), v / TWO_PI, (float)i/curvePoints.size() * v / TWO_PI);
			glVertex3f(newPoint.x, newPoint.y, newPoint.z);

			rotationMatrix = glm::rotate(glm::mat4(1.0f), (float)v + vInc, glm::vec3(surfacePoint2, 0.0f) - glm::vec3(surfacePoint1, 0.0f));
			newPoint = rotationMatrix * glm::vec4(curvePoints[i], 0.0f, 1.0f);
			glColor3f((float)i/curvePoints.size(), v / TWO_PI, (float)i/curvePoints.size() * v / TWO_PI);
			glVertex3f(newPoint.x, newPoint.y, newPoint.z);
			glEnd();
		}
	}

}

void keyboardSurface (GLFWwindow *sender, int key, int scancode, int action, int mods) {
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
		glfwSetWindowShouldClose(surfaceWindow, 1);
	}
}

void mouse_button_callback_surface(GLFWwindow* window, int button, int action, int mods)
{
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		performingRotation = true;
		glfwGetCursorPos(window, &prevMouseX, &prevMouseY);
	} else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
		performingRotation = false;
    }
}

void cursor_position_callback_surface(GLFWwindow* window, double xpos, double ypos)
{
	if (performingRotation) {
		double xpos, ypos;
		glfwGetCursorPos(window, &xpos, &ypos);

		rotZ -= (xpos - prevMouseX);
		rotX -= (ypos - prevMouseY);
		prevMouseX = xpos;
		prevMouseY = ypos;
	}
}

void createSurfaceWindow() {
	surfaceWindow = glfwCreateWindow (800, 800, "CPSC 589 - 2", NULL, NULL);

	glfwMakeContextCurrent (surfaceWindow);
	glfwSetKeyCallback (surfaceWindow, keyboardSurface);
	glfwSetMouseButtonCallback(surfaceWindow, mouse_button_callback_surface);
	glfwSetCursorPosCallback(surfaceWindow, cursor_position_callback_surface);
}

void createSurface() {
	if (surfaceWindow != NULL) {
		glfwDestroyWindow (surfaceWindow);
		surfaceWindow = NULL;
	}

	createSurfaceWindow();
	rotX = 0;
	rotZ = 0;
}

void validateVariables() {
	if (curveOrder < 2) {
		curveOrder = 2;
	} else if (curveOrder > points.size()) {
		curveOrder = points.size();
	}

	if (parameterIncrement < 0.01) {
		parameterIncrement = 0.01;
	}

	if (currentUForGeometricFashion < 0) {
		currentUForGeometricFashion = 0;
	} else if (currentUForGeometricFashion > 1) {
		currentUForGeometricFashion = 1;
	}
}

void printVariableInfo() {
	std::cout << "Currently changing ";
	if (changeMode == CHANGE_MODE_CURVE_ORDER) {
		std::cout << "curve order" << std::endl;
	} else if (changeMode == CHANGE_MODE_PARAMETER_INCREMENT) {
		std::cout << "parameter increment" << std::endl;
	}

	std::cout << "Curve order: " << curveOrder << std::endl;
	std::cout << "Parameter increment: " << parameterIncrement << std::endl;
	std::cout << "U for geometric fashion: " << currentUForGeometricFashion << std::endl;
}

void keyboard (GLFWwindow *sender, int key, int scancode, int action, int mods) {
	//set what to change
	if (key == GLFW_KEY_K && action == GLFW_PRESS) {
		changeMode = CHANGE_MODE_CURVE_ORDER;
	} else if (key == GLFW_KEY_U && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		changeMode = CHANGE_MODE_PARAMETER_INCREMENT;
	} else if (key == GLFW_KEY_S && (action == GLFW_PRESS)) {
		creatingSurface = !creatingSurface;
	} else if (key == GLFW_KEY_G && action == GLFW_PRESS) {
		displayGeomatricFashion = !displayGeomatricFashion;
	} else if (key == GLFW_KEY_LEFT && action == GLFW_PRESS) {
		currentUForGeometricFashion -= GEOMETRIC_FASHION_U_INCREMENT;
	} else if (key == GLFW_KEY_RIGHT && action == GLFW_PRESS) {
		currentUForGeometricFashion += GEOMETRIC_FASHION_U_INCREMENT;
	}

	//perform the actual change
	if (key == GLFW_KEY_KP_ADD && action == GLFW_PRESS) {
		curveOrder += (changeMode == CHANGE_MODE_CURVE_ORDER) ? CURVE_ORDER_INCREMENT : 0;
		parameterIncrement += (changeMode == CHANGE_MODE_PARAMETER_INCREMENT) ? PARAMETER_INCREMENT_INCREMENT : 0;
	} else if (key == GLFW_KEY_KP_SUBTRACT && action == GLFW_PRESS) {
		curveOrder -= (changeMode == CHANGE_MODE_CURVE_ORDER) ? CURVE_ORDER_INCREMENT : 0;
		parameterIncrement -= (changeMode == CHANGE_MODE_PARAMETER_INCREMENT) ? PARAMETER_INCREMENT_INCREMENT : 0;
	}
	validateVariables();

	if (action == GLFW_PRESS) {
		printVariableInfo();
	}
	createKnots();
}

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
	if (movingAPoint) {
		double xpos, ypos;
		glfwGetCursorPos(window, &xpos, &ypos);

		int w, h;
		glfwGetFramebufferSize(window, &w, &h);

		xpos = 2 * xpos / w - 1;
		ypos = 1 - 2 * ypos / h;

		points[selectedPoint].first.x = xpos;
		points[selectedPoint].first.y = ypos;
	}
}

void placePoint() {
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);

	int w, h;
	glfwGetFramebufferSize(window, &w, &h);

	xpos = 2 * xpos / w - 1;
	ypos = 1 - 2 * ypos / h;

	if (!creatingSurface) {
		bool pointExistsInProximity = false;
		for (unsigned int i=0; i<points.size(); i++) {
			double xDiff = points[i].first.x - xpos;
			double yDiff = points[i].first.y - ypos;
			if ( sqrtf(xDiff * xDiff + yDiff * yDiff) <= points[i].second * RADIUS_RATIO ) {
				pointExistsInProximity = true;
			}
		}

		if (!pointExistsInProximity) {
			points.push_back(std::pair<glm::vec2, float>(glm::vec2(xpos, ypos), 1));
			createKnots();
		}
	} else {
		if (!point1Created) {
			surfacePoint1 = glm::vec2(xpos, ypos);
			point1Created = true;
		} else {
			surfacePoint2 = glm::vec2(xpos, ypos);
			point1Created = false;
			creatingSurface = false;
			createSurface();
		}
	}

}

void removePoint() {
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);

	int w, h;
	glfwGetFramebufferSize(window, &w, &h);

	xpos = 2 * xpos / w - 1;
	ypos = 1 - 2 * ypos / h;

	int toRemove = -1;
	for (unsigned int i=0; i<points.size(); i++) {
		double xDiff = points[i].first.x - xpos;
		double yDiff = points[i].first.y - ypos;
		if ( sqrtf(xDiff * xDiff + yDiff * yDiff) <= points[i].second * RADIUS_RATIO ) {
			toRemove = i;
		}
	}

	if (toRemove != -1) {
		points.erase(points.begin()+toRemove);
		validateVariables();
		createKnots();
	}
}

void selectPoint() {
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);

	int w, h;
	glfwGetFramebufferSize(window, &w, &h);

	xpos = 2 * xpos / w - 1;
	ypos = 1 - 2 * ypos / h;

	bool pointFound = false;
	for (unsigned int i=0; i<points.size(); i++) {
		double xDiff = points[i].first.x - xpos;
		double yDiff = points[i].first.y - ypos;
		if ( sqrtf(xDiff * xDiff + yDiff * yDiff) <= points[i].second * RADIUS_RATIO ) {
			selectedPoint = i;
			pointFound = true;
		}
	}

	if (!pointFound) {
		selectedPoint = -1;
	}
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		selectPoint();
		if (selectedPoint != -1) {
			movingAPoint = true;
		}
	} else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
        placePoint();
        movingAPoint = false;
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE) {
    	removePoint();
    }
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	if (selectedPoint >= 0) {
		points[selectedPoint].second += RADIUS_INCREMENT * yoffset;
	}

	if (points[selectedPoint].second < RADIUS_INCREMENT * 4) {
		points[selectedPoint].second = RADIUS_INCREMENT * 4;
	}
}

void setup() {
	curveOrder = 2;
	parameterIncrement = 0.01;
	selectedPoint = -1;
	movingAPoint = false;
	creatingSurface = false;
	point1Created = false;
	rotX = 0;
	rotZ = 0;
	currentUForGeometricFashion = 0.5;
}

int main () {
	if (!glfwInit())
		return 1;

	window = glfwCreateWindow (800, 800, "CPSC 589 - 2", NULL, NULL);
	if (!window)
		return 1;

	setup();

	glfwSetKeyCallback (window, keyboard);
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	glfwSetCursorPosCallback(window, cursor_position_callback);
	glfwSetScrollCallback(window, scroll_callback);
	int w, h;
	while (!glfwWindowShouldClose (window)) {
		glfwMakeContextCurrent (window);
		glfwGetFramebufferSize (window, &w, &h);
		glViewport (0, 0, w, h);

		render ();

		glfwSwapBuffers (window);
		glfwPollEvents();

		if (surfaceWindow != NULL) {
			glfwMakeContextCurrent (surfaceWindow);
			glfwGetFramebufferSize (surfaceWindow, &w, &h);
			glViewport (0, 0, w, h);

			renderSurface ();

			glfwSwapBuffers (surfaceWindow);
			glfwPollEvents();

			if (glfwWindowShouldClose(surfaceWindow)) {
				glfwDestroyWindow (surfaceWindow);
				surfaceWindow = NULL;
			}
		}
	}

	if (surfaceWindow != NULL) {
		glfwDestroyWindow (surfaceWindow);
	}
	glfwDestroyWindow (window);
	glfwTerminate();
	return 0;
}

