#include "Render.h"
#include <sstream>
#include <iostream>
#include <windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>
#include "MyOGL.h"
#include "Camera.h"
#include "Light.h"
#include "Primitives.h"
#include "GUItextRectangle.h"

//===============================================================================
int lineSteps = 100; // Количество шагов для линий
int animSpeed = 10; // Скорость анимации: [1; INT_MAX] (меньше -> быстрее)
double figureSize = 0.1; // Размер фигуры (кубика)
int numberOfLines[] = { 1, 1, 2 }; // Количество линий: Безье 2-го порядка, Безье 3-го порядка, Эрмита

int N = 4; // Ширина плоскости (N * max_X)
int M = 4; // Длина плоскости (M * max_Y)
int surfaceSteps = 10; // Количество шагов для плоскостей
int numberOfSurfaces = 1; // Количество поверхностей Безье
bool drawMainPointsAndLines = true; // Рисовать опорные точки и линии?
bool drawSurfaceWithLines = true; // Рисовать поверхность с линиями?
bool drawWithQuads = true; // Рисовать поверхность квадратами? (перекрывает drawMainPointsAndLines и drawSurfaceWithLines)

// --- Do not touch params ---
int timer = 0;
double t_max = 1.0001;
double t_step = 0.0;
double surfaceStep = 0.0;

unsigned int textureId = 0;
bool textureMode = true;
bool lightMode = true;
bool isGrowing = true;

double start_X = 0.0;
double start_Y = 0.0;
double start_Z = 0.0;
double max_X = 5.0;
double max_Y = 5.0;
double max_Z = 5.0;

int lineNumber = 0;
double colors[3][3];
double *dragPoint;
std::vector <std::vector <std::vector <double>>> pointsForSurfaces;
std::vector <std::vector <std::vector <double>>> surfacesPoints;

std::vector <std::vector <std::vector <double>>> pointsForLines;
std::vector <std::vector <std::vector <double>>> linePoints;
std::vector <std::vector <std::vector <double>>> translationPoints;
std::vector <std::vector <std::vector <double>>> rotationPoints;

void getTangentVector(double *vector, std::vector <double> line_1, std::vector <double> line_2);
// --- Do not touch params ---

void loadTexture(std::string textureName) {
	RGBTRIPLE *textureArray;
	char *textureCharArray;
	int textureWidth, textureHeight;

	OpenGL::LoadBMP(textureName.c_str(), &textureWidth, &textureHeight, &textureArray);
	OpenGL::RGBtoChar(textureArray, textureWidth, textureHeight, &textureCharArray);

	glGenTextures(1, &textureId);
	glBindTexture(GL_TEXTURE_2D, textureId);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, textureWidth, textureHeight, 0,
				 GL_RGBA, GL_UNSIGNED_BYTE, textureCharArray);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	free(textureCharArray);
	free(textureArray);
}

void calculateNormal(double *firstPoint, double *secondPoint, double *generalPoint) {
	double x_1 = firstPoint[0];
	double y_1 = firstPoint[1];
	double z_1 = firstPoint[2];

	double x_2 = secondPoint[0];
	double y_2 = secondPoint[1];
	double z_2 = secondPoint[2];

	double x_3 = generalPoint[0];
	double y_3 = generalPoint[1];
	double z_3 = generalPoint[2];

	double line_1[] = { x_3 - x_1, y_3 - y_1, z_3 - z_1 };
	double line_2[] = { x_3 - x_2, y_3 - y_2, z_3 - z_2 };

	double n_x = line_1[1] * line_2[2] - line_2[1] * line_1[2];
	double n_y = line_1[0] * line_2[2] - line_2[0] * line_1[2];
	double n_z = line_1[0] * line_2[1] - line_2[0] * line_1[1];

	double length = sqrt(pow(n_x, 2) + pow(n_y, 2) + pow(n_z, 2));
	double result[] = { n_x / length, -n_y / length, n_z / length };

	glNormal3dv(result);
}

// Кривая Безье 2-го порядка
void curveBezierSecond(std::vector <std::vector <double>> points, std::vector <double> &point, double t) {
	for (int i = 0; i < 3; i++) {
		point.push_back(
			points[0][i] * pow(1.0 - t, 2.0) * 1.0 * pow(t, 0.0) +
			points[1][i] * pow(1.0 - t, 1.0) * 2.0 * pow(t, 1.0) +
			points[2][i] * pow(1.0 - t, 0.0) * 1.0 * pow(t, 2.0));
	}
}

// Кривая Безье 3-го порядка
void curveBezierThird(std::vector <std::vector <double>> points, std::vector <double> &point, double t) {
	for (int i = 0; i < 3; i++) {
		point.push_back(
			points[0][i] * pow(1.0 - t, 3.0) + 0.0 * pow(t, 0.0) +
			points[1][i] * pow(1.0 - t, 2.0) * 3.0 * pow(t, 1.0) +
			points[2][i] * pow(1.0 - t, 1.0) * 3.0 * pow(t, 2.0) +
			points[3][i] * pow(1.0 - t, 0.0) * 1.0 * pow(t, 3.0));
	}
}

// Кривая Эрмита
void curveHermite(std::vector <std::vector <double>> points, std::vector <double> &point, double t) {
	double tangentVector_1[3], tangentVector_2[3];
	getTangentVector(tangentVector_1, points[1], points[0]);
	getTangentVector(tangentVector_2, points[2], points[3]);

	for (int i = 0; i < 3; i++) {
		point.push_back(
			points[0][i] * (2.0 * pow(t, 3.0) - 3.0 * pow(t, 2.0) + 1.0) +
			points[3][i] * (3.0 * pow(t, 2.0) - 2.0 * pow(t, 3.0) + 0.0) +
			tangentVector_1[i] * (1.0 * pow(t, 3.0) - 2.0 * pow(t, 2.0) + t) +
			tangentVector_2[i] * (1.0 * pow(t, 3.0) - 1.0 * pow(t, 2.0) + 0.0));
	}
}


// Расчет факториала
double factorial(int n) {
	if (n == 0) {
		return 1.0;
	} else return n * factorial(n - 1);
}

// Расчет многочлена Бернштейна
double calculateBernsteinPolynomials(int i, int n, double u) {
	n -= 1;
	return (factorial(n) / (factorial(i) * factorial(n - i))) * pow(u, i) * pow(1.0 - u, (double) n - i);
}

// Поверхность Безье
void surfaceBezier(int numberOfSurface, double *point, double u, double v) {
	for (int i = 0; i < N * M; i++) {
		double bernsteinPolynomials_N = calculateBernsteinPolynomials(i / M, N, u);
		double bernsteinPolynomials_M = calculateBernsteinPolynomials(i - (i / M) * M, M, v);

		for (int j = 0; j < 3; j++) {
			point[j] += bernsteinPolynomials_N * bernsteinPolynomials_M * pointsForSurfaces[numberOfSurface][i][j];
		}
	}
}

void recalculateSurfaceBezier(int numberOfSurface) {
	surfacesPoints[numberOfSurface].clear();

	for (double u = 0; u < 1.0; u += surfaceStep) {
		for (double v = 0; v < 1.0; v += surfaceStep) {
			double point[3] = { 0, 0, 0 };
			surfaceBezier(numberOfSurface, point, u, v);
			surfacesPoints[numberOfSurface].push_back({ point[0], point[1], point[2] });
		}
	}
}


void setColorForLines(int lineIndex) {
	if (lineIndex < numberOfLines[0]) {
		glColor3dv(colors[0]);
	} else if (lineIndex >= numberOfLines[0] && lineIndex < (numberOfLines[1] + numberOfLines[0])) {
		glColor3dv(colors[1]);
	} else glColor3dv(colors[2]);
}

// Нарисовать линии с точками
void drawLine(double *point_1, double *point_2) {
	glLineWidth(3);
	glBegin(GL_LINES);
	glVertex3dv(point_1);
	glVertex3dv(point_2);
	glEnd();
	glLineWidth(1);

	glPointSize(10);
	glBegin(GL_POINTS);
	glVertex3dv(point_1);
	glVertex3dv(point_2);
	glEnd();
}

// Нарисовать все линии в векторе
void drawLines(std::vector <std::vector <std::vector <double>>> vector) {
	for (int i = 0; i < (int) vector.size(); i++) {
		setColorForLines(i);
		for (int j = 0; j < (int) vector[i].size() - 1; j++) {
			drawLine(&vector[i][j][0], &vector[i][j + 1][0]);
			if (i > numberOfLines[0] + numberOfLines[1] - 1) j++;
		}
	}
}


// Получить случайное число
double getRandomNumber(double min, int width) {
	return double(min + rand() % width);
}

// Генерация случайных линий
void generatePoints(bool isForSurface, int numberOfLines, int numberOfPoints, std::vector <std::vector <std::vector <double>>> &vector) {
	// Генерация линий
	for (int i = 0; i < numberOfLines; i++) {
		// Генерация точек
		std::vector <std::vector <double>> vertexs;
		for (int j = 0; j < numberOfPoints; j++) {
			// Генерация точки
			if (isForSurface) {
				std::vector <double> vertex = { start_X, start_Y, getRandomNumber(start_Z, (int) max_Z) };
				start_X -= max_X;
				if ((j + 1) % M == 0) {
					start_X = 0;
					start_Y -= max_Y;
				}
				vertexs.push_back(vertex);
			} else {
				std::vector <double> vertex = {
					getRandomNumber(start_X, (int) max_X),
					getRandomNumber(start_Y, (int) max_Y),
					getRandomNumber(start_Z, (int) max_Z) };
				vertexs.push_back(vertex);
			}
		}
		vector.push_back(vertexs);
		if (isForSurface) {
			start_Z += max_Z;
		} else start_X += max_X;
	}
}

// Получить вектор касательной по двум линиям
void getTangentVector(double *vector, std::vector <double> line_1, std::vector <double> line_2) {
	vector[0] = line_2[0] - line_1[0];
	vector[1] = line_2[1] - line_1[1];
	vector[2] = line_2[2] - line_1[2];
}


// Расчитать прямые
void calculateLines() {
	for (int curveNumber = 0; curveNumber < 3; curveNumber++) {
		int startIndex = 0;
		int endIndex = 0;

		switch (curveNumber) {
			case 0:
				startIndex = 0;
				endIndex = pointsForLines.size() - numberOfLines[1] - numberOfLines[2];
				break;
			case 1:
				startIndex = numberOfLines[0];
				endIndex = pointsForLines.size() - numberOfLines[2];
				break;
			case 2:
				startIndex = numberOfLines[0] + numberOfLines[1];
				endIndex = pointsForLines.size();
				break;
		}

		for (int i = startIndex; i < endIndex; i++) {
			std::vector <std::vector <double>> linePoints_local;
			for (double t = 0.0; t <= t_max; t += t_step) {
				std::vector <double> point;
				switch (curveNumber) {
					case 0:
						curveBezierSecond(pointsForLines[i], point, t);
						break;
					case 1:
						curveBezierThird(pointsForLines[i], point, t);
						break;
					case 2:
						curveHermite(pointsForLines[i], point, t);
						break;
				}
				linePoints_local.push_back(point);
			}
			linePoints.push_back(linePoints_local);

			std::vector <std::vector <double>> translationPoints_local;
			for (int j = 0; j < (int) linePoints_local.size() - 1; j++) {
				double tangentVector[3];
				getTangentVector(tangentVector, linePoints_local[j], linePoints_local[j + 1]);

				for (int k = 0; k < animSpeed; k++) {
					std::vector <double> point = {
						linePoints_local[j][0] + tangentVector[0] / animSpeed * k,
						linePoints_local[j][1] + tangentVector[1] / animSpeed * k,
						linePoints_local[j][2] + tangentVector[2] / animSpeed * k,
					};
					translationPoints_local.push_back(point);
				}
			}
			translationPoints.push_back(translationPoints_local);

			std::vector <std::vector <double>> rotationPoints_local;
			for (int j = 0; j < (int) translationPoints_local.size() - 1; j++) {
				double vector[3];
				getTangentVector(vector, translationPoints_local[j], translationPoints_local[j + 1]);
				double rotationVector_XY[3] = { vector[0], vector[1] };

				double rotationVector_XY_Length = sqrt(pow(rotationVector_XY[0], 2) + pow(rotationVector_XY[1], 2));
				double vector_Length = sqrt(pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2));

				rotationVector_XY[0] /= rotationVector_XY_Length;
				rotationVector_XY[1] /= rotationVector_XY_Length;

				double sinSign = rotationVector_XY[1] > 0 ? 1.0 : -1.0;
				double rotationAngle_X = acos(rotationVector_XY[0]) * 180.0 / M_PI * sinSign;
				double rotationAngle_YZ = acos(vector[2] / vector_Length) * 180.0 / M_PI - 90.0;

				std::vector <double> point = { rotationAngle_X, rotationAngle_YZ };
				rotationPoints_local.push_back(point);
			}
			rotationPoints_local.push_back(rotationPoints_local[rotationPoints_local.size() - 1]);
			rotationPoints.push_back(rotationPoints_local);

			lineNumber++;
		}
	}
}

// Расчитать поверхности
void calculateSurfaces() {
	for (int k = 0; k < (int) pointsForSurfaces.size(); k++) {
		std::vector <std::vector <double>> surfacePoints_local;
		for (double u = 0; u < 1.0; u += surfaceStep) {
			for (double v = 0; v < 1.0; v += surfaceStep) {
				double point[3] = { 0, 0, 0 };
				surfaceBezier(k, point, u, v);
				surfacePoints_local.push_back({ point[0], point[1], point[2] });
			}
		}
		surfacesPoints.push_back(surfacePoints_local);
	}
}


// Нарисовать фигуру
void drawFigure(int lineIndex) {
	glPushMatrix();
	glTranslated(translationPoints[lineIndex][timer][0],
				 translationPoints[lineIndex][timer][1],
				 translationPoints[lineIndex][timer][2]);

	if (isGrowing) {
		glRotated(rotationPoints[lineIndex][timer][0], 0, 0, 1);
		glRotated(rotationPoints[lineIndex][timer][1], 0, 1, 0);
	} else {
		glRotated(rotationPoints[lineIndex][timer][0] - 180.0, 0, 0, 1);
		glRotated(-rotationPoints[lineIndex][timer][1], 0, 1, 0);
	}

	// Направляющие
	figureSize *= 2.0;
	glDisable(GL_LIGHTING);
	glBegin(GL_LINES);
	glColor3d(1, 0, 0);
	glVertex3d(0, 0, 0);
	glVertex3d(figureSize, 0, 0);

	glColor3d(0, 1, 0);
	glVertex3d(0, 0, 0);
	glVertex3d(0, figureSize, 0);

	glColor3d(0, 0, 1);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 0, figureSize);
	glEnd();

	// Кубик
	figureSize /= 2.0;
	lightMode ? glEnable(GL_LIGHTING) : glDisable(GL_LIGHTING);
	glBegin(GL_QUADS);
	// Оранжевая грань
	glColor3d(1, 0.5, 0);
	glNormal3d(-1, 0, 0);
	glVertex3d(-figureSize, -figureSize, -figureSize);
	glVertex3d(-figureSize, figureSize, -figureSize);
	glVertex3d(-figureSize, figureSize, figureSize);
	glVertex3d(-figureSize, -figureSize, figureSize);

	// Красная грань
	glColor3d(1, 0, 0);
	glNormal3d(1, 0, 0);
	glVertex3d(figureSize, -figureSize, -figureSize);
	glVertex3d(figureSize, -figureSize, figureSize);
	glVertex3d(figureSize, figureSize, figureSize);
	glVertex3d(figureSize, figureSize, -figureSize);

	// Зеленая грань
	glColor3d(0, 1, 0);
	glNormal3d(0, -1, 0);
	glVertex3d(-figureSize, -figureSize, -figureSize);
	glVertex3d(-figureSize, -figureSize, figureSize);
	glVertex3d(figureSize, -figureSize, figureSize);
	glVertex3d(figureSize, -figureSize, -figureSize);

	// Синяя грань
	glColor3d(0, 0, 1);
	glNormal3d(0, 1, 0);
	glVertex3d(-figureSize, figureSize, -figureSize);
	glVertex3d(-figureSize, figureSize, figureSize);
	glVertex3d(figureSize, figureSize, figureSize);
	glVertex3d(figureSize, figureSize, -figureSize);

	// Желтая грань
	glColor3d(1, 1, 0);
	glNormal3d(0, 0, -1);
	glVertex3d(-figureSize, -figureSize, -figureSize);
	glVertex3d(figureSize, -figureSize, -figureSize);
	glVertex3d(figureSize, figureSize, -figureSize);
	glVertex3d(-figureSize, figureSize, -figureSize);

	// Белая грань
	glColor3d(1, 1, 1);
	glNormal3d(0, 0, 1);
	glVertex3d(-figureSize, -figureSize, figureSize);
	glVertex3d(figureSize, -figureSize, figureSize);
	glVertex3d(figureSize, figureSize, figureSize);
	glVertex3d(-figureSize, figureSize, figureSize);
	glEnd();
	glPopMatrix();
}


// Нарисовать линии
void drawLines() {
	for (int lineIndex = 0; lineIndex < lineNumber; lineIndex++) {
		setColorForLines(lineIndex);

		glLineWidth(1);
		glDisable(GL_LIGHTING);
		glBegin(GL_LINE_STRIP);
		for (int partIndex = 0; partIndex < (int) linePoints[lineIndex].size(); partIndex++) {
			glVertex3dv(&linePoints[lineIndex][partIndex][0]);
		}
		glEnd();

		drawLines(pointsForLines);
		drawFigure(lineIndex);
		lightMode ? glEnable(GL_LIGHTING) : glDisable(GL_LIGHTING);
	}
}


// Нарисовать опорные точки поверхности
void drawMainPointsOfSurface(int surfaceNumber) {
	glPointSize(10);
	glColor3d(1, 0, 0);
	glBegin(GL_POINTS);
	for (int i = 0; i < (int) pointsForSurfaces[surfaceNumber].size(); i++) {
		if (dragPoint != NULL && &pointsForSurfaces[surfaceNumber][i][0] == dragPoint) {
			glColor3d(0, 1, 0);
			glVertex3dv(&pointsForSurfaces[surfaceNumber][i][0]);
			glColor3d(1, 0, 0);
		}
		glVertex3dv(&pointsForSurfaces[surfaceNumber][i][0]);
	}
	glEnd();
}

// Нарисовать опорные линии поверхности
void drawMainLinesOfSurface(int surfaceNumber) {
	glLineWidth(3);
	glColor3d(0, 0, 1);
	glBegin(GL_LINES);
	for (int i = 0; i < (int) pointsForSurfaces[surfaceNumber].size() - (M + 1); i++) {
		if ((i + 1) % M != 0) {
			glVertex3dv(&pointsForSurfaces[surfaceNumber][i][0]);
			glVertex3dv(&pointsForSurfaces[surfaceNumber][i + 1][0]);

			glVertex3dv(&pointsForSurfaces[surfaceNumber][i + 1][0]);
			glVertex3dv(&pointsForSurfaces[surfaceNumber][i + M + 1][0]);

			if (i >= (int) pointsForSurfaces[surfaceNumber].size() - (2 * M + 1)) {
				glVertex3dv(&pointsForSurfaces[surfaceNumber][i + M + 1][0]);
				glVertex3dv(&pointsForSurfaces[surfaceNumber][i + M][0]);
			}

			if (i % M == 0) {
				glVertex3dv(&pointsForSurfaces[surfaceNumber][i + M][0]);
				glVertex3dv(&pointsForSurfaces[surfaceNumber][i][0]);
			}
		}
	}
	glEnd();
}

// Нарисовать точки поверхности
void drawSurfacePoints(int surfaceNumber) {
	glPointSize(5);
	glColor3d(0, 0, 0);
	glBegin(GL_POINTS);
	for (int i = 0; i < (int) surfacesPoints[surfaceNumber].size(); i++) {
		glVertex3dv(&surfacesPoints[surfaceNumber][i][0]);
	}
	glEnd();
}

// Нарисовать линии поверхности
void drawSurfaceLines(int surfaceNumber) {
	glLineWidth(1);
	glColor3d(0, 0, 0);
	glBegin(GL_LINES);
	for (int i = 0; i < (int) surfacesPoints[surfaceNumber].size() - (surfaceSteps + 2); i++) {
		if ((i + 1) % (surfaceSteps + 1) != 0) {
			glVertex3dv(&surfacesPoints[surfaceNumber][i][0]);
			glVertex3dv(&surfacesPoints[surfaceNumber][i + 1][0]);

			glVertex3dv(&surfacesPoints[surfaceNumber][i + 1][0]);
			glVertex3dv(&surfacesPoints[surfaceNumber][i + surfaceSteps + 2][0]);

			if (i >= (int) surfacesPoints[surfaceNumber].size() - (2 * surfaceSteps + 2)) {
				glVertex3dv(&surfacesPoints[surfaceNumber][i + surfaceSteps + 2][0]);
				glVertex3dv(&surfacesPoints[surfaceNumber][i + surfaceSteps + 1][0]);
			}

			if (i % (surfaceSteps + 1) == 0) {
				glVertex3dv(&surfacesPoints[surfaceNumber][i + surfaceSteps + 1][0]);
				glVertex3dv(&surfacesPoints[surfaceNumber][i][0]);
			}
		}
	}
	glEnd();
}

// Нарисовать поверхность квадратами
void drawSurfaceQuads(int surfaceNumber) {
	double texStep = 1.0 / surfaceSteps;
	double color = 0.0;
	double texCoord_X = 0.0;
	double texCoord_Y = 0.0;
	textureMode ? glEnable(GL_TEXTURE_2D) : glDisable(GL_TEXTURE_2D);
	lightMode ? glEnable(GL_LIGHTING) : glDisable(GL_LIGHTING);

	glBegin(GL_QUADS);
	for (int i = 0; i < (int) surfacesPoints[surfaceNumber].size() - (surfaceSteps + 2); i++) {
		textureMode ? glColor3d(1, 1, 1) : glColor3d(color, 0, 1.0 - color);
		if ((i + 1) % (surfaceSteps + 1) != 0) {
			calculateNormal(
				&surfacesPoints[surfaceNumber][i][0],
				&surfacesPoints[surfaceNumber][i + 1][0],
				&surfacesPoints[surfaceNumber][i + surfaceSteps + 2][0]);

			glTexCoord2d(texCoord_X, texCoord_Y);
			glVertex3dv(&surfacesPoints[surfaceNumber][i][0]);
			glTexCoord2d(texCoord_X + texStep, texCoord_Y);
			glVertex3dv(&surfacesPoints[surfaceNumber][i + 1][0]);
			glTexCoord2d(texCoord_X + texStep, texCoord_Y + texStep);
			glVertex3dv(&surfacesPoints[surfaceNumber][i + surfaceSteps + 2][0]);
			glTexCoord2d(texCoord_X, texCoord_Y + texStep);
			glVertex3dv(&surfacesPoints[surfaceNumber][i + surfaceSteps + 1][0]);
		} else {
			texCoord_X = -texStep;
			texCoord_Y += texStep;
		}
		color += 1.0 / pow(surfaceSteps, 2);
		texCoord_X += texStep;
	}
	glEnd();
}

// Нарисовать поверхности
void drawSurfaces() {
	for (int k = 0; k < (int) surfacesPoints.size(); k++) {
		if (drawMainPointsAndLines) {
			glDisable(GL_LIGHTING);
			drawMainPointsOfSurface(k);
			drawMainLinesOfSurface(k);
			lightMode ? glEnable(GL_LIGHTING) : glDisable(GL_LIGHTING);
		}
		if (drawWithQuads) {
			drawSurfaceQuads(k);
		} else {
			drawSurfacePoints(k);
			if (drawSurfaceWithLines) drawSurfaceLines(k);
		}
	}
}


void firstLaunchingFunction() {
	loadTexture("texture.bmp");
	t_step = (1.0 + pow(0.1, 10)) / lineSteps;
	surfaceStep = (1.0 - pow(0.1, 10)) / surfaceSteps;
	if (animSpeed < 1) animSpeed = 1;

	for (int i = 0; i < 3; i++) {
		colors[i][i] = 1.0;
		generatePoints(false, numberOfLines[i], i == 0 ? 3 : 4, pointsForLines);
	}

	for (int i = 0; i < numberOfSurfaces; i++) {
		start_X = 0;
		start_Y = 0;
		generatePoints(true, 1, N * M, pointsForSurfaces);
	}

	calculateLines();
	calculateSurfaces();
}

void animationController() {
	if (timer < lineSteps * animSpeed - 1 && isGrowing) {
		timer++;
	} else {
		isGrowing = false;
		timer--;
		if (timer <= 0) isGrowing = true;
	}
}


void myRender() {
	drawLines();
	drawSurfaces();
	animationController();
}
//===============================================================================

// Класс для настройки камеры
class CustomCamera : public Camera {
public:
	//дистанция камеры
	double camDist;
	//углы поворота камеры
	double fi1, fi2;

	//значния масеры по умолчанию
	CustomCamera() {
		camDist = 15;
		fi1 = 1;
		fi2 = 1;
	}

	//считает позицию камеры, исходя из углов поворота, вызывается движком
	void SetUpCamera() {
		//отвечает за поворот камеры мышкой
		lookPoint.setCoords(0, 0, 0);

		pos.setCoords(camDist * cos(fi2) * cos(fi1),
					  camDist * cos(fi2) * sin(fi1),
					  camDist * sin(fi2));

		if (cos(fi2) <= 0) {
			normal.setCoords(0, 0, -1);
		} else normal.setCoords(0, 0, 1);

		LookAt();
	}

	void CustomCamera::LookAt() {
		//функция настройки камеры
		gluLookAt(pos.X(), pos.Y(), pos.Z(), lookPoint.X(), lookPoint.Y(), lookPoint.Z(), normal.X(), normal.Y(), normal.Z());
	}
}  camera; // создаем объект камеры

// Класс для настройки света
class CustomLight : public Light {
public:
	CustomLight() {
		//начальная позиция света
		pos = Vector3(1, 1, 3);
	}

	//рисует сферу и линии под источником света, вызывается движком
	void  DrawLightGhismo() {
		glDisable(GL_LIGHTING);

		glColor3d(0.9, 0.8, 0);
		Sphere s;
		s.pos = pos;
		s.scale = s.scale * 0.08;
		s.Show();

		if (OpenGL::isKeyPressed('G')) {
			glColor3d(0, 0, 0);
			//линия от источника света до окружности
			glBegin(GL_LINES);
			glVertex3d(pos.X(), pos.Y(), pos.Z());
			glVertex3d(pos.X(), pos.Y(), 0);
			glEnd();

			//рисуем окруность
			Circle c;
			c.pos.setCoords(pos.X(), pos.Y(), 0);
			c.scale = c.scale * 1.5;
			c.Show();
		}
	}

	void SetUpLight() {
		GLfloat amb[] = { 0.2, 0.2, 0.2, 0 };
		GLfloat dif[] = { 1.0, 1.0, 1.0, 0 };
		GLfloat spec[] = { .7, .7, .7, 0 };
		GLfloat position[] = { pos.X(), pos.Y(), pos.Z(), 1. };

		// параметры источника света
		glLightfv(GL_LIGHT0, GL_POSITION, position);
		// характеристики излучаемого света
		// фоновое освещение (рассеянный свет)
		glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
		// диффузная составляющая света
		glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);
		// зеркально отражаемая составляющая света
		glLightfv(GL_LIGHT0, GL_SPECULAR, spec);

		glEnable(GL_LIGHT0);
	}


} light; // создаем источник света


// старые координаты мыши
int mouseX = 0, mouseY = 0;

void mouseEvent(OpenGL *ogl, int mX, int mY) {
	int dx = mouseX - mX;
	int dy = mouseY - mY;
	mouseX = mX;
	mouseY = mY;

	//меняем углы камеры при нажатой правой кнопке мыши
	if (OpenGL::isKeyPressed(VK_RBUTTON)) {
		camera.fi1 += 0.01 * dx;
		camera.fi2 += -0.01 * dy;
	}

	if (OpenGL::isKeyPressed(VK_LBUTTON)) {
		double matModelView[16], matProjection[16];
		int viewport[4];
		double x, y, z;
		glGetDoublev(GL_MODELVIEW_MATRIX, matModelView);
		glGetDoublev(GL_PROJECTION_MATRIX, matProjection);
		glGetIntegerv(GL_VIEWPORT, viewport);

		LPPOINT POINT = new tagPOINT();
		GetCursorPos(POINT);
		ScreenToClient(ogl->getHwnd(), POINT);
		POINT->y = ogl->getHeight() - POINT->y;

		for (int k = 0; k < (int) pointsForSurfaces.size(); k++) {
			for (int i = 0; i < N * M; i++) {
				gluProject(pointsForSurfaces[k][i][0],
						   pointsForSurfaces[k][i][1],
						   pointsForSurfaces[k][i][2],
						   matModelView, matProjection, viewport, &x, &y, &z);

				if (pow(x - POINT->x, 2) + pow(y - POINT->y, 2) < 200) {
					dragPoint = &pointsForSurfaces[k][i][0];
					dragPoint[2] += 0.02 * dy;
					recalculateSurfaceBezier(k);
					break;
				}
			}
		}
		delete POINT;
	}

	if (!OpenGL::isKeyPressed(VK_LBUTTON)) dragPoint = NULL;

	//двигаем свет по плоскости, в точку где мышь
	if (OpenGL::isKeyPressed('G') && !OpenGL::isKeyPressed(VK_LBUTTON)) {
		LPPOINT POINT = new tagPOINT();
		GetCursorPos(POINT);
		ScreenToClient(ogl->getHwnd(), POINT);
		POINT->y = ogl->getHeight() - POINT->y;

		Ray r = camera.getLookRay(POINT->x, POINT->y);

		double z = light.pos.Z();
		double k = 0, x = 0, y = 0;

		if (r.direction.Z() == 0) {
			k = 0;
		} else k = (z - r.origin.Z()) / r.direction.Z();

		x = k * r.direction.X() + r.origin.X();
		y = k * r.direction.Y() + r.origin.Y();

		light.pos = Vector3(x, y, z);
	}

	if (OpenGL::isKeyPressed('G') && OpenGL::isKeyPressed(VK_LBUTTON)) {
		light.pos = light.pos + Vector3(0, 0, 0.02 * dy);
	}
}

void mouseWheelEvent(OpenGL *ogl, int delta) {
	if (delta < 0 && camera.camDist <= 1) return;
	if (delta > 0 && camera.camDist >= 100)	return;

	camera.camDist += 0.01 * delta;
}


void keyDownEvent(OpenGL *ogl, int key) {
	if (key == 'T') textureMode = !textureMode;
	if (key == 'R') {
		camera.fi1 = 1;
		camera.fi2 = 1;
		camera.camDist = 15;

		light.pos = Vector3(1, 1, 3);
	}
	if (key == 'L') lightMode = !lightMode;
	if (key == 'F') light.pos = camera.pos;
}

void keyUpEvent(OpenGL *ogl, int key) {

}


// Выполняется перед первым рендером
void initRender(OpenGL *ogl) {
	firstLaunchingFunction();

	//камеру и свет привязываем к "движку"
	ogl->mainCamera = &camera;
	ogl->mainLight = &light;

	// нормализация нормалей : их длины будет равна 1
	glEnable(GL_NORMALIZE);

	// устранение ступенчатости для линий
	glEnable(GL_LINE_SMOOTH);

	//  задать параметры освещения
	//  параметр GL_LIGHT_MODEL_TWO_SIDE - 
	//                0 -  лицевые и изнаночные рисуются одинаково(по умолчанию), 
	//                1 - лицевые и изнаночные обрабатываются разными режимами       
	//                соответственно лицевым и изнаночным свойствам материалов.    
	//  параметр GL_LIGHT_MODEL_AMBIENT - задать фоновое освещение, 
	//                не зависящее от сточников
	//  по умолчанию (0.2, 0.2, 0.2, 1.0)

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);

	camera.fi1 = -1.3;
	camera.fi2 = 0.8;
}

void Render(OpenGL *ogl) {
	glEnable(GL_DEPTH_TEST);

	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// настройка материала
	GLfloat amb[] = { 0.2, 0.2, 0.1, 1. };
	GLfloat dif[] = { 0.4, 0.65, 0.5, 1. };
	GLfloat spec[] = { 0.9, 0.8, 0.3, 1. };
	GLfloat sh = 0.1f * 256;

	// фоновая
	glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
	// дифузная
	glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
	// зеркальная
	glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
	// размер блика
	glMaterialf(GL_FRONT, GL_SHININESS, sh);

	// чтоб было красиво, без квадратиков (сглаживание освещения)
	glShadeModel(GL_SMOOTH);
	//===================================
	// Прогать тут
	glBindTexture(GL_TEXTURE_2D, textureId);
	myRender();
	//===================================

   // Cообщение вверху экрана
	glMatrixMode(GL_PROJECTION);	// Делаем активной матрицу проекций. 
									// (всек матричные операции, будут ее видоизменять.)
	glPushMatrix();					// сохраняем текущую матрицу проецирования (которая описывает перспективную проекцию) в стек 				    
	glLoadIdentity();				// Загружаем единичную матрицу
	glOrtho(0, ogl->getWidth(), 0, ogl->getHeight(), 0, 1);	 // врубаем режим ортогональной проекции
	glMatrixMode(GL_MODELVIEW);		// переключаемся на модел-вью матрицу
	glPushMatrix();					// сохраняем текущую матрицу в стек (положение камеры, фактически)
	glLoadIdentity();				// сбрасываем ее в дефолт
	glDisable(GL_LIGHTING);


	GuiTextRectangle rec;			// классик моего авторства для удобной работы с рендером текста.
	rec.setSize(300, 400);
	rec.setPosition(10, ogl->getHeight() - 400 - 10);

	std::stringstream ss;
	std::string textureModeState = textureMode ? "[вкл]" : "[выкл]";
	std::string lightModeState = lightMode ? "[вкл]" : "[выкл]";

	ss << "T -> вкл/выкл текстур " << textureModeState << std::endl;
	ss << "L -> вкл/выкл освещение " << lightModeState << std::endl;
	ss << "F -> Свет из камеры" << std::endl;
	ss << "G -> двигать свет по горизонтали" << std::endl;
	ss << "G + ЛКМ -> двигать свет по вертекали" << std::endl;
	ss << "B -> переключить текстуру" << std::endl;
	ss << std::endl << std::endl;

	ss << "Light Coords: (" << light.pos.X() << ", " << light.pos.Y() << ", " << light.pos.Z() << ")" << std::endl;
	ss << "Cam Coords: (" << camera.pos.X() << ", " << camera.pos.Y() << ", " << camera.pos.Z() << ")" << std::endl;
	ss << "Cam Params: R = " << camera.camDist << ", fi1 = " << camera.fi1 << ", fi2 = " << camera.fi2 << std::endl;

	rec.setText(ss.str().c_str());
	rec.Draw();

	glMatrixMode(GL_PROJECTION);	  // восстанавливаем матрицы проекции и модел-вью обратьно из стека.
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}