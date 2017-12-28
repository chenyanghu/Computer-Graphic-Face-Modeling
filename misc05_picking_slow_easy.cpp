// Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <array>
#include <stack>
#include <sstream>
// Include GLEW
#include <GL/glew.h>
// Include GLFW
#include <glfw3.h>
// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
using namespace glm;
// Include AntTweakBar
#include <AntTweakBar.h>

#include <common/shader.hpp>
#include <common/controls.hpp>
#include <common/objloader.hpp>
#include <common/vboindexer.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <common/tga.h>
#include <common/ray_casting.h>
using namespace std;
#define then
#ifndef TRUE
#define TRUE  (1)
#define FALSE (0)
#endif

const int window_width = 600, window_height = 600;

typedef struct Vertex {
	float Position[4];
	float Color[4];
	float Normal[3];
	float UV[2];
	void SetPosition(float *coords) {
		Position[0] = coords[0];
		Position[1] = coords[1];
		Position[2] = coords[2];
		Position[3] = 1.0;
	}
	void SetColor(float *color) {
		Color[0] = color[0];
		Color[1] = color[1];
		Color[2] = color[2];
		Color[3] = color[3];
	}
	void SetNormal(float *coords) {
		Normal[0] = coords[0];
		Normal[1] = coords[1];
		Normal[2] = coords[2];
	}
	void SetUV(float *coords) {
		UV[0] = coords[0];
		UV[1] = coords[1];
	}
} Vertex;

// function prototypes
int initWindow(void);
void initOpenGL(void);
void SMVerts(void);
void loadObject(char*, glm::vec4, Vertex * &, GLuint* &, int);
void createVAOs(Vertex[], GLuint[], int);
void createObjects(void);
void pickObject(void);
void renderScene(void);
void cleanup(void);
static void keyCallback(GLFWwindow*, int, int, int, int);
static void mouseCallback(GLFWwindow*, int, int, int);
void fileexport();
void fileload();
float pickingColor[441];

// GLOBAL VARIABLES
GLFWwindow* window;

glm::mat4 gProjectionMatrix;
glm::mat4 gViewMatrix;

GLuint programID;
GLuint pickingProgramID;
GLuint gPickedIndex = -1;
GLuint pickingColorID;
std::string gMessage;
bool Ani_smile = false;
bool Ani_angry = false;
float time_smile = 0.0;
float time_angry = 0.0;
bool isPick = false;
bool is_zDirection = false;
glm::vec3 curCoor;
glm::vec3 oriCoor;
float zpos = 0.0;
int inputMode = 0;
// For camera rotating
double CAngleXZ = 0.0;
double CAngleXY = 0.0;
double CRadius = 0.0;
glm::vec3 cameraCoords;
float cameraUp;
int direction = 0;
// For rendering face model
bool showFace = false;
// For drawing control mesh
bool showControlMesh = false;
float meshOffset[3] = { 0.5, -0.75, -0.25 };
Vertex CtlPoints[441];
Vertex CtlMeshVerts[2400];
Vertex SmMeshVerts[19494];
GLuint texture;
// For refresh/init control mesh
void updateControlMesh(bool isInit);
void opepoints(void);
// For picking
GLuint pickedIndex = -1;
//bool holdingShift = false;
bool enableDrag = false;
double lastXPos = -1, lastYPos = -1;
bool smooth = false;
bool hidesverts = false;
const GLuint NumObjects = 7;	// ATTN: THIS NEEDS TO CHANGE AS YOU ADD NEW OBJECTS
/*
0 = XYZ axes
1 = face model
2 = XZ grid
3 = Ctl Mesh Points
4 = Ctl Mesh Lines
5 = Triangles of the green grids
*/

GLuint VertexArrayId[NumObjects] = { 0 };
GLuint VertexBufferId[NumObjects] = { 0 };
GLuint IndexBufferId[NumObjects] = { 0 };

size_t NumIndices[NumObjects] = { 0 };
size_t VertexBufferSize[NumObjects] = { 0 };
size_t IndexBufferSize[NumObjects] = { 0 };

Vertex SMeshPoints[3249];
Vertex GreenGridVerts[1680];
GLuint MatrixID;
GLuint ModelMatrixID;
GLuint ViewMatrixID;
GLuint ProjMatrixID;
GLuint PickingMatrixID;
GLuint pickingColorArrayID;
GLuint LightID;
GLuint useLightModelID;
GLuint useTextureID;
GLuint TextureID;

bool hide = false;// hide the control points and mesh

void opepoints()
{
	if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
			is_zDirection = true;
	if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_RELEASE)
			is_zDirection = false;
}

int load_TGA(tTGA *tga, const char *filename) {

#define SIZEOF_TGA_HEADER 18

	unsigned char   buffer[256];

	int             size_of_image_id;
	int             is_colormap;
	int             targa_type;
	int             colormap_origin;
	unsigned int    colormap_length;
	int             colormap_entry_size;
	int             image_pixel_size;
	int             image_descriptor;
	int             is_inverted;

	int             image_width;
	int             image_height;

	unsigned char  *colormap;
	FILE           *f;
	unsigned char  *data;
	int             x, y, i, j;
	int             raster_width;

	if ((f = fopen(filename, "rb")) == NULL) return FALSE;

	/* header info */
	if (fread(buffer, 1, SIZEOF_TGA_HEADER, f) != SIZEOF_TGA_HEADER) return FALSE;

	size_of_image_id = buffer[0];
	is_colormap = buffer[1];
	targa_type = buffer[2];

	colormap_origin = buffer[3] + ((int)(buffer[4]) << 8);
	colormap_length = buffer[5] + ((int)(buffer[6]) << 8);
	colormap_entry_size = buffer[7];

	image_width = buffer[12] + ((unsigned)(buffer[13]) << 8);
	image_height = buffer[14] + ((unsigned)(buffer[15]) << 8);
	image_pixel_size = buffer[16];
	image_descriptor = buffer[17];

	/* valid type? */
	if ((targa_type != 1) && (targa_type != 2)) return FALSE;

	/* colormap required but missing? */
	if ((targa_type == 1) && !is_colormap) return FALSE;

	/* cannot load direct-color images */
	if ((targa_type == 2) && is_colormap) return FALSE;

	/* image id */
	if (size_of_image_id)
	if ((int)fread(buffer, 1, size_of_image_id, f) != size_of_image_id)
		return FALSE;

	is_inverted = (image_descriptor & 0x10) != 0;

	/* cannot handle interlacing */
	if ((image_descriptor & 0xC0))
		return FALSE;

	/* assume that targa 32 contains alpha (image_descriptor bits 0..3) */

	/* load colormap, if any */
	if (is_colormap)
	{
		/* must be targa 24 or targa 32 */
		if ((colormap_entry_size != 24) && (colormap_entry_size != 32)) return FALSE;

		/* convert to number of bytes/color entry */
		colormap_entry_size >>= 3;

		colormap = (unsigned char*)malloc(colormap_length *colormap_entry_size);
		if (colormap == NULL) return FALSE;

		if (fread(colormap, colormap_entry_size, colormap_length, f) != colormap_length)
		{
		lerror:
			free(colormap);
			return FALSE;
		}

		/* initializations */
		image_pixel_size = (image_pixel_size + 7) >> 3;
		raster_width = image_width *colormap_entry_size;
	}
	else {
		/* must be targa 24 or targa 32 */
		if ((image_pixel_size != 24) && (image_pixel_size != 32)) return FALSE;
		image_pixel_size >>= 3;
		raster_width = image_width *image_pixel_size;
	}

	data = (unsigned char*)malloc(raster_width *image_height);
	if (data == NULL)
		goto lerror;

	/* load image data */
	for (y = (is_inverted ? (image_height - 1) : 0);
		(is_inverted ? (y >= 0) : (y < (int)image_height));
		(is_inverted ? (--y) : (++y)))
	for (x = 0; x < image_width; x++) {

		/* get the next pixel */
		if ((int)fread(buffer, 1, image_pixel_size, f) != image_pixel_size)
			goto lerror;

		/* store it */
		if (is_colormap)
		{
			/* colormapped */
			i = ((buffer[0] + ((unsigned)(buffer[1]) << 8)) - colormap_origin)
				*colormap_entry_size;
			j = (y *raster_width) + (x *colormap_entry_size);

			data[j] = colormap[i + 2];
			data[j + 1] = colormap[i + 1];
			data[j + 2] = colormap[i];

			if (colormap_entry_size > 3)
				data[j + 3] = colormap[i + 3];
		}
		else {
			/* non-colormapped */
			j = (y *raster_width) + (x *image_pixel_size);

			data[j] = buffer[2];
			data[j + 1] = buffer[1];
			data[j + 2] = buffer[0];

			if (image_pixel_size > 3)
				data[j + 3] = buffer[3];
		}
	}

	/* free the colormap if we had loaded it */
	if (is_colormap)
		free(colormap);

	/* store the result */
	tga->width = image_width;
	tga->height = image_height;
	tga->data = data;
	tga->alpha = (is_colormap ? (colormap_entry_size > 3) : (image_pixel_size > 3));

#undef SIZEOF_TGA_HEADER

	return TRUE;
}

/*--------------------------------------------------------------------------+/
free_TGA
/+--------------------------------------------------------------------------*/
void free_TGA(tTGA *tga) {

	if (tga->data)
		free(tga->data);

	tga->data = NULL;
	tga->height =
		tga->width = 0;
	tga->alpha = 0;
}

/*--------------------------------------------------------------------------+/
load_texture_TGA
/+--------------------------------------------------------------------------*/
GLuint load_texture_TGA(const char *filename, long *width, long *height, GLint wrap_s, GLint wrap_t) {

	GLuint result;
	tTGA   tga;

	glEnable(GL_TEXTURE_2D);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	glGenTextures(1, &result);

	if (!load_TGA(&tga, filename))
	{
		glDeleteTextures(1, &result);
		return 0;
	}


	glBindTexture(GL_TEXTURE_2D, result);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrap_s);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrap_t);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexImage2D(
		GL_TEXTURE_2D,
		0,
		(tga.alpha) ? (GL_RGBA) : (GL_RGB),
		tga.width,
		tga.height,
		0,
		(tga.alpha) ? (GL_RGBA) : (GL_RGB),
		GL_UNSIGNED_BYTE,
		tga.data
		);

	if (width)
		*width = tga.width;
	if (height)
		*height = tga.height;

	free_TGA(&tga);

	return result;
}
void fileexport(){
	fstream file;
	file.open("cm.p3", ios::out | ios::binary);
	for (int i = 0; i < 441; i++){
		file << CtlPoints[i].Position[0] << " " << CtlPoints[i].Position[1] << " " << CtlPoints[i].Position[2]<<" ";
	}
	file.close();
}
void fileload(){
	ifstream in("cm.p3");
	string line; string s; int i = 0; stringstream ss; int j = 0; float a, b, c;
	while (in >> s)
	{
		j++;
		ss << s;
		if (j % 3 == 1) {
			ss >> CtlPoints[i].Position[0]; ss.clear();
		}
		else if (j % 3 == 2) {
			ss >> CtlPoints[i].Position[1]; ss.clear();
		}
		else {
			ss >> CtlPoints[i].Position[2]; ss.clear(); i++;
		}
	}
	VertexBufferSize[3] = sizeof(CtlPoints);	// ATTN: this needs to be done for each hand-made object with the ObjectID (subscript)
	createVAOs(CtlPoints, NULL, 3);
	updateControlMesh(1);
}
void loadObject(char* file, glm::vec4 color, Vertex * &out_Vertices, GLuint* &out_Indices, int ObjectId)
{
	// Read our .obj file
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec2> uvs;
	std::vector<glm::vec3> normals;
	if (!loadOBJ(file, vertices, normals)) return;

	std::vector<GLushort> indices;
	std::vector<glm::vec3> indexed_vertices;
	std::vector<glm::vec2> indexed_uvs;
	std::vector<glm::vec3> indexed_normals;
	indexVBO(vertices, normals, indices, indexed_vertices, indexed_normals);

	const size_t vertCount = indexed_vertices.size();
	const size_t idxCount = indices.size();

	// populate output arrays
	out_Vertices = new Vertex[vertCount];
	for (int i = 0; i < vertCount; i++) {
		out_Vertices[i].SetPosition(&indexed_vertices[i].x);
		out_Vertices[i].SetNormal(&indexed_normals[i].x);
		out_Vertices[i].SetColor(&color[0]);
	}
	out_Indices = new GLuint[idxCount];
	for (int i = 0; i < idxCount; i++) {
		out_Indices[i] = indices[i];
	}

	// set global variables!!
	NumIndices[ObjectId] = idxCount;
	VertexBufferSize[ObjectId] = sizeof(out_Vertices[0]) * vertCount;
	IndexBufferSize[ObjectId] = sizeof(GLuint)* idxCount;
}


void createObjects(void)
{
	//-- COORDINATE AXES --//
	Vertex CoordVerts[] =
	{
		{ { 0.0, 0.0, 0.0, 1.0 }, { 1.0, 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 5.0, 0.0, 0.0, 1.0 }, { 1.0, 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 0.0, 0.0, 0.0, 1.0 }, { 0.0, 1.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 0.0, 5.0, 0.0, 1.0 }, { 0.0, 1.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 0.0, 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 0.0, 0.0, 5.0, 1.0 }, { 0.0, 0.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
	};

	VertexBufferSize[0] = sizeof(CoordVerts);	// ATTN: this needs to be done for each hand-made object with the ObjectID (subscript)
	createVAOs(CoordVerts, NULL, 0);

	//-- GRID --//

	Vertex GRID[] =
	{
		{ { -5.0, 0.0, -5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -5.0, 0.0, 5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -4.0, 0.0, 5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -4.0, 0.0, -5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -3.0, 0.0, -5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -3.0, 0.0, 5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -2.0, 0.0, 5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -2.0, 0.0, -5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -1.0, 0.0, -5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -1.0, 0.0, 5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 0.0, 0.0, 5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 0.0, 0.0, -5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 1.0, 0.0, -5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 1.0, 0.0, 5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 2.0, 0.0, 5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 2.0, 0.0, -5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 3.0, 0.0, -5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 3.0, 0.0, 5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 4.0, 0.0, 5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 4.0, 0.0, -5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 5.0, 0.0, -5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 5.0, 0.0, 5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -5.0, 0.0, -5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 5.0, 0.0, -5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 5.0, 0.0, -4.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -5.0, 0.0, -4.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -5.0, 0.0, -3.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 5.0, 0.0, -3.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 5.0, 0.0, -2.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -5.0, 0.0, -2.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -5.0, 0.0, -1.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 5.0, 0.0, -1.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 5.0, 0.0, 0.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -5.0, 0.0, 0.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -5.0, 0.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 5.0, 0.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 5.0, 0.0, 2.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -5.0, 0.0, 2.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -5.0, 0.0, 3.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 5.0, 0.0, 3.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 5.0, 0.0, 4.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -5.0, 0.0, 4.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { -5.0, 0.0, 5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 5.0, 0.0, 5.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
	};
	VertexBufferSize[2] = sizeof(GRID);	// ATTN: this needs to be done for each hand-made object with the ObjectID (subscript)
	createVAOs(GRID, NULL, 2);

	// ATTN: create your grid vertices here!

	// draw Control Mesh
	float Pos[4] = { 0.0, 0.0, -5.0, 1.0 };
	float tempColor[4] = { 0.0, 1.0, 0.0, 1.0 };
	float tempNormal[3] = { 0.0, 0.0, 1.0 };
	for (int i = 0; i < 441; i++)
	{
		Pos[0] = ((-0.5)*(i / 21 - 10) + meshOffset[0] - 0.5);
		Pos[1] = (0.5*(i % 21 - 10) + 5.0f + meshOffset[1] + 0.5);
		CtlPoints[i].SetPosition(Pos);
		CtlPoints[i].SetColor(tempColor);
		// Set colors for picking
		CtlPoints[i].SetNormal(tempNormal);
	}
	VertexBufferSize[3] = sizeof(CtlPoints);
	createVAOs(CtlPoints, NULL, 3);

	updateControlMesh(1);

	for (int i = 0; i < 3249; i++)
	{
		float posi[4] = { 1.0, 0.0, 0.0, 1.0 };
		SMeshPoints[i].SetColor(posi);
		SMeshPoints[i].SetPosition(posi);
	}


	// TGA Import Functions
	long width = 0;
	long height = 0;
	texture = load_texture_TGA("Hu_Chenyang.tga", &width, &height, GL_CLAMP, GL_CLAMP);

	//-- .OBJs --//
	// ATTN: load your models here
	Vertex* Verts;
	GLuint* Idcs;
	//loadObject("models/base.obj", glm::vec4(1.0, 0.0, 0.0, 1.0), Verts, Idcs, ObjectID);
	//createVAOs(Verts, Idcs, ObjectID);

	loadObject("HCY.obj", glm::vec4(210.0f / 255.0f, 205.0f / 255.0f, 205.0f / 255.0f, 1.0), Verts, Idcs, 1);
	createVAOs(Verts, Idcs, 1);

}

void updateControlMesh(bool isInit) {
	//horizontal mesh lines
	for (int i = 0; i < 21; i++) {
		for (int j = 0; j < 20; j++) {
			GreenGridVerts[i * 40 + j * 2] = CtlPoints[i * 21 + j];
			GreenGridVerts[i * 40 + j * 2 + 1] = CtlPoints[i * 21 + j + 1];
		}
	}
	//vertical mesh lines
	for (int i = 0; i < 20; i++) {
		for (int j = 0; j < 21; j++) {
			GreenGridVerts[840 + i * 42 + j * 2] = CtlPoints[j + (i * 21)];
			GreenGridVerts[840 + i * 42 + j * 2 + 1] = CtlPoints[j + ((i + 1) * 21)];
		}
	}
	if (isInit) {
		VertexBufferSize[4] = sizeof(GreenGridVerts);
		createVAOs(GreenGridVerts, NULL, 4);
	}
	else {
		glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[4]);
		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(GreenGridVerts), GreenGridVerts);
	}

	int count = 0;
	int vert = 0;
	int hor = 0;
	float UV[2], u, v;
	for (int i = 0; i < 20; i++) {
		for (int j = 0; j < 20; j++) {

			vert = i * 21;
			hor = j;

			//triangle 1
			CtlMeshVerts[count].SetPosition(CtlPoints[vert + hor + 21].Position);
			if (isInit) {
				u = (CtlMeshVerts[count].Position[0] - meshOffset[0] + 5) / 10;
				v = (CtlMeshVerts[count].Position[1] - meshOffset[1]) / 10;
				UV[0] = u;
				UV[1] = v;
				CtlMeshVerts[count].SetUV(UV);
			}

			Vertex vertex2 = CtlMeshVerts[count + 1];
			vertex2.SetPosition(CtlPoints[vert + hor].Position);
			if (isInit) {
				u = (vertex2.Position[0] - meshOffset[0] + 5) / 10;
				v = (vertex2.Position[1] - meshOffset[1]) / 10;
				UV[0] = u;
				UV[1] = v;
				vertex2.SetUV(UV);
			}
			CtlMeshVerts[count + 1] = vertex2;

			Vertex vertex3 = CtlMeshVerts[count + 2];
			vertex3.SetPosition(CtlPoints[vert + hor + 22].Position);
			if (isInit) {
				u = (vertex3.Position[0] - meshOffset[0] + 5) / 10;
				v = (vertex3.Position[1] - meshOffset[1]) / 10;
				UV[0] = u;
				UV[1] = v;
				vertex3.SetUV(UV);
			}
			CtlMeshVerts[count + 2] = vertex3;

			//triangle 2
			CtlMeshVerts[count + 3].SetPosition(CtlPoints[vert + hor + 1].Position);
			if (isInit) {
				u = (CtlMeshVerts[count + 3].Position[0] - meshOffset[0] + 5) / 10;
				v = (CtlMeshVerts[count + 3].Position[1] - meshOffset[1]) / 10;
				UV[0] = u;
				UV[1] = v;
				CtlMeshVerts[count + 3].SetUV(UV);
			}

			CtlMeshVerts[count + 4] = vertex3;
			CtlMeshVerts[count + 5] = vertex2;

			count += 6;
		}
	}
	if (isInit) {
		VertexBufferSize[5] = sizeof(CtlMeshVerts);
		createVAOs(CtlMeshVerts, NULL, 5);
	}
	else {
		glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[5]);
		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(CtlMeshVerts), CtlMeshVerts);
	}
}

void showSMFace(bool isInit)
{
	int count = 0;
	int vert = 0;
	int hor = 0;
	float UV[2], u, v;
	for (int i = 0; i < 56; i++) {
		for (int j = 0; j < 56; j++) {

			vert = i * 57;
			hor = j;

			//triangle 1
			SmMeshVerts[count].SetPosition(SMeshPoints[vert + hor + 57].Position);
			if (isInit) {
				u = (SmMeshVerts[count].Position[0] - meshOffset[0] + 5) / 10;
				v = (SmMeshVerts[count].Position[1] - meshOffset[1]) / 10;
				UV[0] = u;
				UV[1] = v;
				SmMeshVerts[count].SetUV(UV);
			}

			Vertex vertex2 = SmMeshVerts[count + 1];
			vertex2.SetPosition(SMeshPoints[vert + hor].Position);
			if (isInit) {
				u = (vertex2.Position[0] - meshOffset[0] + 5) / 10;
				v = (vertex2.Position[1] - meshOffset[1]) / 10;
				UV[0] = u;
				UV[1] = v;
				vertex2.SetUV(UV);
			}
			SmMeshVerts[count + 1] = vertex2;

			Vertex vertex3 = SmMeshVerts[count + 2];
			vertex3.SetPosition(SMeshPoints[vert + hor + 58].Position);
			if (isInit) {
				u = (vertex3.Position[0] - meshOffset[0] + 5) / 10;
				v = (vertex3.Position[1] - meshOffset[1]) / 10;
				UV[0] = u;
				UV[1] = v;
				vertex3.SetUV(UV);
			}
			SmMeshVerts[count + 2] = vertex3;

			//triangle 2
			SmMeshVerts[count + 3].SetPosition(SMeshPoints[vert + hor + 1].Position);
			if (isInit) {
				u = (SmMeshVerts[count + 3].Position[0] - meshOffset[0] + 5) / 10;
				v = (SmMeshVerts[count + 3].Position[1] - meshOffset[1]) / 10;
				UV[0] = u;
				UV[1] = v;
				SmMeshVerts[count + 3].SetUV(UV);
			}

			SmMeshVerts[count + 4] = vertex3;
			SmMeshVerts[count + 5] = vertex2;

			count += 6;
		}
	}
	if (isInit) {
		VertexBufferSize[7] = sizeof(SmMeshVerts);
		createVAOs(SmMeshVerts, NULL, 7);
	}
	else {
		glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[7]);
		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(SmMeshVerts), SmMeshVerts);
	}
}

// Moving Camera
void moveCamera(double delta) {
	switch (direction) {
	case GLFW_KEY_LEFT:
		CAngleXZ -= delta;
		break;
	case GLFW_KEY_RIGHT:
		CAngleXZ += delta;
		break;
	case GLFW_KEY_UP:
		CAngleXY += delta;
		break;
	case GLFW_KEY_DOWN:
		CAngleXY -= delta;
		break;
	default:
		break;
	}
	cameraUp = 1.0;
	// Horizontal limit
	if (CAngleXZ > 2 * pi<double>())
		CAngleXZ -= 2 * pi<double>();
	if (CAngleXZ < -2 * pi<double>())
		CAngleXZ += 2 * pi<double>();

	// Vertical limit
	if (CAngleXY > 2 * pi<double>())
		CAngleXY -= 2 * pi<double>();
	if (CAngleXY < -2 * pi<double>())
		CAngleXY += 2 * pi<double>();

	if (CAngleXY > pi<double>() / 2 && CAngleXY < 3 * pi<double>() / 2) cameraUp = -1.0;  // For W
	if (CAngleXY < -pi<double>() / 2 && CAngleXY > -3 * pi<double>() / 2) cameraUp = -1.0;  // For S

	cameraCoords.x = CRadius * cos(CAngleXY) * cos(CAngleXZ);
	cameraCoords.y = CRadius * sin(CAngleXY);
	cameraCoords.z = CRadius * cos(CAngleXY) * sin(CAngleXZ);

	gViewMatrix = glm::lookAt(cameraCoords,	// eye
		glm::vec3(0.0, 0.0, 0.0),	// center
		glm::vec3(0.0, cameraUp, 0.0));	// up
}


void reset() {
	// Clear rotation
	cameraCoords = glm::vec3(10.0, 10.0, 10.0f);
	cameraUp = 1.0;
	// Init camera angle & radius for later rotation
	CRadius = sqrt(3 * (10 * 10));
	CAngleXZ = atan(1);
	CAngleXY = atan(1 / sqrt(2));

	gViewMatrix = glm::lookAt(cameraCoords,	// eye
		glm::vec3(0.0, 0.0, 0.0),	// center
		glm::vec3(0.0, cameraUp, 0.0));	// up

	// Activate face model
	//showFace = false;
	// Reset control points on X-Y plane
	//createObjects();
	float Pos[4] = { 0.0, 0.0, 0.0, 1.0 };
	float tempColor[4] = { 0.0, 1.0, 0.0, 1.0 };
	float tempNormal[3] = { 0.0, 0.0, 1.0 };
	for (int i = 0; i < 441; i++)
	{
		Pos[0] = ((-0.5)*(i / 21 - 10) + meshOffset[0] - 0.5);
		Pos[1] = (0.5*(i % 21 - 10) + 5.0f + meshOffset[1] + 0.5);
		CtlPoints[i].SetPosition(Pos);
		CtlPoints[i].SetColor(tempColor);
		// Set colors for picking
		CtlPoints[i].SetNormal(tempNormal);
	}
	VertexBufferSize[3] = sizeof(CtlPoints);
	createVAOs(CtlPoints, NULL, 3);
	updateControlMesh(1);
	hidesverts = true;
	smooth = false;
	//for (int i = 0; i < 441; i++)
	//{
	//	CtlPoints[i].Position[2] = 0.0;
	//}
	//VertexBufferSize[3] = sizeof(CtlPoints);
	//createVAOs(CtlPoints, NULL, 3);
	//updateControlMesh(1);
}


void renderScene(void)
{
	//ATTN: DRAW YOUR SCENE HERE. MODIFY/ADAPT WHERE NECESSARY!


	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.2f, 0.0f);
	// Re-clear the screen for real rendering
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(programID);
	{
		glPointSize(7.0);
		glm::vec3 lightPos = glm::vec3(0, 5, 5);
		glm::mat4x4 ModelMatrix = glm::mat4(1.0);
		glUniform3f(LightID, lightPos.x, lightPos.y, lightPos.z);
		glUniformMatrix4fv(ViewMatrixID, 1, GL_FALSE, &gViewMatrix[0][0]);
		glUniformMatrix4fv(ProjMatrixID, 1, GL_FALSE, &gProjectionMatrix[0][0]);
		glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);

		glUniform1i(useLightModelID, false);
		glUniform1i(useTextureID, false);
		glBindVertexArray(VertexArrayId[0]);	// draw CoordAxes
		glDrawArrays(GL_LINES, 0, 6);

		if (Ani_angry)
		{
			if (time_angry > 1.5)
				Ani_angry = false;
			else
			{
				CtlPoints[241].Position[1] -= 0.002f;
				CtlPoints[241].Position[0] += 0.001f;
				CtlPoints[240].Position[1] -= 0.002f;
				CtlPoints[240].Position[0] += 0.001f;
				CtlPoints[199].Position[1] -= 0.001f;
				CtlPoints[199].Position[0] -= 0.001f;
				VertexBufferSize[3] = sizeof(CtlPoints);	// ATTN: this needs to be done for each hand-made object with the ObjectID (subscript)
				createVAOs(CtlPoints, NULL, 3);
				updateControlMesh(0);
				//printf("angry!!");
			}
		}

		if (Ani_smile)
		{
			if (time_smile > 1.5)
				Ani_smile = false;
			CtlPoints[235].Position[1] += 0.002f;
			CtlPoints[235].Position[0] -= 0.001f;
			CtlPoints[193].Position[1] += 0.002f;
			CtlPoints[193].Position[0] -= 0.001f;
			CtlPoints[214].Position[1] += 0.001f;
			CtlPoints[214].Position[0] += 0.001f;
			VertexBufferSize[3] = sizeof(CtlPoints);	// ATTN: this needs to be done for each hand-made object with the ObjectID (subscript)
			createVAOs(CtlPoints, NULL, 3);
			updateControlMesh(0);
		}

		if (showFace) {
			glUniform1i(useLightModelID, true);
			glBindVertexArray(VertexArrayId[1]);
			glDrawElements(GL_TRIANGLES, VertexBufferSize[1], GL_UNSIGNED_INT, NULL);
		}
		// draw grid
		glUniform1i(useLightModelID, false);
		glBindVertexArray(VertexArrayId[2]);
		glDrawArrays(GL_LINES, 0, 44);




		if (showControlMesh) {
			// Draw Mesh Verticies
			if (!hide)
			{
				glBindVertexArray(VertexArrayId[3]);
				glDrawArrays(GL_POINTS, 0, 441);

				// Draw lines of Mesh
				glBindVertexArray(VertexArrayId[4]);
				glDrawArrays(GL_LINES, 0, 1680);
			}

			if (smooth)
			{
				if (!hidesverts)
				{
					glBindVertexArray(VertexArrayId[6]);
					glDrawArrays(GL_POINTS, 0, 3249);
					// Draw Triangles for Mesh Using Texture
				}
				glUniform1i(useTextureID, true);
				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, texture);
				glUniform1i(TextureID, 0);
				glBindVertexArray(VertexArrayId[7]);
				glDrawArrays(GL_TRIANGLES, 0, 19494);
			}
			else
			{
				// Draw Triangles for Mesh Using Texture
				glUniform1i(useTextureID, true);
				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, texture);
				glUniform1i(TextureID, 0);
				glBindVertexArray(VertexArrayId[5]);
				glDrawArrays(GL_TRIANGLES, 0, 2400);
			}

		}
		glBindVertexArray(0);

	}
	glUseProgram(0);

	glfwSwapBuffers(window);
	glfwPollEvents();
}

void pickObject(void)
{
	// Clear the screen in white
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(pickingProgramID);
	{
		glm::mat4 ModelMatrix = glm::mat4(1.0); // TranslationMatrix * RotationMatrix;
		glm::mat4 MVP = gProjectionMatrix * gViewMatrix * ModelMatrix;

		// Send our transformation to the currently bound shader, in the "MVP" uniform
		glUniformMatrix4fv(PickingMatrixID, 1, GL_FALSE, &MVP[0][0]);

		// ATTN: DRAW YOUR PICKING SCENE HERE. REMEMBER TO SEND IN A DIFFERENT PICKING COLOR FOR EACH OBJECT BEFOREHAND
		glUniform1fv(pickingColorArrayID, 441, pickingColor);	// here we pass in the picking marker array

		// Draw the points
		glBindVertexArray(VertexArrayId[3]);
		glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[3]);
		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(CtlPoints), CtlPoints);	// update buffer data
		glDrawArrays(GL_POINTS, 0, 441);
		glBindVertexArray(0);

	}
	glUseProgram(0);
	// Wait until all the pending drawing commands are really done.
	// Ultra-mega-over slow !
	// There are usually a long time between glDrawElements() and
	// all the fragments completely rasterized.
	glFlush();
	glFinish();

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	// Read the pixel at the center of the screen.
	// You can also use glfwGetMousePos().
	// Ultra-mega-over slow too, even for 1 pixel,
	// because the framebuffer is on the GPU.
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);
	unsigned char data[4];
	// Convert the color back to an integer ID
	glReadPixels(xpos, window_height - ypos, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, data); // OpenGL renders with (0,0) on bottom, mouse reports with (0,0) on top
	glReadPixels(xpos, window_height - ypos, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &zpos);
	// Convert the color back to an integer ID
	gPickedIndex = int(data[0]) + int(data[1]);
	printf("%d %d %d\n", data[0], data[1], gPickedIndex);
	if (gPickedIndex >= 0 && gPickedIndex <= 440)
	{
		std::ostringstream oss;
		oss << "point " << gPickedIndex;
		gMessage = oss.str();
		isPick = true;
		oriCoor.x = CtlPoints[gPickedIndex].Position[0];
		oriCoor.y = CtlPoints[gPickedIndex].Position[1];
		oriCoor.z = CtlPoints[gPickedIndex].Position[2];
	}
	else
	{
		gMessage = "background";
	}
}

void moveVertex(void)
{
	glm::mat4 ModelMatrix = glm::mat4(1.0);
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	glm::vec4 vp = glm::vec4(viewport[0], viewport[1], viewport[2], viewport[3]);
	double xpos = 0, ypos = 0;
	glfwGetCursorPos(window, &xpos, &ypos);
	glm::vec3 oriCoordinate = glm::vec3(xpos, window_height - ypos, zpos);
	glm::vec3 newCoordinate = glm::unProject(oriCoordinate, gViewMatrix * ModelMatrix, gProjectionMatrix, vp);

	if (!is_zDirection)
	{
		curCoor.x = newCoordinate.x;
		curCoor.y = newCoordinate.y;
		curCoor.z = oriCoor.z;
	}
	else
	{
		curCoor.x = oriCoor.x;
		curCoor.y = oriCoor.y;
		curCoor.z = oriCoor.z + (oriCoor.x - newCoordinate.x);
	}
	float coord[3] = { curCoor.x, curCoor.y, curCoor.z };
	CtlPoints[gPickedIndex].SetPosition(coord);
	VertexBufferSize[3] = sizeof(CtlPoints);	// ATTN: this needs to be done for each hand-made object with the ObjectID (subscript)
	createVAOs(CtlPoints, NULL, 3);

	updateControlMesh(0);
	// CONTROL MESH

}

void ray_casting() {
	int triIndex = 0;
	float rayDirection[3] = { 0.0, 0.0, 1.0 };
	float bc[3];
	int change = 0;
	float tempPos[4];
	Vertex* Verts;
	GLuint* Idcs;
	loadObject("HCY.obj", glm::vec4(1.0), Verts, Idcs, 1);
	for (int ctlVertCount = 0; ctlVertCount<441; ctlVertCount++) {
		for (triIndex = 0; triIndex < NumIndices[1] / 3; triIndex++) {

			ray_cast(Verts[Idcs[triIndex * 3]].Position, Verts[Idcs[triIndex * 3 + 1]].Position, Verts[Idcs[triIndex * 3 + 2]].Position, CtlPoints[ctlVertCount].Position, rayDirection, bc);

			if (bc[0] >= 0 && bc[1] >= 0 && bc[2] >= 0) {
				float tempZ = bc[0] * Verts[Idcs[triIndex * 3]].Position[2] +
					bc[1] * Verts[Idcs[triIndex * 3 + 1]].Position[2] +
					bc[2] * Verts[Idcs[triIndex * 3 + 2]].Position[2];
				if (tempZ < CtlPoints[ctlVertCount].Position[2])
					continue;
				for (int i = 0; i<2; i++) {
					tempPos[i] = bc[0] * Verts[Idcs[triIndex * 3]].Position[i] +
						bc[1] * Verts[Idcs[triIndex * 3 + 1]].Position[i] +
						bc[2] * Verts[Idcs[triIndex * 3 + 2]].Position[i];
				}
				tempPos[2] = tempZ + 1;
				tempPos[3] = 1.0;
				CtlPoints[ctlVertCount].SetPosition(tempPos);
				change = 1;
			}
		}
		if (change == 0)
			CtlPoints[ctlVertCount].Position[2] = 0.0;
		else
			change = 0;
	}
	glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[3]);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(CtlPoints), CtlPoints);
	updateControlMesh(0);
}

void 
ray_casting2() {
	for (int i = 0; i < 441; i++)
	{
			//if (CtlPoints[i].Position[1]>7)
			//{
			//	CtlPoints[i].Position[2] = -1;
			//	CtlPoints[i].Position[1] = 5;
			//}
			//if (CtlPoints[i].Position[0] < -3 && CtlPoints[i].Position[1] <= 7)
			//{
			//	float a = CtlPoints[i].Position[0] + 3;
			//	CtlPoints[i].Position[0] = -3;
			//	CtlPoints[i].Position[2] = a;
			//}
			if (CtlPoints[i].Position[2]<0)
			{
				//float a = CtlPoints[i].Position[0] - 3.5;
				CtlPoints[i].Position[0] = 0;
				CtlPoints[i].Position[1] = 3;
				CtlPoints[i].Position[2] = 0;
			}
	}
	VertexBufferSize[3] = sizeof(CtlPoints);	// ATTN: this needs to be done for each hand-made object with the ObjectID (subscript)
	createVAOs(CtlPoints, NULL, 3);
	updateControlMesh(0);
}

int initWindow(void)
{
	// Initialise GLFW
	if (!glfwInit()) {
		fprintf(stderr, "Failed to initialize GLFW\n");
		return -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

	// Open a window and create its OpenGL context
	window = glfwCreateWindow(window_width, window_height, "Hu,Chenyang(UFID:94699563)", NULL, NULL);
	if (window == NULL) {
		fprintf(stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n");
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		return -1;
	}

	// Initialize the GUI
	//TwInit(TW_OPENGL_CORE, NULL);
	//TwWindowSize(window_width, window_height);
	//TwBar * GUI = TwNewBar("Picking");
	//TwSetParam(GUI, NULL, "refresh", TW_PARAM_CSTRING, 1, "0.1");
	//TwAddVarRW(GUI, "Last picked object", TW_TYPE_STDSTRING, &gMessage, NULL);

	// Set up inputs
	glfwSetCursorPos(window, window_width / 2, window_height / 2);
	glfwSetKeyCallback(window, keyCallback);
	glfwSetMouseButtonCallback(window, mouseCallback);

	return 0;
}

void initOpenGL(void)
{

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);
	// Cull triangles which normal is not towards the camera
	glEnable(GL_CULL_FACE);

	// Projection matrix : 458 Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
	gProjectionMatrix = glm::perspective(45.0f, 1.0f, 0.1f, 100.0f);
	// Or, for an ortho camera :
	//gProjectionMatrix = glm::ortho(-4.0f, 4.0f, -3.0f, 3.0f, 0.0f, 100.0f); // In world coordinates

	// Camera matrix
	cameraCoords = glm::vec3(10.0, 10.0, 10.0f);
	cameraUp = 1.0;
	// Init camera angle & radius for later rotation
	CRadius = sqrt(3 * (10 * 10));
	CAngleXZ = atan(1);
	CAngleXY = atan(1 / sqrt(2));

	gViewMatrix = glm::lookAt(cameraCoords,	// eye
		glm::vec3(0.0, 0.0, 0.0),	// center
		glm::vec3(0.0, cameraUp, 0.0));	// up


	// Create and compile our GLSL program from the shaders
	programID = LoadShaders("StandardShading.vertexshader", "StandardShading.fragmentshader");
	pickingProgramID = LoadShaders("Picking.vertexshader", "Picking.fragmentshader");

	// Get a handle for our "MVP" uniform
	MatrixID = glGetUniformLocation(programID, "MVP");
	ModelMatrixID = glGetUniformLocation(programID, "M");
	ViewMatrixID = glGetUniformLocation(programID, "V");
	ProjMatrixID = glGetUniformLocation(programID, "P");

	PickingMatrixID = glGetUniformLocation(pickingProgramID, "MVP");
	// Get a handle for our "pickingColorID" uniform
	//pickingColorID = glGetUniformLocation(pickingProgramID, "PickingColor");

	//Get a handle for our "pickingColorArrayID" uniform
	pickingColorArrayID = glGetUniformLocation(pickingProgramID, "PickingColorArray");
	pickingColorID = glGetUniformLocation(pickingProgramID, "PickingColor");

	// Get a handle for our "LightPosition" uniform
	LightID = glGetUniformLocation(programID, "LightPosition_worldspace");

	// Get a handle for our "useLightModel" uniform
	useLightModelID = glGetUniformLocation(programID, "useLightModel");
	useTextureID = glGetUniformLocation(programID, "useTexture");
	TextureID = glGetUniformLocation(programID, "textureSampler");


	for (int i = 0; i <= 20; i++)
	{
		for (int j = 0; j <= 20; j++)
		{
			int index = i * 21 + j;
			pickingColor[index] = float(index) / 255.0f;
			//printf("%d %f\n", index, pickingColor[index]);
		}
	}
	createObjects();
}

void createVAOs(Vertex Vertices[], unsigned int Indices[], int ObjectId) {

	GLenum ErrorCheckValue = glGetError();
	const size_t VertexSize = sizeof(Vertices[0]);
	const size_t RgbOffset = sizeof(Vertices[0].Position);
	const size_t Normaloffset = sizeof(Vertices[0].Color) + RgbOffset;
	const size_t UVoffset = sizeof(Vertices[0].Normal) + Normaloffset;

	// Create Vertex Array Object
	glGenVertexArrays(1, &VertexArrayId[ObjectId]);	//
	glBindVertexArray(VertexArrayId[ObjectId]);		//

	// Create Buffer for vertex data
	glGenBuffers(1, &VertexBufferId[ObjectId]);
	glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[ObjectId]);
	glBufferData(GL_ARRAY_BUFFER, VertexBufferSize[ObjectId], Vertices, GL_STATIC_DRAW);

	// Create Buffer for indices
	if (Indices != NULL) {
		glGenBuffers(1, &IndexBufferId[ObjectId]);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IndexBufferId[ObjectId]);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, IndexBufferSize[ObjectId], Indices, GL_STATIC_DRAW);
	}

	// Assign vertex attributes
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, VertexSize, 0);
	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, VertexSize, (GLvoid*)RgbOffset);
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, VertexSize, (GLvoid*)Normaloffset);
	glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, VertexSize, (GLvoid*)UVoffset);


	glEnableVertexAttribArray(0);	// position
	glEnableVertexAttribArray(1);	// color
	glEnableVertexAttribArray(2);	// normal
	glEnableVertexAttribArray(3);	// UVs


	// Disable our Vertex Buffer Object
	glBindVertexArray(0);

	ErrorCheckValue = glGetError();
	if (ErrorCheckValue != GL_NO_ERROR)
	{
		fprintf(
			stderr,
			"ERROR: Could not create a VBO: %s \n",
			gluErrorString(ErrorCheckValue)
			);
	}
}

void SMVerts()
{
	int index = 0;
	for (int i = 1; i < 20; i++)
	{
		for (int j = 1; j < 20; j++)
		{
			Vertex s11 = CtlPoints[i * 21 + j];
			Vertex s12 = CtlPoints[(i + 1) * 21 + j];
			Vertex s21 = CtlPoints[i * 21 + j + 1];
			Vertex s10 = CtlPoints[(i - 1) * 21 + j];
			Vertex s01 = CtlPoints[i * 21 + j - 1];
			Vertex s00 = CtlPoints[(i - 1) * 21 + j - 1];
			Vertex s20 = CtlPoints[(i - 1) * 21 + j + 1];
			Vertex s02 = CtlPoints[(i + 1) * 21 + j - 1];
			Vertex s22 = CtlPoints[(i + 1) * 21 + j + 1];
			float color[3] = { 1.0f, 0.0f, 0.0f };
			SMeshPoints[171 * (i - 1) + 3 * (j - 1)].SetColor(color);
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 1].SetColor(color);
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 2].SetColor(color);
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 57].SetColor(color);
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 58].SetColor(color);
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 59].SetColor(color);
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 114].SetColor(color);
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 115].SetColor(color);
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 116].SetColor(color);
			SMeshPoints[171 * (i - 1) + 3 * (j - 1)].Position[0] = (16 * s11.Position[0] + 4 * (s21.Position[0] + s12.Position[0] + s01.Position[0] + s10.Position[0]) + s22.Position[0] + s02.Position[0] + s00.Position[0] + s20.Position[0]) / 36;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1)].Position[1] = (16 * s11.Position[1] + 4 * (s21.Position[1] + s12.Position[1] + s01.Position[1] + s10.Position[1]) + s22.Position[1] + s02.Position[1] + s00.Position[1] + s20.Position[1]) / 36;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1)].Position[2] = (16 * s11.Position[2] + 4 * (s21.Position[2] + s12.Position[2] + s01.Position[2] + s10.Position[2]) + s22.Position[2] + s02.Position[2] + s00.Position[2] + s20.Position[2]) / 36;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 1].Position[0] = (8 * s11.Position[0] + 2 * (s10.Position[0] + s12.Position[0]) + 4 * s21.Position[0] + s22.Position[0] + s20.Position[0]) / 18;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 1].Position[1] = (8 * s11.Position[1] + 2 * (s10.Position[1] + s12.Position[1]) + 4 * s21.Position[1] + s22.Position[1] + s20.Position[1]) / 18;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 1].Position[2] = (8 * s11.Position[2] + 2 * (s10.Position[2] + s12.Position[2]) + 4 * s21.Position[2] + s22.Position[2] + s20.Position[2]) / 18;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 2].Position[0] = (8 * s21.Position[0] + 2 * (s22.Position[0] + s20.Position[0]) + 4 * s11.Position[0] + s12.Position[0] + s10.Position[0]) / 18;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 2].Position[1] = (8 * s21.Position[1] + 2 * (s22.Position[1] + s20.Position[1]) + 4 * s11.Position[1] + s12.Position[1] + s10.Position[1]) / 18;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 2].Position[2] = (8 * s21.Position[2] + 2 * (s22.Position[2] + s20.Position[2]) + 4 * s11.Position[2] + s12.Position[2] + s10.Position[2]) / 18;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 57].Position[0] = (8 * s11.Position[0] + 2 * (s01.Position[0] + s21.Position[0]) + 4 * s12.Position[0] + s02.Position[0] + s22.Position[0]) / 18;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 57].Position[1] = (8 * s11.Position[1] + 2 * (s01.Position[1] + s21.Position[1]) + 4 * s12.Position[1] + s02.Position[1] + s22.Position[1]) / 18;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 57].Position[2] = (8 * s11.Position[2] + 2 * (s01.Position[2] + s21.Position[2]) + 4 * s12.Position[2] + s02.Position[2] + s22.Position[2]) / 18;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 58].Position[0] = (4 * s11.Position[0] + 2 * (s21.Position[0] + s12.Position[0]) + s22.Position[0]) / 9;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 58].Position[1] = (4 * s11.Position[1] + 2 * (s21.Position[1] + s12.Position[1]) + s22.Position[1]) / 9;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 58].Position[2] = (4 * s11.Position[2] + 2 * (s21.Position[2] + s12.Position[2]) + s22.Position[2]) / 9;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 59].Position[0] = (4 * s21.Position[0] + 2 * (s11.Position[0] + s22.Position[0]) + s12.Position[0]) / 9;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 59].Position[1] = (4 * s21.Position[1] + 2 * (s11.Position[1] + s22.Position[1]) + s12.Position[1]) / 9;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 59].Position[2] = (4 * s21.Position[2] + 2 * (s11.Position[2] + s22.Position[2]) + s12.Position[2]) / 9;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 114].Position[0] = (8 * s12.Position[0] + 2 * (s02.Position[0] + s22.Position[0]) + 4 * s11.Position[0] + s01.Position[0] + s21.Position[0]) / 18;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 114].Position[1] = (8 * s12.Position[1] + 2 * (s02.Position[1] + s22.Position[1]) + 4 * s11.Position[1] + s01.Position[1] + s21.Position[1]) / 18;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 114].Position[2] = (8 * s12.Position[2] + 2 * (s02.Position[2] + s22.Position[2]) + 4 * s11.Position[2] + s01.Position[2] + s21.Position[2]) / 18;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 115].Position[0] = (4 * s12.Position[0] + 2 * (s11.Position[0] + s22.Position[0]) + s21.Position[0]) / 9;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 115].Position[1] = (4 * s12.Position[1] + 2 * (s11.Position[1] + s22.Position[1]) + s21.Position[1]) / 9;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 115].Position[2] = (4 * s12.Position[2] + 2 * (s11.Position[2] + s22.Position[2]) + s21.Position[2]) / 9;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 116].Position[0] = (4 * s22.Position[0] + 2 * (s12.Position[0] + s21.Position[0]) + s11.Position[0]) / 9;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 116].Position[1] = (4 * s22.Position[1] + 2 * (s12.Position[1] + s21.Position[1]) + s11.Position[1]) / 9;
			SMeshPoints[171 * (i - 1) + 3 * (j - 1) + 116].Position[2] = (4 * s22.Position[2] + 2 * (s12.Position[2] + s21.Position[2]) + s11.Position[2]) / 9;
		}
	}

	VertexBufferSize[6] = sizeof(SMeshPoints);
	createVAOs(SMeshPoints, NULL, 6);
}

void cleanup(void)
{
	// Cleanup VBO and shader
	for (int i = 0; i < NumObjects; i++) {
		glDeleteBuffers(1, &VertexBufferId[i]);
		glDeleteBuffers(1, &IndexBufferId[i]);
		glDeleteVertexArrays(1, &VertexArrayId[i]);
	}
	glDeleteProgram(programID);
	glDeleteProgram(pickingProgramID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();
}

static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	// ATTN: MODIFY AS APPROPRIATE
	if (action == GLFW_RELEASE) {
		inputMode = 0;
		direction = 0;
	}
	if (action == GLFW_PRESS) {
		switch (key)
		{
			// Direction Keys
		case GLFW_KEY_LEFT:
		case GLFW_KEY_RIGHT:
		case GLFW_KEY_UP:
		case GLFW_KEY_DOWN:
			direction = key;
			break;
		case GLFW_KEY_Q:
			ray_casting2();

			break;
		case GLFW_KEY_H:
			if (hide == 0)
				hide = 1;
			else
				hide = 0;
			break;
		case GLFW_KEY_F:
			if (showFace)
				showFace = false;
			else
				showFace = true;
			break;
		case GLFW_KEY_C:
		{
			if (showControlMesh)
			   showControlMesh = false;
			else
			   showControlMesh = true;
			SMVerts();
		}break;
		case GLFW_KEY_J:
		{
			SMVerts();
			showSMFace(true);
			if (smooth)
			   smooth = false;
			else
			   smooth = true;
		}break;
		case GLFW_KEY_K:
		{
			if (hidesverts)
			   hidesverts = false;
			else
			   hidesverts = true;
		}break;
		case GLFW_KEY_A:
		{
			ray_casting();
			SMVerts();
			showSMFace(true);
		}break;
		case GLFW_KEY_P:
		{
			if (Ani_smile)
			   Ani_smile = false;
			else
			   Ani_smile = true;
		}break;
		case GLFW_KEY_O:
		{
			if (Ani_angry)
			   Ani_angry = false;
			else
			   Ani_angry = true;
		}break;
		case GLFW_KEY_S:
		{
			fileexport();
		}break;
		case GLFW_KEY_L:
		{
			fileload();
			SMVerts();
			showSMFace(true);
		}break;
		case GLFW_KEY_R:
		{
			reset();
		}break;

		default:
			break;
		}
	}
}

static void mouseCallback(GLFWwindow* window, int button, int action, int mods)
{
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		{
			enableDrag = true;
			pickObject();
		}
	}

	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
		enableDrag = false;
		pickedIndex = -1;
		lastXPos = lastYPos = -1;
	}
}

int main(void)
{
	// initialize window
	int errorCode = initWindow();
	if (errorCode != 0)
		return errorCode;

	// initialize OpenGL pipeline
	initOpenGL();

	do {

		if (direction) moveCamera(glm::pi<float>() / 72.0);

		if (enableDrag) {
			moveVertex();
			updateControlMesh(0);
			SMVerts();
		}
		opepoints();
		// DRAWING POINTS
		renderScene();
		if (Ani_smile)
		{
			time_smile += 0.01;
		}
		if (Ani_angry)
		{
			time_angry += 0.01;
		}
	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
	glfwWindowShouldClose(window) == 0);

	cleanup();

	return 0;
}
