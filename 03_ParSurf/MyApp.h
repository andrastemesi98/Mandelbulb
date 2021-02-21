#pragma once

// GLEW
#include <GL/glew.h>

// SDL
#include <SDL.h>
#include <SDL_opengl.h>

// GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform2.hpp>

#include "gCamera.h"

class CMyApp
{
public:
	CMyApp(void);
	~CMyApp(void);

	bool Init();
	void Clean();

	void Update();
	void Render();

	void KeyboardDown(SDL_KeyboardEvent&);
	void KeyboardUp(SDL_KeyboardEvent&);
	void MouseMove(SDL_MouseMotionEvent&);
	void MouseDown(SDL_MouseButtonEvent&);
	void MouseUp(SDL_MouseButtonEvent&);
	void MouseWheel(SDL_MouseWheelEvent&);
	void Resize(int, int);
protected:

	// shaderekhez szükséges változók
	GLuint m_programID; // shaderek programja

	// OpenGL-es dolgok
	GLuint m_vaoID; // vertex array object erõforrás azonosító
	GLuint m_vboID; // vertex buffer object erõforrás azonosító
	GLuint m_ibID;  // index buffer object erõforrás azonosító

	// transzformációs mátrixok
	glm::mat4 m_matWorld;
	glm::mat4 m_matView;
	glm::mat4 m_matProj;

	// mátrixok helye a shaderekben
	GLuint	m_loc_mvp; // a három mátrixunk szorzatát adjuk át a hatékonyság érdekében
	GLuint m_loc_win;
	GLuint m_loc_eye;
	GLuint m_loc_at;
	GLuint m_loc_up;
	GLuint m_loc_time;

	gCamera	m_camera;

	struct Vertex
	{
		glm::vec2 p;
		glm::vec2 c;
	};

	float width = 640;
	float height = 480;


};

